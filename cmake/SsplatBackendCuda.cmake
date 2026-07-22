# SSPLAT_BACKEND=cuda (default): the CUDA kernels plus, when Torch and Python
# are available, the `csrc` Python extension module.
#
# Defines, for the shared app targets in SsplatApps.cmake:
#   SSPLAT_WITH_TORCH   ON when the Python extension is being built
#   SSPLAT_APP_LIBS     the engine library the apps link
#   SPLAT_CXX_FLAGS     host flags the app targets reuse

enable_language(CUDA)

# ---------------------------------------------------------------------------
# Find Torch (optional)
# ---------------------------------------------------------------------------
set(SSPLAT_WITH_TORCH OFF)
if(NOT SSPLAT_NO_TORCH)
    find_package(Python3 QUIET COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
        execute_process(
            COMMAND ${Python3_EXECUTABLE} -c "import torch; print(torch.utils.cmake_prefix_path)"
            OUTPUT_VARIABLE PYTORCH_CMAKE_PREFIX_PATH
            RESULT_VARIABLE SSPLAT_TORCH_PROBE_RESULT
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        message(STATUS "PyTorch CMake prefix path: ${PYTORCH_CMAKE_PREFIX_PATH}")

        execute_process(
            COMMAND ${Python3_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_paths()['include'])"
            OUTPUT_VARIABLE PYTHON_ADDITION_INCLUDE_PATH
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        message(STATUS "Python additional include path: ${PYTHON_ADDITION_INCLUDE_PATH}")

        if(SSPLAT_TORCH_PROBE_RESULT EQUAL 0)
            set(CMAKE_PREFIX_PATH ${PYTORCH_CMAKE_PREFIX_PATH})
            find_package(Torch QUIET)
        endif()
    endif()

    if(Torch_FOUND)
        message(STATUS "Found Torch: ${TORCH_VERSION}")
        find_library(TORCH_PYTHON_LIBRARY torch_python PATHS "${TORCH_INSTALL_PREFIX}/lib")
        set(SSPLAT_WITH_TORCH ON)
    endif()
endif()

if(NOT SSPLAT_WITH_TORCH)
    message(STATUS "Torch/Python3 not available (or SSPLAT_NO_TORCH=ON) -- "
        "skipping the Python extension; building the standalone ssplat-train CLI only.")
    set(SSPLAT_BUILD_CLI ON)
endif()

find_package(CUDAToolkit REQUIRED)

# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------
ssplat_collect_sources(SPLAT_SOURCES)

if(NOT SSPLAT_WITH_TORCH)
    # bindings/ext.cpp (the pybind11 module) is the only Torch-dependent source.
    list(REMOVE_ITEM SPLAT_SOURCES ${SSPLAT_SRC}/bindings/ext.cpp)
endif()

# ---------------------------------------------------------------------------
# Detect CUDA architectures
# ---------------------------------------------------------------------------
if(NOT DEFINED TORCH_CUDA_ARCH_LIST)
    if(SSPLAT_WITH_TORCH)
        execute_process(
            COMMAND ${Python3_EXECUTABLE} -c "import torch; print(' '.join(f'{a[0]}{a[1]}' for a in [torch.cuda.get_device_capability(i) for i in range(torch.cuda.device_count())]))"
            OUTPUT_VARIABLE TORCH_CUDA_ARCH_LIST
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()

    if(NOT TORCH_CUDA_ARCH_LIST)
        # No Torch (or Torch saw no device): ask the driver directly.
        execute_process(
            COMMAND nvidia-smi --query-gpu=compute_cap --format=csv,noheader
            OUTPUT_VARIABLE SSPLAT_SMI_CAPS
            RESULT_VARIABLE SSPLAT_SMI_RESULT
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(SSPLAT_SMI_RESULT EQUAL 0)
            string(REPLACE "\r" "" SSPLAT_SMI_CAPS "${SSPLAT_SMI_CAPS}")
            string(REPLACE "\n" " " TORCH_CUDA_ARCH_LIST "${SSPLAT_SMI_CAPS}")
        endif()
    endif()

    if(NOT TORCH_CUDA_ARCH_LIST)
        message(FATAL_ERROR "Could not detect a CUDA GPU. "
            "Pass -DTORCH_CUDA_ARCH_LIST=\"8.6\" (or similar) to set the target architecture.")
    endif()
endif()

# normalize: "8.6" / "8.6 7.5" / "86;75" all become a deduped list like "86;75"
string(REGEX REPLACE "\\." "" TORCH_CUDA_ARCH_LIST "${TORCH_CUDA_ARCH_LIST}")
string(REPLACE " " ";" TORCH_CUDA_ARCH_LIST "${TORCH_CUDA_ARCH_LIST}")
list(REMOVE_DUPLICATES TORCH_CUDA_ARCH_LIST)

message(STATUS "CUDA architecture(s): ${TORCH_CUDA_ARCH_LIST}")

set(CMAKE_CUDA_ARCHITECTURES "")
foreach(arch ${TORCH_CUDA_ARCH_LIST})
    list(APPEND CMAKE_CUDA_ARCHITECTURES "${arch}")
endforeach()

# ---------------------------------------------------------------------------
# Compiler flags
#
# These intentionally mirror, but do not share, setup.py's flags: the pip build
# targets a released wheel (its own WITH_SYMBOLS / LINE_INFO env gates) while
# this one targets local development. Only the *source list* is shared, via
# cmake/sources.txt.
# ---------------------------------------------------------------------------
if(MSVC)
    set(SPLAT_CXX_FLAGS "/O2")
else()
    set(SPLAT_CXX_FLAGS "-O3")
    list(APPEND SPLAT_CXX_FLAGS "-Wno-sign-compare")
endif()

if(MSVC)
    add_compile_options(
        $<$<COMPILE_LANGUAGE:CXX>:/Zc:preprocessor>
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/Zc:preprocessor>
    )
endif()

# Detect OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND AND NOT APPLE)
    message(STATUS "Compiling with OpenMP")
    list(APPEND SPLAT_CXX_FLAGS "-DAT_PARALLEL_OPENMP")
    if(SSPLAT_WITH_TORCH)
        # Compile flag only: the OpenMP runtime comes from torch's bundled
        # libgomp. Linking the system one too would shadow it (CMake warns
        # about the unsafe rpath). The no-torch build instead links the
        # OpenMP::OpenMP_CXX target, whose flags are language-guarded so
        # nvcc never sees them.
        list(APPEND SPLAT_CXX_FLAGS ${OpenMP_CXX_FLAGS})
    endif()
else()
    message(STATUS "Compiling without OpenMP...")
endif()

# macOS ARM64 special case
if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
    message(STATUS "Enabling macOS ARM64 flags")
    add_compile_options("-arch" "arm64")
    add_link_options("-arch" "arm64")
endif()

# nvcc flags
set(SPLAT_NVCC_FLAGS
    "-O3" "--use_fast_math"
    "--expt-relaxed-constexpr"
    "-Xcudafe=--diag_suppress=20012"
    "-Xcudafe=--diag_suppress=550"
    "--threads" "0"
    "-Xptxas" "--warn-on-double-precision-use"
)

# cubin/PTX line info (adds embedded source; needed for profiling but bloats
# the fatbin). Off unless SSPLAT_DEBUG_SYMBOLS.
if(SSPLAT_DEBUG_SYMBOLS)
    list(APPEND SPLAT_NVCC_FLAGS "-lineinfo" "--generate-line-info" "--source-in-ptx")
endif()

# Add gencode flags
foreach(arch ${TORCH_CUDA_ARCH_LIST})
    list(APPEND SPLAT_NVCC_FLAGS "-gencode" "arch=compute_${arch},code=sm_${arch}")
endforeach()

# Host side optimizations
if(NOT WIN32)
    # list(APPEND SPLAT_NVCC_FLAGS "-Xcompiler=-O3,-march=native")
    # Host debug symbols (nvcc -Xcompiler=-g). Off unless SSPLAT_DEBUG_SYMBOLS.
    if(SSPLAT_DEBUG_SYMBOLS)
        list(APPEND SPLAT_NVCC_FLAGS "-Xcompiler=-g")
    endif()
endif()

# ---------------------------------------------------------------------------
# The engine library / Python extension module
# ---------------------------------------------------------------------------
if(SSPLAT_WITH_TORCH)
    add_library(csrc SHARED ${SPLAT_SOURCES})
else()
    # The engine API has no dllexport annotations (Linux .so exports
    # everything by default); a static lib sidesteps that on Windows and
    # makes ssplat-train a self-contained executable.
    add_library(csrc STATIC ${SPLAT_SOURCES})
endif()

target_include_directories(csrc PRIVATE
    ${SSPLAT_SRC}
    ${CUDAToolkit_INCLUDE_DIRS}
)

if(WIN32)
    target_compile_definitions(csrc PRIVATE spirulae_splat_EXPORTS)
endif()

if(SSPLAT_WITH_TORCH)
    target_include_directories(csrc PRIVATE
        ${Python3_INCLUDE_DIRS}
        ${PYTHON_ADDITION_INCLUDE_PATH}
    )

    target_compile_definitions(csrc PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:TORCH_EXTENSION_NAME=csrc>
    )

    target_link_libraries(csrc
        ${TORCH_LIBRARIES}
        ${TORCH_PYTHON_LIBRARY}
        ${Python3_LIBRARIES}
    )
else()
    # Torch normally drags these in for us
    find_package(Threads REQUIRED)
    target_link_libraries(csrc CUDA::cudart Threads::Threads)
    if(OpenMP_CXX_FOUND AND NOT APPLE)
        target_link_libraries(csrc OpenMP::OpenMP_CXX)
    endif()
endif()

target_compile_options(csrc PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:${SPLAT_CXX_FLAGS}>
    $<$<COMPILE_LANGUAGE:CUDA>:${SPLAT_NVCC_FLAGS}>
)

# Optional: strip symbols (setup.py uses '-s' if WITH_SYMBOLS==False)
if(DEFINED WITH_SYMBOLS AND NOT WITH_SYMBOLS)
    if(NOT WIN32)
        target_link_options(csrc PRIVATE "-s")
    endif()
endif()

# Required by PyTorch C++ extensions
set_property(TARGET csrc PROPERTY CXX_STANDARD 17)
set_property(TARGET csrc PROPERTY CUDA_STANDARD 17)

set(SSPLAT_APP_LIBS csrc CUDA::cudart)

# ---------------------------------------------------------------------------
# Cross-backend parity tools -- CUDA side (dump reference outputs; the Vulkan
# build's twin executables compare against them). See backend/vulkan/README.md.
# The Vulkan branch builds its comparing twins unconditionally.
# ---------------------------------------------------------------------------
if(SSPLAT_BUILD_BACKEND_TESTS)
    if(SSPLAT_WITH_TORCH)
        # Same libpython + static-libstdc++ notes as ssplat-train.
        find_package(Python3 REQUIRED COMPONENTS Development.Embed)
    endif()
    file(GLOB SSPLAT_PARITY_TESTS CONFIGURE_DEPENDS
        ${SSPLAT_SRC}/backend/tests/*.cpp
        ${SSPLAT_SRC}/backend/tests/engine/*.cpp)
    foreach(test_src ${SSPLAT_PARITY_TESTS})
        get_filename_component(test_name ${test_src} NAME_WE)
        add_executable(${test_name} ${test_src})
        target_include_directories(${test_name} PRIVATE
            ${SSPLAT_SRC}
            ${CUDAToolkit_INCLUDE_DIRS})
        target_link_libraries(${test_name} PRIVATE csrc CUDA::cudart)
        if(SSPLAT_WITH_TORCH)
            target_link_libraries(${test_name} PRIVATE Python3::Python)
            if(NOT WIN32)
                target_link_options(${test_name} PRIVATE
                    -static-libgcc -static-libstdc++)
            endif()
        endif()
    endforeach()
endif()
