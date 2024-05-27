cmake_minimum_required(VERSION 3.12)  # You can adjust the minimum required version
set(CMAKE_CUDA_ARCHITECTURES 70 75 89)  # Ti 2080 uses 75. V100 uses 70. RTX 4090 uses 89.

project(spirulae_splat CXX CUDA)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CUDA_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# our spirulae_splat library
add_library(spirulae_splat forward.cu backward.cu helpers.cuh ch.cuh)
target_link_libraries(spirulae_splat PUBLIC cuda)
target_include_directories(spirulae_splat PRIVATE
    ${PROJECT_SOURCE_DIR}/glm
    ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
)
set_target_properties(spirulae_splat PROPERTIES CUDA_ARCHITECTURES "70;75;86")

# # To add an executable that uses the gsplat library,
# # follow example in the comments for a script `run_forward.cpp`
# # Add the executable
# add_executable(run_forward run_forward.cpp)

# # Link against CUDA runtime library
# target_link_libraries(run_forward PUBLIC cuda gsplat)

# # Include directories for the header-only library
# target_include_directories(run_forward PRIVATE
#     ${PROJECT_SOURCE_DIR}/glm
# )
