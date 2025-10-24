# run this from the `cuda` folder

# tested with slang-2025.19.1-linux-x86_64

compile_args="-ignore-capabilities -line-directive-mode none"

slangc slang/.slang -target cuda -o csrc/generated/slang.cuh $compile_args

slangc slang/all.slang -target cuda -o csrc/generated/slang_all.cu $compile_args

# slangc slang/projection_3dgs.slang -target cuda -o csrc/generated/projection_3dgs.cu $compile_args
# slangc slang/projection_opaque_triangle.slang -target cuda -o csrc/generated/projection_opaque_triangle.cu $compile_args
slangc slang/projection.slang -target cuda -o csrc/generated/projection.cu $compile_args

python3 slang/build_postprocess.py
