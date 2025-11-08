# run this from the `cuda` folder

# tested with slang-2025.19.1-linux-x86_64

compile_args="-ignore-capabilities -line-directive-mode none"
# compile_args="-ignore-capabilities"  # for debugging

slangc slang/.slang -target cuda -o csrc/generated/slang.cuh $compile_args

slangc slang/all.slang -target cuda -o csrc/generated/slang_all.cu $compile_args

slangc slang/primitive.slang -target cuda -o csrc/generated/primitive.cu $compile_args

python3 slang/build_postprocess.py
