# run this from the `cuda` folder

compile_args="-ignore-capabilities -line-directive-mode none"

slangc slang/.slang -target cuda -o csrc/generated/slang.cuh $compile_args

slangc slang/all.slang -target cuda -o csrc/generated/slang_all.cu $compile_args

python3 slang/build_postprocess.py
