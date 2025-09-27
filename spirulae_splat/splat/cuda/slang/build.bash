# run this from the `cuda` folder

compile_args="-ignore-capabilities -line-directive-mode none"

slangc slang/all.slang -target cuda -o csrc/generated/slang_all.cu $compile_args
