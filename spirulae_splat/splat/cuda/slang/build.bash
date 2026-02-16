# run this from the `cuda` folder

# tested with slang-2026.2.1-linux-x86_64

compile_args="-ignore-capabilities -line-directive-mode none"
# compile_args="-ignore-capabilities"  # for debugging

slangc slang/.slang -target cuda -o csrc/generated/slang.cuh $compile_args

slangc slang/all.slang -target cuda -o csrc/generated/slang_all.cu $compile_args

# slangc slang/primitive.slang -target cuda -o csrc/generated/primitive.cu $compile_args
slangc slang/projection_utils.slang -target cuda -o csrc/generated/projection_utils.cu $compile_args

slangc slang/primitive_3dgs.slang -target cuda -o csrc/generated/primitive_3dgs.cu $compile_args
slangc slang/primitive_3dgs_sv.slang -target cuda -o csrc/generated/primitive_3dgs_sv.cu $compile_args
# slangc slang/primitive_opaque_triangle.slang -target cuda -o csrc/generated/primitive_opaque_triangle.cu $compile_args
slangc slang/primitive_opaque_triangle_eval3d.slang -target cuda -o csrc/generated/primitive_opaque_triangle_eval3d.cu $compile_args
slangc slang/primitive_voxel.slang -target cuda -o csrc/generated/primitive_voxel.cu $compile_args

python3 slang/build_postprocess.py

mv csrc/generated/slang_all.cu csrc/generated/slang_all.cuh
mv csrc/generated/projection_utils.cu csrc/generated/projection_utils.cuh
mv csrc/generated/primitive_3dgs.cu csrc/generated/primitive_3dgs.cuh
mv csrc/generated/primitive_3dgs_sv.cu csrc/generated/primitive_3dgs_sv.cuh
mv csrc/generated/primitive_opaque_triangle_eval3d.cu csrc/generated/primitive_opaque_triangle_eval3d.cuh
mv csrc/generated/primitive_voxel.cu csrc/generated/primitive_voxel.cuh

content=$(cat <<'EOF'
using ::FixedArray;
using ::Matrix;

using ::__ldg;
using ::atomicAdd;
using ::atomicCAS;
using ::atomicExch;

#define _IMPORT_GLOBAL_FUNC(dtype) \
    using ::make_##dtype##2; \
    using ::make_##dtype##3; \
    using ::make_##dtype##4;

_IMPORT_GLOBAL_FUNC(int)
//_IMPORT_GLOBAL_FUNC(bool)
_IMPORT_GLOBAL_FUNC(uint)
_IMPORT_GLOBAL_FUNC(short)
_IMPORT_GLOBAL_FUNC(ushort)
_IMPORT_GLOBAL_FUNC(char)
_IMPORT_GLOBAL_FUNC(uchar)
_IMPORT_GLOBAL_FUNC(longlong)
_IMPORT_GLOBAL_FUNC(ulonglong)
_IMPORT_GLOBAL_FUNC(float)
_IMPORT_GLOBAL_FUNC(double)

#undef _IMPORT_GLOBAL_FUNC
EOF
)
echo "$content" > csrc/generated/set_namespace.cuh
