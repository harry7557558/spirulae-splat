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
