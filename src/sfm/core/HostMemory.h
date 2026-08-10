#pragma once
// How much physical memory the machine has, or 0 when the platform will not
// say. Callers that need a working number must supply their own fallback.

#include <cstddef>

#if defined(_WIN32)
#include <windows.h>
#else
#include <unistd.h>
#endif

namespace sfm {

inline size_t physicalRamBytes() {
#if defined(_WIN32)
    MEMORYSTATUSEX st{};
    st.dwLength = sizeof(st);
    if (GlobalMemoryStatusEx(&st)) return (size_t)st.ullTotalPhys;
#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)
    long pages = sysconf(_SC_PHYS_PAGES), page = sysconf(_SC_PAGESIZE);
    if (pages > 0 && page > 0) return (size_t)pages * (size_t)page;
#endif
    return 0;
}

}  // namespace sfm
