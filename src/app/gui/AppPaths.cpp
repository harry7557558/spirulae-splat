// AppPaths.cpp -- see AppPaths.h.

#include "app/gui/AppPaths.h"

#include <cstdlib>
#include <filesystem>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#else
#include <unistd.h>
#endif

namespace fs = std::filesystem;

namespace gui {

namespace {

fs::path env_path(const char* var) {
    const char* v = std::getenv(var);
    return v && *v ? fs::path(v) : fs::path();
}

fs::path ensure(fs::path dir) {
    std::error_code ec;
    fs::create_directories(dir, ec);
    return dir;
}

}  // namespace

std::string config_dir() {
#ifdef _WIN32
    fs::path dir = env_path("APPDATA");
#else
    fs::path dir = env_path("XDG_CONFIG_HOME");
    if (dir.empty()) {
        fs::path home = env_path("HOME");
        if (!home.empty()) dir = home / ".config";
    }
#endif
    if (dir.empty()) dir = ".";
    return ensure(dir / "spirulae-splat").string();
}

std::string cache_dir() {
#ifdef _WIN32
    fs::path dir = env_path("LOCALAPPDATA");
#else
    fs::path dir = env_path("XDG_CACHE_HOME");
    if (dir.empty()) {
        fs::path home = env_path("HOME");
        if (!home.empty()) dir = home / ".cache";
    }
#endif
    if (dir.empty()) dir = ".";
    return ensure(dir / "spirulae-splat").string();
}

std::string exe_dir() {
    static const std::string cached = [] {
        std::error_code ec;
#ifdef _WIN32
        wchar_t buf[32768];
        DWORD n = GetModuleFileNameW(nullptr, buf, (DWORD)std::size(buf));
        if (n == 0 || n >= std::size(buf)) return std::string();
        return fs::path(std::wstring(buf, n)).parent_path().string();
#elif defined(__APPLE__)
        uint32_t n = 0;
        _NSGetExecutablePath(nullptr, &n);
        std::string buf(n, '\0');
        if (_NSGetExecutablePath(buf.data(), &n) != 0) return std::string();
        fs::path p = fs::weakly_canonical(fs::path(buf.c_str()), ec);
        return ec ? std::string() : p.parent_path().string();
#else
        fs::path p = fs::read_symlink("/proc/self/exe", ec);
        return ec ? std::string() : p.parent_path().string();
#endif
    }();
    return cached;
}

std::string sibling_tool(const std::string& name) {
    const std::string dir = exe_dir();
    if (dir.empty()) return "";
#ifdef _WIN32
    fs::path p = fs::path(dir) / (name + ".exe");
#else
    fs::path p = fs::path(dir) / name;
#endif
    std::error_code ec;
    return fs::is_regular_file(p, ec) ? p.string() : std::string();
}

}  // namespace gui
