// Subprocess.cpp -- see Subprocess.h.

#include "app/gui/Subprocess.h"

#include <algorithm>
#include <cstring>
#include <filesystem>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#ifndef NOMINMAX
#define NOMINMAX   // keep windows.h from defining min()/max() macros
#endif
#include <windows.h>
#else
#include <fcntl.h>
#include <poll.h>
#include <signal.h>
#include <sys/wait.h>
#include <unistd.h>
#include <cstdlib>
#endif

namespace fs = std::filesystem;

namespace gui {

namespace {

// Feed a chunk into the line splitter. Both '\n' and a bare '\r' end a line
// -- ffmpeg/curl emit progress as '\r'-terminated updates that would
// otherwise never surface in the log (and look like a hang).
void emit_lines(std::string& acc, const char* buf, size_t n,
                const std::function<void(const std::string&)>& on_line) {
    acc.append(buf, n);
    size_t pos = 0, nl;
    while ((nl = acc.find_first_of("\r\n", pos)) != std::string::npos) {
        if (acc[nl] == '\r' && nl + 1 >= acc.size())
            break;   // might be the first half of a \r\n; wait for more
        size_t next = nl + 1;
        if (acc[nl] == '\r' && acc[next] == '\n') next++;
        if (nl > pos && on_line) on_line(acc.substr(pos, nl - pos));
        pos = next;
    }
    acc.erase(0, pos);
}

void emit_tail(std::string& acc, const std::function<void(const std::string&)>& on_line) {
    if (!acc.empty() && on_line) on_line(acc);
    acc.clear();
}

}  // namespace


#ifdef _WIN32

// ---------------------------------------------------------------------------
// Windows: CreateProcess + anonymous pipe, PeekNamedPipe polling.
// ---------------------------------------------------------------------------

namespace {

// Quote one argument per MSVCRT rules.
std::string quote_arg(const std::string& a) {
    if (!a.empty() && a.find_first_of(" \t\"") == std::string::npos) return a;
    std::string out = "\"";
    size_t bs = 0;
    for (char c : a) {
        if (c == '\\') { bs++; continue; }
        if (c == '"') { out.append(bs * 2 + 1, '\\'); out += '"'; bs = 0; continue; }
        out.append(bs, '\\'); bs = 0;
        out += c;
    }
    out.append(bs * 2, '\\');
    out += '"';
    return out;
}

}  // namespace

int run_process(const std::vector<std::string>& argv,
                const std::string& cwd,
                const std::function<void(const std::string&)>& on_line,
                const std::atomic<bool>& cancel) {
    if (argv.empty()) return kSpawnFailed;
    std::string cmdline;
    for (size_t i = 0; i < argv.size(); i++)
        cmdline += (i ? " " : "") + quote_arg(argv[i]);

    SECURITY_ATTRIBUTES sa{};
    sa.nLength = sizeof sa;
    sa.bInheritHandle = TRUE;
    HANDLE rd = nullptr, wr = nullptr;
    if (!CreatePipe(&rd, &wr, &sa, 0)) return kSpawnFailed;
    SetHandleInformation(rd, HANDLE_FLAG_INHERIT, 0);

    STARTUPINFOA si{};
    si.cb = sizeof si;
    si.dwFlags = STARTF_USESTDHANDLES;
    si.hStdOutput = wr;
    si.hStdError  = wr;
    si.hStdInput  = INVALID_HANDLE_VALUE;
    PROCESS_INFORMATION pi{};
    std::vector<char> cmd(cmdline.begin(), cmdline.end());
    cmd.push_back('\0');
    BOOL ok = CreateProcessA(nullptr, cmd.data(), nullptr, nullptr, TRUE,
                             CREATE_NO_WINDOW, nullptr,
                             cwd.empty() ? nullptr : cwd.c_str(), &si, &pi);
    CloseHandle(wr);
    if (!ok) { CloseHandle(rd); return kSpawnFailed; }
    CloseHandle(pi.hThread);

    std::string acc;
    char buf[4096];
    bool killed = false;
    for (;;) {
        if (cancel.load() && !killed) {
            TerminateProcess(pi.hProcess, 1);
            killed = true;
        }
        DWORD avail = 0;
        if (!PeekNamedPipe(rd, nullptr, 0, nullptr, &avail, nullptr)) break;
        if (avail > 0) {
            DWORD got = 0;
            if (!ReadFile(rd, buf, (DWORD)std::min<size_t>(sizeof buf, avail), &got, nullptr) || !got)
                break;
            emit_lines(acc, buf, got, on_line);
        } else {
            if (WaitForSingleObject(pi.hProcess, 100) == WAIT_OBJECT_0) {
                // Drain whatever is left.
                while (PeekNamedPipe(rd, nullptr, 0, nullptr, &avail, nullptr) && avail) {
                    DWORD got = 0;
                    if (!ReadFile(rd, buf, (DWORD)std::min<size_t>(sizeof buf, avail), &got, nullptr) || !got)
                        break;
                    emit_lines(acc, buf, got, on_line);
                }
                break;
            }
        }
    }
    emit_tail(acc, on_line);
    WaitForSingleObject(pi.hProcess, INFINITE);
    DWORD code = 1;
    GetExitCodeProcess(pi.hProcess, &code);
    CloseHandle(pi.hProcess);
    CloseHandle(rd);
    return killed ? kCancelled : (int)code;
}

bool command_exists(const std::string& exe) {
    if (exe.find('\\') != std::string::npos || exe.find('/') != std::string::npos) {
        std::error_code ec;
        return fs::exists(exe, ec);
    }
    char found[MAX_PATH];
    return SearchPathA(nullptr, exe.c_str(), ".exe", MAX_PATH, found, nullptr) > 0;
}

#else

// ---------------------------------------------------------------------------
// POSIX: fork/execvp + pipe, poll(2), process-group kill on cancel.
// ---------------------------------------------------------------------------

int run_process(const std::vector<std::string>& argv,
                const std::string& cwd,
                const std::function<void(const std::string&)>& on_line,
                const std::atomic<bool>& cancel) {
    if (argv.empty()) return kSpawnFailed;
    int fds[2];
    if (pipe(fds) != 0) return kSpawnFailed;

    pid_t pid = fork();
    if (pid < 0) { close(fds[0]); close(fds[1]); return kSpawnFailed; }
    if (pid == 0) {
        // Child: own process group so cancel kills helpers too.
        setpgid(0, 0);
        dup2(fds[1], STDOUT_FILENO);
        dup2(fds[1], STDERR_FILENO);
        close(fds[0]);
        close(fds[1]);
        if (!cwd.empty() && chdir(cwd.c_str()) != 0) _exit(127);
        std::vector<char*> args;
        for (const auto& a : argv) args.push_back(const_cast<char*>(a.c_str()));
        args.push_back(nullptr);
        execvp(args[0], args.data());
        _exit(127);
    }
    close(fds[1]);

    std::string acc;
    char buf[4096];
    bool killed = false;
    for (;;) {
        if (cancel.load() && !killed) {
            kill(-pid, SIGKILL);
            killed = true;
        }
        struct pollfd pfd{fds[0], POLLIN, 0};
        int pr = poll(&pfd, 1, 100);
        if (pr > 0) {
            ssize_t got = read(fds[0], buf, sizeof buf);
            if (got <= 0) break;   // EOF or error: child exited
            emit_lines(acc, buf, (size_t)got, on_line);
        } else if (pr < 0 && errno != EINTR) {
            break;
        }
    }
    emit_tail(acc, on_line);
    close(fds[0]);
    int status = 0;
    waitpid(pid, &status, 0);
    if (killed) return kCancelled;
    if (WIFEXITED(status)) {
        int code = WEXITSTATUS(status);
        return code == 127 ? kSpawnFailed : code;
    }
    return 1;
}

bool command_exists(const std::string& exe) {
    if (exe.find('/') != std::string::npos)
        return access(exe.c_str(), X_OK) == 0;
    const char* path = std::getenv("PATH");
    if (!path) return false;
    std::string p = path;
    size_t pos = 0;
    while (pos <= p.size()) {
        size_t colon = p.find(':', pos);
        std::string dir = p.substr(pos, colon == std::string::npos
                                            ? std::string::npos : colon - pos);
        if (!dir.empty() && access((dir + "/" + exe).c_str(), X_OK) == 0)
            return true;
        if (colon == std::string::npos) break;
        pos = colon + 1;
    }
    return false;
}

#endif

}  // namespace gui
