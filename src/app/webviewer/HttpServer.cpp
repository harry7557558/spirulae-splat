// HttpServer.cpp -- see HttpServer.h.

#include "app/webviewer/HttpServer.h"

#include <cstdio>
#include <cstring>
#include <stdexcept>

#ifdef _WIN32
  #include <winsock2.h>
  #include <ws2tcpip.h>
  #pragma comment(lib, "ws2_32.lib")
  typedef SOCKET socket_t;
  static const socket_t kInvalidSocket = INVALID_SOCKET;
  static void close_socket(socket_t s) { closesocket(s); }
  static bool socket_startup() {
      static bool done = false;
      if (!done) {
          WSADATA wsa;
          if (WSAStartup(MAKEWORD(2, 2), &wsa) != 0) return false;
          done = true;
      }
      return true;
  }
#else
  #include <arpa/inet.h>
  #include <netinet/in.h>
  #include <netinet/tcp.h>
  #include <sys/select.h>
  #include <sys/socket.h>
  #include <unistd.h>
  typedef int socket_t;
  static const socket_t kInvalidSocket = -1;
  static void close_socket(socket_t s) { ::close(s); }
  static bool socket_startup() { return true; }
#endif


// ===========================================================================
// Request helpers
// ===========================================================================

double HttpRequest::get_double(const std::string& key, double def) const {
    auto it = query.find(key);
    if (it == query.end()) return def;
    try { return std::stod(it->second); } catch (...) { return def; }
}
int HttpRequest::get_int(const std::string& key, int def) const {
    auto it = query.find(key);
    if (it == query.end()) return def;
    try { return std::stoi(it->second); } catch (...) { return def; }
}
bool HttpRequest::get_bool(const std::string& key, bool def) const {
    auto it = query.find(key);
    if (it == query.end()) return def;
    const std::string& v = it->second;
    return v == "1" || v == "true" || v == "True" || v == "yes";
}

namespace {

std::string url_decode(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (size_t i = 0; i < s.size(); i++) {
        if (s[i] == '+') { out += ' '; continue; }
        if (s[i] == '%' && i + 2 < s.size()) {
            auto hex = [](char c) -> int {
                if (c >= '0' && c <= '9') return c - '0';
                if (c >= 'a' && c <= 'f') return c - 'a' + 10;
                if (c >= 'A' && c <= 'F') return c - 'A' + 10;
                return -1;
            };
            int hi = hex(s[i+1]), lo = hex(s[i+2]);
            if (hi >= 0 && lo >= 0) {
                out += (char)((hi << 4) | lo);
                i += 2;
                continue;
            }
        }
        out += s[i];
    }
    return out;
}

// Read until the end of the request headers (\r\n\r\n). We only need the
// request line; the rest is drained and discarded (GET-only protocol).
bool read_request_head(socket_t fd, std::string& head) {
    char buf[2048];
    head.clear();
    while (head.find("\r\n\r\n") == std::string::npos) {
        if (head.size() > 65536) return false;
#ifdef _WIN32
        int n = recv(fd, buf, sizeof buf, 0);
#else
        ssize_t n = recv(fd, buf, sizeof buf, 0);
#endif
        if (n <= 0) return false;
        head.append(buf, (size_t)n);
    }
    return true;
}

bool send_all(socket_t fd, const void* data, size_t len) {
    const char* p = (const char*)data;
    while (len > 0) {
#ifdef _WIN32
        int n = send(fd, p, (int)len, 0);
#else
        ssize_t n = send(fd, p, len, MSG_NOSIGNAL);
#endif
        if (n <= 0) return false;
        p += n;
        len -= (size_t)n;
    }
    return true;
}

const char* status_text(int code) {
    switch (code) {
        case 200: return "OK";
        case 400: return "Bad Request";
        case 404: return "Not Found";
        case 500: return "Internal Server Error";
        default:  return "";
    }
}

}  // namespace


// ===========================================================================
// HttpServer
// ===========================================================================

void HttpServer::start(const std::string& host, int port) {
    if (!socket_startup())
        throw std::runtime_error("HttpServer: socket subsystem init failed");

    socket_t fd = ::socket(AF_INET, SOCK_STREAM, 0);
    if (fd == kInvalidSocket)
        throw std::runtime_error("HttpServer: socket() failed");

    int one = 1;
    setsockopt(fd, SOL_SOCKET, SO_REUSEADDR, (const char*)&one, sizeof one);

    sockaddr_in addr{};
    addr.sin_family = AF_INET;
    addr.sin_port = htons((uint16_t)port);
    if (host.empty() || host == "0.0.0.0") addr.sin_addr.s_addr = INADDR_ANY;
    else if (inet_pton(AF_INET, host.c_str(), &addr.sin_addr) != 1) {
        close_socket(fd);
        throw std::runtime_error("HttpServer: bad host " + host);
    }
    if (bind(fd, (sockaddr*)&addr, sizeof addr) != 0) {
        close_socket(fd);
        throw std::runtime_error("HttpServer: cannot bind port " + std::to_string(port));
    }
    if (listen(fd, 8) != 0) {
        close_socket(fd);
        throw std::runtime_error("HttpServer: listen() failed");
    }

    _listen_fd = (intptr_t)fd;
    _running = true;
    _thread = std::thread(&HttpServer::serve_loop, this);
}

void HttpServer::stop() {
    if (!_running.exchange(false)) return;
    // Do NOT close the listening socket to break accept(): closing an fd that
    // another thread is blocked in accept() on is undefined by POSIX, and on
    // Linux the blocked accept() simply never returns -- which hung every
    // ssplat-train run that had the viewer enabled and keep_viewer_alive off.
    // serve_loop() select()s with a timeout instead, so clearing _running is
    // enough; it exits within one poll interval and we close after the join.
    if (_thread.joinable()) _thread.join();
    if (_listen_fd != -1) {
        close_socket((socket_t)_listen_fd);
        _listen_fd = -1;
    }
}

void HttpServer::serve_loop() {
    // Poll interval: the worst-case delay between stop() and this thread
    // noticing. Short enough to feel instant, long enough to be free.
    static constexpr long kPollUs = 100000;   // 100 ms

    while (_running) {
        fd_set readfds;
        FD_ZERO(&readfds);
        FD_SET((socket_t)_listen_fd, &readfds);
        timeval tv{};
        tv.tv_sec = 0;
        tv.tv_usec = kPollUs;

        int ready = select((int)_listen_fd + 1, &readfds, nullptr, nullptr, &tv);
        if (ready <= 0) continue;   // timeout, or interrupted -- re-check _running

        socket_t client = accept((socket_t)_listen_fd, nullptr, nullptr);
        if (client == kInvalidSocket) {
            if (_running)
                continue;    // transient accept error
            break;           // shutting down
        }

        std::string head;
        HttpResponse resp;
        if (!read_request_head(client, head)) {
            close_socket(client);
            continue;
        }

        // Request line: METHOD SP target SP version.
        size_t eol = head.find("\r\n");
        std::string line = head.substr(0, eol);
        size_t sp1 = line.find(' '), sp2 = line.rfind(' ');
        if (sp1 == std::string::npos || sp2 == sp1) {
            resp = HttpResponse::text(400, "bad request");
        } else {
            HttpRequest req;
            req.method = line.substr(0, sp1);
            std::string target = line.substr(sp1 + 1, sp2 - sp1 - 1);
            size_t qpos = target.find('?');
            req.path = url_decode(target.substr(0, qpos));
            if (qpos != std::string::npos) {
                std::string qs = target.substr(qpos + 1);
                size_t pos = 0;
                while (pos <= qs.size()) {
                    size_t amp = qs.find('&', pos);
                    std::string kv = qs.substr(
                        pos, amp == std::string::npos ? std::string::npos : amp - pos);
                    size_t eq = kv.find('=');
                    if (eq != std::string::npos)
                        req.query[url_decode(kv.substr(0, eq))] = url_decode(kv.substr(eq + 1));
                    else if (!kv.empty())
                        req.query[url_decode(kv)] = "";
                    if (amp == std::string::npos) break;
                    pos = amp + 1;
                }
            }

            auto it = _routes.find(req.path);
            if (req.method != "GET") {
                resp = HttpResponse::text(400, "GET only");
            } else if (it == _routes.end()) {
                resp = HttpResponse::text(404, "not found");
            } else {
                try {
                    resp = it->second(req);
                } catch (const std::exception& e) {
                    std::fprintf(stderr, "[viewer] handler error: %s\n", e.what());
                    resp = HttpResponse::text(500, e.what());
                }
            }
        }

        char hdr[256];
        int n = std::snprintf(hdr, sizeof hdr,
            "HTTP/1.0 %d %s\r\nContent-Type: %s\r\nContent-Length: %zu\r\n"
            "Access-Control-Allow-Origin: *\r\nConnection: close\r\n\r\n",
            resp.status, status_text(resp.status),
            resp.content_type.c_str(), resp.body.size());
        // Client disconnects mid-response (tab refresh, latest-wins drop) are
        // routine; send_all just fails and we move on -- parity with the
        // Python handler's BrokenPipeError swallowing.
        if (send_all(client, hdr, (size_t)n) && !resp.body.empty())
            send_all(client, resp.body.data(), resp.body.size());
        close_socket(client);
    }
}
