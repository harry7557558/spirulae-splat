#pragma once

// HttpServer -- minimal dependency-free HTTP/1.0 server for the CLI viewer.
//
// Deliberately tiny: GET only, serial request handling on one background
// thread (the viewer protocol is polling + latest-wins, so serialization is
// part of the contract), Connection: close per request. Handlers run on the
// server thread; anything slow (rendering) hands off and waits.
//
// POSIX sockets on Linux/macOS; winsock on Windows (compiles, untested --
// cloud training boxes are Linux).

#include <atomic>
#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <thread>
#include <vector>

struct HttpRequest {
    std::string method;                          // "GET"
    std::string path;                            // "/render" (percent-decoded)
    std::map<std::string, std::string> query;    // decoded key -> value

    std::string get(const std::string& key, const std::string& def = "") const {
        auto it = query.find(key);
        return it == query.end() ? def : it->second;
    }
    double get_double(const std::string& key, double def) const;
    int    get_int(const std::string& key, int def) const;
    bool   get_bool(const std::string& key, bool def) const;
};

struct HttpResponse {
    int status = 200;
    std::string content_type = "text/plain";
    std::vector<uint8_t> body;

    static HttpResponse text(int status, const std::string& s,
                             const std::string& type = "text/plain") {
        HttpResponse r;
        r.status = status;
        r.content_type = type;
        r.body.assign(s.begin(), s.end());
        return r;
    }
    static HttpResponse json(const std::string& s) {
        return text(200, s, "application/json");
    }
};

class HttpServer {
public:
    using Handler = std::function<HttpResponse(const HttpRequest&)>;

    // Exact-path routes; unmatched paths get a 404.
    void route(const std::string& path, Handler h) { _routes[path] = std::move(h); }

    // Bind + start the accept loop on a background thread. Throws on bind
    // failure (port in use etc.).
    void start(const std::string& host, int port);
    void stop();

    ~HttpServer() { stop(); }

private:
    void serve_loop();
    std::map<std::string, Handler> _routes;
    std::atomic<bool> _running{false};
    std::thread _thread;
    intptr_t _listen_fd = -1;
};
