#pragma once

// Json.h -- minimal dependency-free JSON parser for the standalone CLI's
// dataset readers (transforms.json etc.). Header-only, recursive descent,
// no external deps. Supports the full JSON grammar plus two Python
// json.dump extensions that appear in this project's files: Infinity /
// -Infinity and NaN literals.
//
// Not a general-purpose library: values are held by copy (fine for
// transforms.json-sized inputs), numbers are always double, and parse
// errors throw std::runtime_error with a byte offset.

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct JsonValue {
    enum class Type { Null, Bool, Number, String, Array, Object };
    Type type = Type::Null;

    bool        b   = false;
    double      num = 0.0;
    std::string str;
    std::vector<JsonValue> arr;
    // Insertion-ordered; transforms.json keys are few, linear find is fine.
    std::vector<std::pair<std::string, JsonValue>> obj;

    bool is_null()   const { return type == Type::Null; }
    bool is_object() const { return type == Type::Object; }
    bool is_array()  const { return type == Type::Array; }

    // Object lookup; nullptr when absent (or not an object).
    const JsonValue* find(const std::string& key) const {
        if (type != Type::Object) return nullptr;
        for (const auto& [k, v] : obj)
            if (k == key) return &v;
        return nullptr;
    }
    bool has(const std::string& key) const { return find(key) != nullptr; }

    // Typed accessors with defaults. Numbers stored in JSON as ints are
    // doubles here; as_int truncates.
    double as_double(double def = 0.0) const {
        return type == Type::Number ? num : def;
    }
    int64_t as_int(int64_t def = 0) const {
        return type == Type::Number ? (int64_t)num : def;
    }
    bool as_bool(bool def = false) const {
        return type == Type::Bool ? b : def;
    }
    const std::string& as_string() const {
        static const std::string empty;
        return type == Type::String ? str : empty;
    }
    double get_double(const std::string& key, double def) const {
        const JsonValue* v = find(key);
        return v ? v->as_double(def) : def;
    }
};


namespace json_detail {

struct Parser {
    const char* p;
    const char* end;
    const char* begin;

    [[noreturn]] void fail(const std::string& why) const {
        throw std::runtime_error("JSON parse error at byte " +
                                 std::to_string(p - begin) + ": " + why);
    }
    void skip_ws() {
        while (p < end && (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r')) p++;
    }
    char peek() {
        if (p >= end) fail("unexpected end of input");
        return *p;
    }
    bool consume(const char* lit) {
        size_t n = std::strlen(lit);
        if ((size_t)(end - p) >= n && std::memcmp(p, lit, n) == 0) { p += n; return true; }
        return false;
    }

    void append_utf8(std::string& s, uint32_t cp) {
        if (cp < 0x80) s += (char)cp;
        else if (cp < 0x800) {
            s += (char)(0xC0 | (cp >> 6));
            s += (char)(0x80 | (cp & 0x3F));
        } else if (cp < 0x10000) {
            s += (char)(0xE0 | (cp >> 12));
            s += (char)(0x80 | ((cp >> 6) & 0x3F));
            s += (char)(0x80 | (cp & 0x3F));
        } else {
            s += (char)(0xF0 | (cp >> 18));
            s += (char)(0x80 | ((cp >> 12) & 0x3F));
            s += (char)(0x80 | ((cp >> 6) & 0x3F));
            s += (char)(0x80 | (cp & 0x3F));
        }
    }
    uint32_t parse_hex4() {
        if (end - p < 4) fail("truncated \\u escape");
        uint32_t v = 0;
        for (int i = 0; i < 4; i++) {
            char c = *p++;
            v <<= 4;
            if (c >= '0' && c <= '9') v |= (uint32_t)(c - '0');
            else if (c >= 'a' && c <= 'f') v |= (uint32_t)(c - 'a' + 10);
            else if (c >= 'A' && c <= 'F') v |= (uint32_t)(c - 'A' + 10);
            else fail("bad \\u escape");
        }
        return v;
    }

    std::string parse_string() {
        if (*p != '"') fail("expected string");
        p++;
        std::string s;
        for (;;) {
            if (p >= end) fail("unterminated string");
            char c = *p++;
            if (c == '"') return s;
            if (c != '\\') { s += c; continue; }
            if (p >= end) fail("truncated escape");
            char e = *p++;
            switch (e) {
                case '"': s += '"'; break;
                case '\\': s += '\\'; break;
                case '/': s += '/'; break;
                case 'b': s += '\b'; break;
                case 'f': s += '\f'; break;
                case 'n': s += '\n'; break;
                case 'r': s += '\r'; break;
                case 't': s += '\t'; break;
                case 'u': {
                    uint32_t cp = parse_hex4();
                    if (cp >= 0xD800 && cp <= 0xDBFF) {   // surrogate pair
                        if (end - p < 2 || p[0] != '\\' || p[1] != 'u')
                            fail("unpaired surrogate");
                        p += 2;
                        uint32_t lo = parse_hex4();
                        if (lo < 0xDC00 || lo > 0xDFFF) fail("bad low surrogate");
                        cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00);
                    }
                    append_utf8(s, cp);
                    break;
                }
                default: fail("unknown escape");
            }
        }
    }

    JsonValue parse_number() {
        char* num_end = nullptr;
        double v = std::strtod(p, &num_end);
        if (num_end == p) fail("bad number");
        p = num_end;
        JsonValue out;
        out.type = JsonValue::Type::Number;
        out.num = v;
        return out;
    }

    JsonValue parse_value(int depth) {
        if (depth > 256) fail("nesting too deep");
        skip_ws();
        char c = peek();
        JsonValue out;
        if (c == '{') {
            p++;
            out.type = JsonValue::Type::Object;
            skip_ws();
            if (peek() == '}') { p++; return out; }
            for (;;) {
                skip_ws();
                std::string key = parse_string();
                skip_ws();
                if (peek() != ':') fail("expected ':'");
                p++;
                out.obj.emplace_back(std::move(key), parse_value(depth + 1));
                skip_ws();
                char d = peek();
                if (d == ',') { p++; continue; }
                if (d == '}') { p++; return out; }
                fail("expected ',' or '}'");
            }
        }
        if (c == '[') {
            p++;
            out.type = JsonValue::Type::Array;
            skip_ws();
            if (peek() == ']') { p++; return out; }
            for (;;) {
                out.arr.push_back(parse_value(depth + 1));
                skip_ws();
                char d = peek();
                if (d == ',') { p++; continue; }
                if (d == ']') { p++; return out; }
                fail("expected ',' or ']'");
            }
        }
        if (c == '"') {
            out.type = JsonValue::Type::String;
            out.str = parse_string();
            return out;
        }
        if (consume("true"))  { out.type = JsonValue::Type::Bool; out.b = true;  return out; }
        if (consume("false")) { out.type = JsonValue::Type::Bool; out.b = false; return out; }
        if (consume("null"))  { return out; }
        // Python json.dump extensions.
        if (consume("Infinity"))  { out.type = JsonValue::Type::Number; out.num = std::numeric_limits<double>::infinity();  return out; }
        if (consume("-Infinity")) { out.type = JsonValue::Type::Number; out.num = -std::numeric_limits<double>::infinity(); return out; }
        if (consume("NaN"))       { out.type = JsonValue::Type::Number; out.num = std::nan("");  return out; }
        if (c == '-' || (c >= '0' && c <= '9')) return parse_number();
        fail("unexpected character");
    }
};

}  // namespace json_detail


inline JsonValue json_parse(const std::string& text) {
    json_detail::Parser ps{text.data(), text.data() + text.size(), text.data()};
    JsonValue v = ps.parse_value(0);
    ps.skip_ws();
    if (ps.p != ps.end) ps.fail("trailing content");
    return v;
}

inline JsonValue json_parse_file(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) throw std::runtime_error("cannot open " + path);
    std::fseek(f, 0, SEEK_END);
    long n = std::ftell(f);
    std::fseek(f, 0, SEEK_SET);
    std::string text(n > 0 ? (size_t)n : 0, '\0');
    size_t got = n > 0 ? std::fread(text.data(), 1, (size_t)n, f) : 0;
    std::fclose(f);
    text.resize(got);
    return json_parse(text);
}
