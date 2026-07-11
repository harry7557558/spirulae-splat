#pragma once

// Xml.h -- minimal dependency-free XML parser for the standalone CLI's
// Metashape dataset reader (camera-export .xml, .psx, doc.xml inside the
// .files zips). Header-only, recursive descent, no external deps; the
// ElementTree subset the Python parser relies on: elements, attributes,
// text, find (direct child) / findall / iter (recursive descendants).
//
// Not a general-purpose library: no namespaces, no DTD validation
// (DOCTYPE is skipped), input is assumed UTF-8, and parse errors throw
// std::runtime_error with a byte offset. Handles comments, processing
// instructions, CDATA, the five standard entities and numeric character
// references.

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct XmlNode {
    std::string tag;
    // Insertion-ordered; Metashape elements have few attributes.
    std::vector<std::pair<std::string, std::string>> attrs;
    std::string text;   // concatenated character data directly under this node
    std::vector<XmlNode> children;

    // Attribute value; nullptr when absent (ElementTree .get(name)).
    const std::string* attr(const std::string& name) const {
        for (const auto& [k, v] : attrs)
            if (k == name) return &v;
        return nullptr;
    }

    // First DIRECT child with the tag (ElementTree .find(tag)).
    const XmlNode* find(const std::string& t) const {
        for (const auto& c : children)
            if (c.tag == t) return &c;
        return nullptr;
    }

    // All DIRECT children with the tag (ElementTree .findall(tag)).
    std::vector<const XmlNode*> findall(const std::string& t) const {
        std::vector<const XmlNode*> out;
        for (const auto& c : children)
            if (c.tag == t) out.push_back(&c);
        return out;
    }

    // Self + all descendants with the tag, document order (ElementTree
    // .iter(tag)). Metashape nests <camera> inside <group> and <camera_ids>
    // inside <partition>, so the recursive walk is load-bearing.
    std::vector<const XmlNode*> iter(const std::string& t) const {
        std::vector<const XmlNode*> out;
        iter_into(t, out);
        return out;
    }

private:
    void iter_into(const std::string& t, std::vector<const XmlNode*>& out) const {
        if (tag == t) out.push_back(this);
        for (const auto& c : children) c.iter_into(t, out);
    }
};


namespace xml_detail {

struct Parser {
    const char* p;
    const char* end;
    const char* begin;

    [[noreturn]] void fail(const std::string& why) const {
        throw std::runtime_error("XML parse error at byte " +
                                 std::to_string(p - begin) + ": " + why);
    }
    void skip_ws() {
        while (p < end && (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r')) p++;
    }
    bool starts_with(const char* lit) const {
        size_t n = std::strlen(lit);
        return (size_t)(end - p) >= n && std::memcmp(p, lit, n) == 0;
    }
    void skip_until(const char* lit) {
        size_t n = std::strlen(lit);
        while ((size_t)(end - p) >= n) {
            if (std::memcmp(p, lit, n) == 0) { p += n; return; }
            p++;
        }
        fail(std::string("unterminated construct (expected '") + lit + "')");
    }

    static bool is_name_char(char c) {
        return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
               (c >= '0' && c <= '9') || c == '_' || c == '-' ||
               c == '.' || c == ':' || (unsigned char)c >= 0x80;
    }
    std::string parse_name() {
        const char* s = p;
        while (p < end && is_name_char(*p)) p++;
        if (p == s) fail("expected name");
        return std::string(s, p);
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

    // '&' already consumed.
    void parse_entity(std::string& out) {
        if (starts_with("amp;"))  { p += 4; out += '&'; return; }
        if (starts_with("lt;"))   { p += 3; out += '<'; return; }
        if (starts_with("gt;"))   { p += 3; out += '>'; return; }
        if (starts_with("quot;")) { p += 5; out += '"'; return; }
        if (starts_with("apos;")) { p += 5; out += '\''; return; }
        if (p < end && *p == '#') {
            p++;
            uint32_t cp = 0;
            if (p < end && (*p == 'x' || *p == 'X')) {
                p++;
                const char* s = p;
                while (p < end && *p != ';') {
                    char c = *p++;
                    cp <<= 4;
                    if (c >= '0' && c <= '9') cp |= (uint32_t)(c - '0');
                    else if (c >= 'a' && c <= 'f') cp |= (uint32_t)(c - 'a' + 10);
                    else if (c >= 'A' && c <= 'F') cp |= (uint32_t)(c - 'A' + 10);
                    else fail("bad hex character reference");
                }
                if (p == s) fail("empty character reference");
            } else {
                const char* s = p;
                while (p < end && *p != ';') {
                    char c = *p++;
                    if (c < '0' || c > '9') fail("bad character reference");
                    cp = cp * 10 + (uint32_t)(c - '0');
                }
                if (p == s) fail("empty character reference");
            }
            if (p >= end) fail("unterminated character reference");
            p++;  // ';'
            append_utf8(out, cp);
            return;
        }
        fail("unknown entity");
    }

    std::string parse_attr_value() {
        if (p >= end || (*p != '"' && *p != '\'')) fail("expected quoted attribute value");
        char q = *p++;
        std::string v;
        for (;;) {
            if (p >= end) fail("unterminated attribute value");
            char c = *p++;
            if (c == q) return v;
            if (c == '&') { parse_entity(v); continue; }
            v += c;
        }
    }

    // At '<' of an opening tag.
    XmlNode parse_element(int depth) {
        if (depth > 256) fail("nesting too deep");
        p++;  // '<'
        XmlNode node;
        node.tag = parse_name();
        for (;;) {
            skip_ws();
            if (p >= end) fail("unterminated tag");
            if (*p == '/') {
                p++;
                if (p >= end || *p != '>') fail("expected '>' after '/'");
                p++;
                return node;  // self-closing
            }
            if (*p == '>') { p++; break; }
            std::string key = parse_name();
            skip_ws();
            if (p >= end || *p != '=') fail("expected '=' after attribute name");
            p++;
            skip_ws();
            node.attrs.emplace_back(std::move(key), parse_attr_value());
        }
        // Content until the matching close tag.
        for (;;) {
            if (p >= end) fail("unterminated element <" + node.tag + ">");
            if (*p == '<') {
                if (starts_with("<!--")) { p += 4; skip_until("-->"); continue; }
                if (starts_with("<![CDATA[")) {
                    p += 9;
                    const char* s = p;
                    skip_until("]]>");
                    node.text.append(s, p - 3);
                    continue;
                }
                if (starts_with("<?")) { p += 2; skip_until("?>"); continue; }
                if (starts_with("</")) {
                    p += 2;
                    std::string close = parse_name();
                    if (close != node.tag)
                        fail("mismatched close tag </" + close + "> for <" + node.tag + ">");
                    skip_ws();
                    if (p >= end || *p != '>') fail("expected '>' in close tag");
                    p++;
                    return node;
                }
                node.children.push_back(parse_element(depth + 1));
                continue;
            }
            char c = *p++;
            if (c == '&') { parse_entity(node.text); continue; }
            node.text += c;
        }
    }

    // Prolog: BOM, XML declaration, comments, PIs, DOCTYPE.
    void skip_prolog() {
        if ((size_t)(end - p) >= 3 && (unsigned char)p[0] == 0xEF &&
            (unsigned char)p[1] == 0xBB && (unsigned char)p[2] == 0xBF) p += 3;
        for (;;) {
            skip_ws();
            if (starts_with("<?"))   { p += 2; skip_until("?>");  continue; }
            if (starts_with("<!--")) { p += 4; skip_until("-->"); continue; }
            if (starts_with("<!")) {   // DOCTYPE, possibly with [...] subset
                p += 2;
                int bracket = 0;
                while (p < end) {
                    char c = *p++;
                    if (c == '[') bracket++;
                    else if (c == ']') bracket--;
                    else if (c == '>' && bracket <= 0) break;
                }
                continue;
            }
            return;
        }
    }
};

}  // namespace xml_detail


inline XmlNode xml_parse(const std::string& text) {
    xml_detail::Parser ps{text.data(), text.data() + text.size(), text.data()};
    ps.skip_prolog();
    if (ps.p >= ps.end || *ps.p != '<') ps.fail("expected root element");
    XmlNode root = ps.parse_element(0);
    // Trailing comments/whitespace after the root are legal; ignore them.
    return root;
}

inline XmlNode xml_parse_file(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) throw std::runtime_error("cannot open " + path);
    std::fseek(f, 0, SEEK_END);
    long n = std::ftell(f);
    std::fseek(f, 0, SEEK_SET);
    std::string text(n > 0 ? (size_t)n : 0, '\0');
    size_t got = n > 0 ? std::fread(text.data(), 1, (size_t)n, f) : 0;
    std::fclose(f);
    text.resize(got);
    return xml_parse(text);
}
