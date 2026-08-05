#pragma once
// main() for the standalone SfM test binaries.
//
// Everything these tests exercise can throw: VkContext::init reports a device
// that lacks fp64, an atomic add or the integer dot product by throwing a
// std::runtime_error naming the missing feature, which is exactly the message a
// cross-vendor run exists to read. glibc's default terminate handler prints
// what() before aborting, so on Linux that message survives an uncaught throw.
// MSVC's does not: it hands the process to Windows Error Reporting, which on a
// headless box leaves a WerFault.exe sitting on the test with nothing printed
// and the batch that launched it wedged behind it. Catch here so the failure
// reads the same on every platform.
//
//   int main(int argc, char** argv) { return sfmTestMain(argc, argv, cmdFoo); }

#include <cstdio>
#include <exception>

inline int sfmTestMain(int argc, char** argv, int (*body)(int, char**)) {
    try {
        return body(argc, argv);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "FAIL: %s\n", e.what());
        return 2;
    } catch (...) {
        std::fprintf(stderr, "FAIL: unknown exception\n");
        return 2;
    }
}
