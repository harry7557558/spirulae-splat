// The macOS half of NativeDialog.h: NSOpenPanel, which is main-thread only
// and runs modal, so the caller blocks until the panel is dismissed.

#import <AppKit/AppKit.h>

#include <string>
#include <vector>

namespace gui {

std::vector<std::string> mac_pick(const std::string& title, bool folder,
                                  const std::vector<std::string>& extensions,
                                  const std::string& start_dir, bool multi) {
    std::vector<std::string> out;
    @autoreleasepool {
        NSOpenPanel* panel = [NSOpenPanel openPanel];
        panel.canChooseFiles = !folder;
        panel.canChooseDirectories = folder;
        panel.canCreateDirectories = folder;
        panel.allowsMultipleSelection = multi && !folder;
        if (!title.empty())
            panel.message = [NSString stringWithUTF8String:title.c_str()];
        if (!start_dir.empty())
            panel.directoryURL = [NSURL fileURLWithPath:
                [NSString stringWithUTF8String:start_dir.c_str()]
                                           isDirectory:YES];

        if (!folder && !extensions.empty()) {
            // allowedFileTypes rather than allowedContentTypes: the newer one
            // needs the 12.0 SDK, and this deprecated spelling builds against
            // every SDK the rest of the program does.
            NSMutableArray* types = [NSMutableArray array];
            for (const std::string& e : extensions)
                [types addObject:[NSString stringWithUTF8String:
                    e.c_str() + (!e.empty() && e[0] == '.' ? 1 : 0)]];
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
            panel.allowedFileTypes = types;
#pragma clang diagnostic pop
        }

        if ([panel runModal] == NSModalResponseOK)
            for (NSURL* url in panel.URLs)
                if (const char* p = url.fileSystemRepresentation) out.push_back(p);
    }
    return out;
}

}  // namespace gui
