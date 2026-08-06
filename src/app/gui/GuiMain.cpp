// GuiMain.cpp -- the GUI entry point: GLFW window + OpenGL 3.2
// core context (lowest version that also satisfies macOS) + Dear ImGui
// bootstrap, then hands every frame to gui::GuiApp.

#include "app/Tools.h"
#include "app/gui/Fonts.h"
#include "app/gui/GuiApp.h"
#include "i18n/catalog/Brand.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include <GLFW/glfw3.h>

#include <cstdio>
#include <string>
#include <vector>

namespace {

void glfw_error_callback(int error, const char* description) {
    std::fprintf(stderr, "[glfw] error %d: %s\n", error, description);
}

void apply_style(float scale) {
    ImGuiStyle& style = ImGui::GetStyle();
    ImGui::StyleColorsDark();
    style.WindowRounding = 4.0f;
    style.FrameRounding = 3.0f;
    style.GrabRounding = 3.0f;
    style.TabRounding = 3.0f;
    style.ScrollbarRounding = 3.0f;
    style.FramePadding = ImVec2(6, 4);
    style.ItemSpacing = ImVec2(8, 6);
    ImVec4* c = style.Colors;
    c[ImGuiCol_WindowBg]         = ImVec4(0.10f, 0.105f, 0.12f, 1.0f);
    c[ImGuiCol_ChildBg]          = ImVec4(0.115f, 0.12f, 0.135f, 1.0f);
    c[ImGuiCol_FrameBg]          = ImVec4(0.19f, 0.20f, 0.23f, 1.0f);
    c[ImGuiCol_FrameBgHovered]   = ImVec4(0.24f, 0.25f, 0.29f, 1.0f);
    c[ImGuiCol_FrameBgActive]    = ImVec4(0.28f, 0.30f, 0.34f, 1.0f);
    c[ImGuiCol_Button]           = ImVec4(0.22f, 0.32f, 0.44f, 1.0f);
    c[ImGuiCol_ButtonHovered]    = ImVec4(0.28f, 0.40f, 0.55f, 1.0f);
    c[ImGuiCol_ButtonActive]     = ImVec4(0.33f, 0.47f, 0.64f, 1.0f);
    c[ImGuiCol_Header]           = ImVec4(0.20f, 0.25f, 0.32f, 1.0f);
    c[ImGuiCol_HeaderHovered]    = ImVec4(0.26f, 0.32f, 0.41f, 1.0f);
    c[ImGuiCol_HeaderActive]     = ImVec4(0.30f, 0.37f, 0.48f, 1.0f);
    c[ImGuiCol_PlotLines]        = ImVec4(0.45f, 0.75f, 1.0f, 1.0f);
    c[ImGuiCol_PlotHistogram]    = ImVec4(0.30f, 0.55f, 0.85f, 1.0f);  // also the ProgressBar fill
    style.ScaleAllSizes(scale);
    style.FontScaleMain = scale;
}

}  // namespace

int spirula_gui_main(int argc, char** argv) {
    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit()) {
        std::fprintf(stderr, "error: failed to initialize GLFW (is a display "
                     "available?)\n");
        return 1;
    }

    // GL 3.2 core: works on Linux/Windows and is the macOS minimum.
    const char* glsl_version = "#version 150";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);
#endif

    // The title is set from the catalog below, once GuiApp has settled the
    // language; this is only what the window is born with.
    GLFWwindow* window = glfwCreateWindow(1600, 950, "Spirula Studio",
                                          nullptr, nullptr);
    if (!window) {
        std::fprintf(stderr, "error: failed to create window / GL context\n");
        glfwTerminate();
        return 1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);   // vsync

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.IniFilename = nullptr;   // fixed layout; nothing worth persisting

    float xscale = 1.0f, yscale = 1.0f;
    glfwGetWindowContentScale(window, &xscale, &yscale);
    apply_style(xscale > 0.0f ? xscale : 1.0f);

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

    {
        gui::GuiApp app;
        // The atlas has to exist before the first frame, and the language it
        // depends on is only settled once GuiApp has read its settings file.
        app.fonts().ensure();
        // Drag-and-drop onto the window: auto-detect dataset folder / photo
        // folder / video file (fires on the main thread via glfwPollEvents).
        // The whole drop is handed over at once -- several videos dropped
        // together are the inputs of one dataset, not four separate answers.
        glfwSetWindowUserPointer(window, &app);
        glfwSetDropCallback(window, [](GLFWwindow* w, int count,
                                       const char** paths) {
            auto* a = (gui::GuiApp*)glfwGetWindowUserPointer(w);
            if (!a || !paths) return;
            std::vector<std::string> dropped;
            for (int i = 0; i < count; i++)
                if (paths[i] && paths[i][0]) dropped.push_back(paths[i]);
            if (!dropped.empty()) a->handle_drop(dropped);
        });
        // `spirula-gui <path>...` -- the same auto-detection as a drop, so a
        // desktop "Open with" or a shell alias lands on the right screen.
        {
            std::vector<std::string> args;
            for (int i = 1; i < argc; i++)
                if (argv[i] && argv[i][0]) args.push_back(argv[i]);
            if (!args.empty()) app.handle_drop(args);
        }
        // The title bar is interface copy too, and the language it is in is
        // only settled after GuiApp has read its settings -- and changes again
        // whenever the picker does. Comparing the string is enough of a
        // change test: it is immortal .rodata that the current language
        // selects, so this is one pointer-deep compare per frame.
        std::string title;
        auto sync_title = [&] {
            const char* want = spirula::i18n::msg::brand::window_title.get();
            if (title == want) return;
            title = want;
            glfwSetWindowTitle(window, title.c_str());
        };
        sync_title();

        while (!app.wants_exit()) {
            glfwPollEvents();
            if (glfwWindowShouldClose(window)) {
                glfwSetWindowShouldClose(window, GLFW_FALSE);
                app.request_close();
            }

            // Between frames, and before NewFrame: swapping a font face
            // invalidates every ImFont pointer, and this is the only point in
            // the loop where none is live. Returns immediately unless the
            // language changed or a font download finished.
            app.fonts().ensure();
            sync_title();

            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();
            app.frame();
            ImGui::Render();

            int w = 0, h = 0;
            glfwGetFramebufferSize(window, &w, &h);
            glViewport(0, 0, w, h);
            glClearColor(0.07f, 0.07f, 0.08f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
            glfwSwapBuffers(window);
        }
        // Joins worker threads (finishing a final checkpoint save if a stop
        // was just requested) and frees GL resources while the context is
        // still current.
        app.shutdown();
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
