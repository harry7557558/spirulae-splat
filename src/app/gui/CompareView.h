#pragma once

// CompareView -- up to four finished models under one camera, laid out in one
// pane, two, three or a 2x2 grid. The viewer screen and the meshing preview
// are both this. Design: docs/notes/compare-view.md.
//
// It owns the engine while it is open (engine_reset at both ends), so opening
// it is a session-destroying action for GuiApp and nothing else may render
// from the engine meanwhile. GUI thread only.

#include "app/gui/SplatViewer.h"
#include "app/gui/ViewportPanel.h"
#include "i18n/Message.h"

#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace gui {

class CompareView {
public:
    // Four is where a pane stops being big enough to judge a render in.
    static constexpr int kMaxModels = 4;

    ~CompareView();

    // Take the engine over and show `path` alone.
    void open(const std::string& path);
    // Beside what is already open; ignored when full. `title` names the pane
    // when the model has a role rather than just a filename (the meshing
    // screen's two sides); otherwise the file's own name is used.
    void add(const std::string& path, const spirula::i18n::Msg* title = nullptr);
    // Hand the engine back. Safe when nothing is open.
    void close();

    // Whether the engine is this object's: true from the first add() until
    // close(), which is what the owner keys the viewport handover on -- the
    // last model being closed does not hand the engine back.
    bool holds_engine() const { return _engine_taken; }
    int count() const { return (int)_models.size(); }
    bool full() const { return count() >= kMaxModels; }
    // Recent paths offered by the "add" menu, refreshed by the owner each
    // frame -- this object has no settings of its own.
    void set_recents(const std::vector<std::string>& r) { _recents = r; }
    // Called when the user asks for a file picker (the owner owns the dialog).
    void set_pick_file(std::function<void()> f) { _pick_file = std::move(f); }

    // Attach panels whose loaders have finished and refresh the placements.
    // Once a frame, before draw().
    void poll();
    // The row above the panes: what is open, the add button, the link.
    void draw_toolbar();
    // The panes themselves, `height` tall (0 = all that is left).
    void draw(float height);
    void destroy_gl();
    std::vector<std::string> drain_log();

private:
    struct Model {
        SplatViewer src;
        ViewportPanel panel;
        const spirula::i18n::Msg* title = nullptr;
        int slot = -1;              // engine scene slot
        bool attached = false;
        // Placement in the shared frame. `align` puts the model in the FIRST
        // model's frame rather than in its own; the rest is the hand
        // adjustment on top of that.
        bool align = true;
        float offset[3] = {0, 0, 0};
        float euler[3] = {0, 0, 0};   // degrees, applied X then Y then Z
        float scale = 1.0f;
    };

    void take_engine();
    void remove(int index);
    void move(int index, int dir);
    int  claim_slot();
    void attach(Model& m);
    // Recompute every pane's model->shared similarity from the placements.
    void update_placements();
    // engine_blit_view refuses to run before engine_viewer_init, and the grid
    // it draws belongs to whichever model defines the shared frame.
    void ensure_viewer_overlay();
    void draw_pane(int index, const ImVec2& size);
    void draw_placement_popup(int index);
    // The panel driving the link this frame: the one being dragged (sticky
    // for the length of a drag), else whichever moved.
    ViewportPanel* link_master();

    // One mutex for every panel's render worker: the engine binds one scene
    // at a time, so two panels rendering at once would race over which.
    // Declared before the models, which hold a pointer to it.
    std::mutex _engine_mutex;
    std::vector<std::unique_ptr<Model>> _models;
    bool _engine_taken = false;
    bool _link = true;
    // Whose frame the axes/grid overlay was built for; "" when there is none.
    std::string _overlay_key;
    uint32_t _slots_used = 0;
    // Pane-menu actions, applied by the next poll(): a pane cannot remove or
    // reorder itself while its own popup is being drawn inside it.
    int _pending_remove = -1;
    int _pending_move = 0;
    int _pending_move_index = -1;
    // Tallest control block across the panes on the last draw; the shorter
    // ones are padded to it so every image comes out the same size.
    float _controls_h = 0.0f;

    std::vector<std::string> _recents;
    std::function<void()> _pick_file;
};

}  // namespace gui
