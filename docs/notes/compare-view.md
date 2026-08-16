# Looking at several models at once

`src/app/gui/CompareView.{h,cpp}` is the GUI's model viewer: up to four
finished models — splats, meshes, point clouds, in any mix — laid out in one
pane, two, three, or a 2×2 grid, on one camera. It is what both the viewer
screen and the meshing preview draw.

Looking at one model answers "is this good?". Looking at two answers "is this
better?", which is the question anyone tuning a reconstruction is actually
asking, and it is only answerable if the two pictures differ in nothing but the
model.

Three pieces make that work.

## Several splat models resident at once

The engine is a process-global singleton with one set of world splats
(`AGENTS.md`, "the engine is a process-global singleton"). Re-uploading a set
per pane is not an option: a 274k-splat model with SH degree 3 is 65 MB on
device, and a camera drag would re-upload every model every frame.

So `Engine.h` grows **scene slots**: `engine_scene_set_data_3dgs(slot, …)`
uploads a complete splat set into its own device buffers, and
`engine_scene_activate(slot)` points `engine().world` at one of them.
Activation is a pointer swap plus four scalars, so `RenderWorker` does it per
render (`ViewerRenderConfig::scene_slot`) rather than leaving it to whoever
touched the engine last.

The buffers go through `DevicePool::acquire_dynamic` with a `scene<i>.` name
prefix rather than through a `PoolSlot`, because the number of live sets is a
runtime property — the compile-time slot table has one `world.means`, and that
one is still the training run's. `engine_reset()` frees them and forgets the
table together, so a slot cannot outlive the memory under it.

A slot carries no optimizer, densify or gradient state. It is for looking at,
not for training on.

Everything else the panes share stays shared: one `ForwardCache`, one set of
screen buffers, one BVH. They are scratch, resized per render, so the pool
settles at the high-water mark of the largest model rather than the sum.

## One camera, four frames

Each pane carries a similarity mapping its model's own normalized frame into
the **shared** frame the camera navigates (`ViewportPanel::set_model_transform`).
It is applied to the CAMERA, not to the geometry — `model_c2w()` maps the nav
pose through the inverse before either renderer sees it — so dragging a model
into place is free, and it works identically for the engine renderer and the
GL preview.

By default that similarity puts every model in the FIRST one's frame: two
reconstructions of one scene share world coordinates, so aligning their fitted
frames (`SplatViewer::frame_center` / `frame_unit`) lines them up exactly. That
is also what places a mesh in its own splats' frame on the meshing screen,
which used to be a special case there. A model with nothing to do with the
first one gets "place in the first model's frame" turned off and is framed on
its own; the position / rotation / size controls adjust it by hand either way.

The axes-and-grid overlay is the first model's, in that model's coordinates:
the models in one view are meant to be one scene, so one grid keeps the panes
comparable.

## One mutex

Every pane's `RenderWorker` takes the SAME engine mutex (`CompareView` owns
it). Two panes rendering concurrently would race over which scene is bound,
which is not a subtle failure — it is one model rendered with another's camera.
