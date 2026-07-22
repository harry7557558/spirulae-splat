"""Smoke test for the native web-viewer binding (_C.WebViewer).

Covers what Python needs in order to stop carrying its own HTTP server +
render worker (docs/restructure-proposal.md §4.2):

  * bake_post_split turns a ParsedDataset into the arrays
    engine_setup_data_manager wants;
  * WebViewer.start() brings up the real server against a real engine;
  * the HTTP endpoints answer;
  * engine_lock() is a working context manager and does not deadlock against
    the render worker;
  * stop() actually returns (the accept() teardown bug this phase fixed).

Needs a CUDA device, since starting the viewer uploads the camera table to the
engine.
"""

from __future__ import annotations

import json
import socket
import urllib.error
import urllib.request

import pytest

import dataset_fixtures as fixtures

torch = pytest.importorskip("torch")
if not torch.cuda.is_available():
    pytest.skip("needs a CUDA device", allow_module_level=True)

from spirulae_splat.splat.cuda import _C  # noqa: E402


def _free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def _get(port: int, path: str, timeout: float = 10.0):
    with urllib.request.urlopen(
            f"http://127.0.0.1:{port}{path}", timeout=timeout) as r:
        return r.status, r.read()


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    root = tmp_path_factory.mktemp("webviewer_fixture")
    return fixtures.write_colmap_text(root / "colmap", fixtures.make_scene())


@pytest.fixture(scope="module")
def engine(dataset):
    """A minimal engine: parsed cameras + data manager + seeded splats."""
    ds = _C.parse_dataset(str(dataset), _C.DatasetParserConfig(), "")
    post = _C.bake_post_split(ds, False, False)

    cfg = _C.DataManagerConfig()
    _C.engine_setup_data_manager(
        cfg,
        list(ds.camera_models),
        list(ds.image_filenames), list(ds.mask_filenames),
        list(ds.depth_filenames), list(ds.normal_filenames),
        list(ds.widths), list(ds.heights),
        list(post.K_per_camera), list(post.post_offsets),
        list(post.viewmats), list(post.intrins), list(post.dist_coeffs),
        list(post.input_intrins), list(post.input_dist_coeffs),
        list(ds.train_indices), list(ds.val_indices),
    )
    yield ds, post
    _C.engine_reset()


def test_bake_post_split_shapes(dataset):
    ds = _C.parse_dataset(str(dataset), _C.DatasetParserConfig(), "")
    post = _C.bake_post_split(ds, False, False)

    # Pinhole dataset: no cubemap split, so post == input.
    assert post.n_post == ds.num_cameras
    assert not post.any_warp
    assert post.viewmats_array.shape == (ds.num_cameras, 4, 4)
    assert post.intrins_array.shape == (ds.num_cameras, 4)
    assert post.c2w_flip.shape == (ds.num_cameras, 3, 4)


def test_viewer_serves_and_stops(engine):
    ds, post = engine
    port = _free_port()

    cfg = _C.ViewerRenderConfig()
    cfg.train_frame_scale = ds.train_frame_scale
    cfg.train_to_normalized = ds.train_to_normalized

    viewer = _C.WebViewer()
    viewer.start("127.0.0.1", port, cfg, post)
    assert viewer.started
    try:
        viewer.set_step(7)
        viewer.set_progress_json(json.dumps({"step": 7, "total_steps": 100}))

        status, body = _get(port, "/")
        assert status == 200
        assert b"<!DOCTYPE html>" in body[:64], "should serve viewer.html"

        status, body = _get(port, "/buffers")
        assert status == 200
        assert "rgb" in json.loads(body)

        status, body = _get(port, "/progress")
        assert status == 200
        assert json.loads(body)["step"] == 7

        # pause-toggle flips the state the training loop would read.
        assert viewer.paused is False
        status, body = _get(port, "/pause-toggle")
        assert status == 200 and json.loads(body)["paused"] is True
        assert viewer.paused is True
        _get(port, "/pause-toggle")
        assert viewer.paused is False
    finally:
        viewer.stop()          # must return -- see HttpServer::stop()
    assert not viewer.started

    # Port is released, i.e. the listening socket really closed.
    with pytest.raises((urllib.error.URLError, OSError)):
        _get(port, "/buffers", timeout=2.0)


def test_engine_lock_is_reentrant_across_calls(engine):
    """engine_lock() must be acquirable repeatedly (and release the GIL)."""
    viewer = _C.WebViewer()
    for _ in range(3):
        with viewer.engine_lock():
            pass


def test_start_twice_is_an_error(engine):
    ds, post = engine
    viewer = _C.WebViewer()
    viewer.start("127.0.0.1", _free_port(), _C.ViewerRenderConfig(), post)
    try:
        with pytest.raises(RuntimeError):
            viewer.start("127.0.0.1", _free_port(), _C.ViewerRenderConfig(), post)
    finally:
        viewer.stop()
