from __future__ import annotations

import json
import threading
import time
from http.server import BaseHTTPRequestHandler, HTTPServer
from typing import Optional, Callable
from urllib.parse import parse_qs, urlparse

import numpy as np

from .render_worker import RenderRequest, RenderWorker, encode_buffer_to_jpeg


class _Handler(BaseHTTPRequestHandler):
    html_content: bytes = b""
    render_worker: Optional[RenderWorker] = None
    progress_fn: Optional[Callable] = None
    pause_toggle_fn: Optional[Callable] = None
    last_keys: list[str] = []

    def do_GET(self) -> None:  # noqa: N802
        parsed = urlparse(self.path)
        path = parsed.path
        query = parse_qs(parsed.query)

        if path in ("/", "/index.html"):
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(self.html_content)))
            self.end_headers()
            self.wfile.write(self.html_content)
        elif path == "/render":
            self._handle_render(query)
        elif path == "/buffers":
            self._handle_buffers()
        elif path == "/progress":
            self._handle_progress()
        elif path == "/pause-toggle":
            self._handle_pause_toggle()
        else:
            self.send_response(404)
            self.end_headers()

    def _handle_render(self, query: dict) -> None:
        try:
            buffer_key = query.get("buffer_key", ["rgb"])[0]
            c2w_str = query.get("c2w", [""])[0]
            c2w_flat = [float(x) for x in c2w_str.split(",")]
            c2w = np.array(c2w_flat, dtype=np.float32).reshape(3, 4)
            fx = float(query.get("fx", ["500"])[0])
            fy = float(query.get("fy", ["500"])[0])
            cx = float(query.get("cx", ["256"])[0])
            cy = float(query.get("cy", ["256"])[0])
            width = int(query.get("width", ["512"])[0])
            height = int(query.get("height", ["512"])[0])
            camera_model = query.get("camera_model", ["PINHOLE"])[0]
            jpeg_quality = int(query.get("jpeg_quality", [75])[0])
            show_training_cameras = query.get("show_training_cameras", ["0"])[0].lower() in ("1", "true", "yes")

            MAX_DIM = 2160
            if max(width, height) > 2160:  # prevent OOM
                raise ValueError(f"Image too large (got {width} x {height}, max {MAX_DIM})")

            req = RenderRequest(
                c2w=c2w,
                fx=fx, fy=fy, cx=cx, cy=cy,
                width=width, height=height,
                camera_model=camera_model,
                buffer_key=buffer_key,
            )
            post_process_params = {
                "show_training_cameras": show_training_cameras,
            }
            if self.render_worker:
                self.render_worker.submit(req)
                # Wait for result
                start_time = time.time()
                while time.time() - start_time < 10:  # timeout 10s
                    result = self.render_worker.get_result(0.1)
                    if result and result.request_id == req.request_id:
                        break
                else:
                    self.send_response(500)
                    self.end_headers()
                    return

                if result.error:
                    self.send_response(500)
                    self.end_headers()
                    return

                # Trainer._render inserts None placeholders for every supported
                # viewer buffer so the dropdown can enumerate them all from a
                # single render. The requested buffer_key is always actually
                # rendered (see Trainer._render's want_keys logic), so a None
                # value here means either an unsupported buffer or a request
                # that didn't trigger rendering -- both 400s.
                if buffer_key not in result.buffers \
                        or result.buffers[buffer_key] is None:
                    self.send_response(400)
                    self.end_headers()
                    return

                # Update last_keys (keep placeholder keys so the dropdown
                # shows them; the encode path filters them out below).
                _Handler.last_keys = list(result.buffers.keys())

                # Encode to JPEG. result.buffers[buffer_key] is already an
                # annotated CUDA uint8 [H,W,3] tensor (engine_blit_view); the
                # post-processor closure threads through show_training_cameras.
                jpeg_bytes = encode_buffer_to_jpeg(
                    result.buffers[buffer_key],
                    result.buffers.get("_post_processor", None),
                    jpeg_quality,
                    post_process_params,
                )
                if '_post_processor' in result.buffers:
                    del result.buffers['_post_processor']

                self.send_response(200)
                self.send_header("Content-Type", "image/jpeg")
                self.send_header("Content-Length", str(len(jpeg_bytes)))
                self.end_headers()
                self.wfile.write(jpeg_bytes)
            else:
                self.send_response(500)
                self.end_headers()
        except BrokenPipeError as e:
            self.send_response(400)
            self.end_headers()
        except Exception as e:
            import traceback
            traceback.print_exc()
            self.send_response(400)
            self.end_headers()

    def _handle_buffers(self) -> None:
        if not _Handler.last_keys:
            # Do a default render to get keys
            self._default_render_for_keys()
        keys = [key for key in _Handler.last_keys if not key.startswith('_')]
    
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.end_headers()
        self.wfile.write(json.dumps(keys).encode("utf-8"))

    def _handle_progress(self) -> None:
        if self.progress_fn:
            progress = self.progress_fn()
        else:
            progress = {
                "step": 0,
                "total_steps": 0,
                "elapsed_time": 0,
                "eta": None,
                "latency_ms": None,
            }
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.end_headers()
        self.wfile.write(json.dumps(progress).encode("utf-8"))

    def _handle_pause_toggle(self) -> None:
        if self.pause_toggle_fn:
            paused = self.pause_toggle_fn()
            response = {"paused": paused}
        else:
            response = {"paused": False, "error": "pause_toggle_fn not available"}
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.end_headers()
        self.wfile.write(json.dumps(response).encode("utf-8"))

    def _default_render_for_keys(self) -> None:
        # Default params. Empty buffer_key signals to Trainer._render that the
        # caller wants the FULL output dict (i.e. don't gate via want_keys),
        # so the dropdown lists every available buffer including the
        # debug-only ones (sh, refinement_score, depth_normal).
        c2w = np.eye(3, 4, dtype=np.float32)
        req = RenderRequest(
            c2w=c2w,
            fx=500, fy=500, cx=256, cy=256,
            width=512, height=512,
            camera_model="PINHOLE",
            buffer_key="",
        )
        if self.render_worker:
            self.render_worker.submit(req)
            start_time = time.time()
            while time.time() - start_time < 10:
                result = self.render_worker.get_result(0.1)
                if result and result.request_id == req.request_id:
                    if not result.error:
                        _Handler.last_keys = list(result.buffers.keys())
                    break

    def log_message(self, fmt: str, *args: object) -> None:  # silence default logging
        pass


class HTTPThread:
    """Serves the viewer HTML and handles requests on a background daemon thread."""

    def __init__(self, html: str, render_worker: RenderWorker, progress_fn: Optional[Callable], pause_toggle_fn: Optional[Callable] = None, host: str = "0.0.0.0", port: int = 8080) -> None:
        _Handler.html_content = html.encode("utf-8")
        _Handler.render_worker = render_worker
        _Handler.progress_fn = progress_fn
        _Handler.pause_toggle_fn = pause_toggle_fn
        self._server = HTTPServer((host, port), _Handler)
        self._thread: Optional[threading.Thread] = None

    def start(self) -> None:
        self._thread = threading.Thread(
            target=self._server.serve_forever,
            daemon=True,
            name="HTTPServer",
        )
        self._thread.start()

    def stop(self) -> None:
        self._server.shutdown()
