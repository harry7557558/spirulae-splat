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

                if buffer_key not in result.buffers:
                    self.send_response(400)
                    self.end_headers()
                    return

                # Update last_keys
                _Handler.last_keys = list(result.buffers.keys())

                # Encode to JPEG
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

    def _default_render_for_keys(self) -> None:
        # Default params
        c2w = np.eye(3, 4, dtype=np.float32)
        req = RenderRequest(
            c2w=c2w,
            fx=500, fy=500, cx=256, cy=256,
            width=512, height=512,
            camera_model="PINHOLE",
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

    def __init__(self, html: str, render_worker: RenderWorker, progress_fn: Optional[Callable], host: str = "0.0.0.0", port: int = 8080) -> None:
        _Handler.html_content = html.encode("utf-8")
        _Handler.render_worker = render_worker
        _Handler.progress_fn = progress_fn
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
