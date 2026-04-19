"""
server.py
~~~~~~~~~
ViewerServer: the main entry point.

Usage::

    from renderer3d import ViewerServer

    def my_render_fn(c2w, fx, fy, cx, cy, w, h, camera_model):
        ...
        return {"rgb": tensor, "depth": depth_tensor}

    server = ViewerServer(my_render_fn)
    server.start()          # opens browser automatically
    server.wait()           # blocks until Ctrl-C
"""

from __future__ import annotations

import json
import os
import pathlib
import threading
import webbrowser
from typing import Any, Callable, Dict, List, Optional

from .render_worker import RenderWorker
from .http_server import HTTPThread


# ── Slider / dropdown descriptor types ────────────────────────────────────────

class SliderDef:
    """Descriptor for an extra slider in the control panel."""
    def __init__(
        self,
        id: str,
        label: str,
        min: float,
        max: float,
        step: float = 1.0,
        value: float = 0.0,
        unit: str = "",
    ) -> None:
        self.id = id
        self.label = label
        self.min = min
        self.max = max
        self.step = step
        self.value = value
        self.unit = unit

    def to_dict(self) -> dict:
        return dict(id=self.id, label=self.label, min=self.min, max=self.max,
                    step=self.step, value=self.value, unit=self.unit)


class DropdownDef:
    """Descriptor for an extra dropdown in the control panel."""
    def __init__(
        self,
        id: str,
        label: str,
        options: List[Dict[str, str]],  # [{"value": ..., "label": ...}]
        default: str = "",
    ) -> None:
        self.id = id
        self.label = label
        self.options = options
        self.default = default or (options[0]["value"] if options else "")

    def to_dict(self) -> dict:
        return dict(id=self.id, label=self.label, options=self.options, default=self.default)


class ViewerServer:
    """
    Lightweight browser-based 3D viewer server.

    Parameters
    ----------
    render_fn:
        Callable with signature::

            render_fn(
                c2w: np.ndarray,          # (3, 4) float32
                fx: float, fy: float,
                cx: float, cy: float,
                width: int, height: int,
                camera_model: str,        # "PINHOLE" | "FISHEYE"
            ) -> dict[str, torch.Tensor]  # float32, shape (1, H, W, C)

        The function is always called from a single dedicated thread.

    http_host / http_port:
        Address for the HTTP server.

    extra_sliders / extra_dropdowns:
        Additional UI controls that appear in the "Visualize" panel section.
        Their current values are forwarded to the render function via the
        ``extra`` field in the WebSocket message (access via the ``extra``
        kwarg if you extend the protocol, or ignore them in the default
        render_fn signature).

    open_browser:
        Automatically open the viewer in the default browser (default True).
    """

    def __init__(
        self,
        render_fn: Callable,
        *,
        http_host: str = "localhost",
        http_port: int = 7007,
        extra_sliders: Optional[List[SliderDef]] = None,
        extra_dropdowns: Optional[List[DropdownDef]] = None,
        open_browser: bool = True,
    ) -> None:
        self._render_fn = render_fn
        self._http_host = http_host
        self._http_port = http_port
        self._extra_sliders = extra_sliders or []
        self._extra_dropdowns = extra_dropdowns or []
        self._open_browser = open_browser

        self._render_worker: Optional[RenderWorker] = None
        self._http_thread: Optional[HTTPThread] = None

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def start(self) -> "ViewerServer":
        """Start all background threads and (optionally) open the browser."""
        # 1. Render worker
        self._render_worker = RenderWorker(
            render_fn=self._render_fn,
        )
        self._render_worker.start()

        # 2. HTTP server
        html = self._build_html()
        self._http_thread = HTTPThread(
            html=html,
            render_worker=self._render_worker,
            host=self._http_host,
            port=self._http_port,
        )
        self._http_thread.start()

        url = f"http://{self._http_host}:{self._http_port}/"
        print(f"[renderer3d] Viewer at  {url}")

        if self._open_browser:
            threading.Timer(0.5, webbrowser.open, args=[url]).start()

        return self

    def stop(self) -> None:
        """Stop all background threads."""
        if self._render_worker:
            self._render_worker.stop()
        if self._http_thread:
            self._http_thread.stop()

    def wait(self) -> None:
        """Block until a KeyboardInterrupt (Ctrl-C)."""
        try:
            threading.Event().wait()
        except KeyboardInterrupt:
            print("\n[renderer3d] Shutting down…")
            self.stop()

    # Allow use as a context manager
    def __enter__(self) -> "ViewerServer":
        return self.start()

    def __exit__(self, *_: Any) -> None:
        self.stop()

    # ------------------------------------------------------------------
    # HTML builder
    # ------------------------------------------------------------------

    def _build_html(self) -> str:
        template_path = pathlib.Path(__file__).parent / "viewer.html"
        html = template_path.read_text(encoding="utf-8")

        extra_sliders_json   = json.dumps([s.to_dict() for s in self._extra_sliders])
        extra_dropdowns_json = json.dumps([d.to_dict() for d in self._extra_dropdowns])

        # Inject server-side config into the JS
        html = html.replace("SERVER_HOST", self._http_host)
        html = html.replace("SERVER_PORT", str(self._http_port))
        html = html.replace("EXTRA_SLIDERS_JSON", extra_sliders_json)
        html = html.replace("EXTRA_DROPDOWNS_JSON", extra_dropdowns_json)

        return html
