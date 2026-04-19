from __future__ import annotations

import json
import os
import pathlib
import threading
import webbrowser
from typing import Any, Callable, Dict, List, Optional

from .render_worker import RenderWorker
from .http_server import HTTPThread


class ViewerServer:

    def __init__(
        self,
        render_fn: Callable,
        *,
        http_host: str = "localhost",
        http_port: int = 7007,
        open_browser: bool = False,
    ) -> None:
        self._render_fn = render_fn
        self._http_host = http_host
        self._http_port = http_port
        self._open_browser = open_browser

        self._render_worker: Optional[RenderWorker] = None
        self._http_thread: Optional[HTTPThread] = None

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

    def _build_html(self) -> str:
        template_path = pathlib.Path(__file__).parent / "viewer.html"
        html = template_path.read_text(encoding="utf-8")

        # Inject server-side config into the JS
        html = html.replace("SERVER_HOST", self._http_host)
        html = html.replace("SERVER_PORT", str(self._http_port))

        return html
