"""
http_server.py
~~~~~~~~~~~~~~
Minimal HTTP server that serves the single-page viewer HTML.
Uses only the stdlib ``http.server`` module.
"""

from __future__ import annotations

import threading
from http.server import BaseHTTPRequestHandler, HTTPServer
from typing import Optional


class _Handler(BaseHTTPRequestHandler):
    html_content: bytes = b""

    def do_GET(self) -> None:  # noqa: N802
        if self.path in ("/", "/index.html"):
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(self.html_content)))
            self.end_headers()
            self.wfile.write(self.html_content)
        else:
            self.send_response(404)
            self.end_headers()

    def log_message(self, fmt: str, *args: object) -> None:  # silence default logging
        pass


class HTTPThread:
    """Serves the viewer HTML on a background daemon thread."""

    def __init__(self, html: str, host: str = "localhost", port: int = 8080) -> None:
        _Handler.html_content = html.encode("utf-8")
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
