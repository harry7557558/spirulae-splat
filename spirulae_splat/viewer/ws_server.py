"""
ws_server.py
~~~~~~~~~~~~
Pure-asyncio WebSocket server (using the ``websockets`` library) that:

  1. Receives camera / control messages from the browser (JSON).
  2. Forwards render requests to the RenderWorker.
  3. Sends encoded JPEG frames back over the same WebSocket connection.

A separate asyncio thread runs the event loop so it never blocks the
caller's main thread.
"""

from __future__ import annotations

import asyncio
import json
import queue
import threading
import time
from typing import Any, Callable, Dict, Optional, Set

import numpy as np

from .render_worker import RenderRequest, RenderResult, RenderWorker, encode_buffer_to_jpeg


class _WSServer:
    """Internal asyncio-based WebSocket server."""

    def __init__(
        self,
        render_worker: RenderWorker,
        host: str,
        port: int,
        jpeg_quality: int,
    ) -> None:
        self._worker = render_worker
        self._host = host
        self._port = port
        self._jpeg_quality = jpeg_quality

        self._connections: Set[Any] = set()
        self._connections_lock = asyncio.Lock()

        # Latest known buffer keys (populated after first render)
        self._last_buffer_keys: list[str] = ["rgb"]

        self._loop: Optional[asyncio.AbstractEventLoop] = None
        self._server = None

    # ------------------------------------------------------------------
    # Start / stop (called from the background thread)
    # ------------------------------------------------------------------

    async def run(self) -> None:
        import websockets  # type: ignore

        self._loop = asyncio.get_running_loop()
        async with websockets.serve(
            self._handle_client,
            self._host,
            self._port,
            max_size=2 * 1024 * 1024,
            ping_interval=20,
            ping_timeout=20,
        ) as server:
            self._server = server
            await asyncio.Future()  # run forever

    # ------------------------------------------------------------------
    # Per-client handler
    # ------------------------------------------------------------------

    async def _handle_client(self, websocket: Any) -> None:
        async with self._connections_lock:
            self._connections.add(websocket)

        sender_task = asyncio.create_task(self._sender(websocket))
        try:
            async for raw in websocket:
                await self._on_message(websocket, raw)
        except Exception:
            pass
        finally:
            sender_task.cancel()
            async with self._connections_lock:
                self._connections.discard(websocket)

    async def _on_message(self, websocket: Any, raw: str | bytes) -> None:
        """Handle an incoming message from the browser."""
        if isinstance(raw, bytes):
            return
        try:
            msg = json.loads(raw)
        except json.JSONDecodeError:
            return

        mtype = msg.get("type")

        if mtype == "render_request":
            c2w_flat = msg["c2w"]          # 12-element list (row-major 3×4)
            c2w = np.array(c2w_flat, dtype=np.float32).reshape(3, 4)
            req = RenderRequest(
                c2w=c2w,
                fx=float(msg["fx"]),
                fy=float(msg["fy"]),
                cx=float(msg["cx"]),
                cy=float(msg["cy"]),
                width=int(msg["width"]),
                height=int(msg["height"]),
                camera_model=msg.get("camera_model", "PINHOLE"),
            )
            self._worker.submit(req)

    # ------------------------------------------------------------------
    # Result sender loop (one per connection)
    # ------------------------------------------------------------------

    async def _sender(self, websocket: Any) -> None:
        """Continuously pull results from the worker and send to this client."""
        loop = asyncio.get_running_loop()
        while True:
            # Poll the render worker result queue in a thread-pool executor
            # so we don't block the event loop.
            result: Optional[RenderResult] = await loop.run_in_executor(
                None, self._worker.get_result, 0.05
            )
            if result is None:
                await asyncio.sleep(0)
                continue

            if result.error:
                err_msg = json.dumps({"type": "error", "message": result.error})
                try:
                    await websocket.send(err_msg)
                except Exception:
                    break
                continue

            # Update known keys
            keys = list(result.buffers.keys())
            if keys != self._last_buffer_keys:
                self._last_buffer_keys = keys
                keys_msg = json.dumps({"type": "buffer_keys", "keys": keys})
                try:
                    await websocket.send(keys_msg)
                except Exception:
                    break

            # Encode and send each buffer as a JSON envelope + binary JPEG
            for key, tensor in result.buffers.items():
                try:
                    jpeg_bytes = await loop.run_in_executor(
                        None,
                        encode_buffer_to_jpeg,
                        tensor,
                        self._jpeg_quality,
                    )
                except Exception as exc:
                    import traceback
                    traceback.print_exc()
                    err_msg = json.dumps({"type": "encode_error", "key": key, "message": str(exc)})
                    try:
                        await websocket.send(err_msg)
                    except Exception:
                        break
                    continue

                # Send metadata frame, then binary frame
                meta = json.dumps({
                    "type": "frame",
                    "key": key,
                    "request_id": result.request_id,
                })
                try:
                    await websocket.send(meta)
                    await websocket.send(jpeg_bytes)
                except Exception:
                    return


class WebSocketThread:
    """
    Runs the asyncio WebSocket server in a daemon background thread.
    """

    def __init__(
        self,
        render_worker: RenderWorker,
        host: str = "localhost",
        port: int = 8765,
        jpeg_quality: int = 85,
    ) -> None:
        self._ws_server = _WSServer(render_worker, host, port, jpeg_quality)
        self._thread: Optional[threading.Thread] = None
        self._loop: Optional[asyncio.AbstractEventLoop] = None

    def start(self) -> None:
        ready = threading.Event()

        def _run() -> None:
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
            self._loop = loop
            ready.set()
            loop.run_until_complete(self._ws_server.run())

        self._thread = threading.Thread(target=_run, daemon=True, name="WSServer")
        self._thread.start()
        ready.wait(timeout=5)

    def stop(self) -> None:
        if self._loop:
            self._loop.call_soon_threadsafe(self._loop.stop)
