from __future__ import annotations

import io
import queue
import threading
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Optional, Tuple

import numpy as np


@dataclass
class RenderRequest:
    c2w: Any                      # 3×4 numpy array
    fx: float
    fy: float
    cx: float
    cy: float
    width: int
    height: int
    camera_model: str             # "PINHOLE" | "FISHEYE"
    request_id: int = 0


@dataclass
class RenderResult:
    request_id: int
    buffers: Dict[str, Any]       # dict[str, torch.Tensor] as returned by the render fn
    error: Optional[str] = None


class RenderWorker:

    def __init__(
        self,
        render_fn: Callable,
    ) -> None:
        self._render_fn = render_fn

        # Latest pending request (None = idle)
        self._pending: Optional[RenderRequest] = None
        self._pending_lock = threading.Lock()
        self._work_event = threading.Event()

        # Completed results, consumed by the WebSocket sender
        self._result_queue: queue.Queue[RenderResult] = queue.Queue(maxsize=4)

        self._running = False
        self._thread: Optional[threading.Thread] = None
        self._request_counter = 0

    def start(self) -> None:
        self._running = True
        self._thread = threading.Thread(target=self._loop, daemon=True, name="RenderWorker")
        self._thread.start()

    def stop(self) -> None:
        self._running = False
        self._work_event.set()
        if self._thread:
            self._thread.join(timeout=5)

    def submit(self, req: RenderRequest) -> None:
        """Replace the pending request with a newer one (latest-wins)."""
        self._request_counter += 1
        req.request_id = self._request_counter
        with self._pending_lock:
            self._pending = req
        self._work_event.set()

    def get_result(self, timeout: float = 0.05) -> Optional[RenderResult]:
        try:
            return self._result_queue.get(timeout=timeout)
        except queue.Empty:
            return None

    def _loop(self) -> None:
        while self._running:
            self._work_event.wait()
            self._work_event.clear()

            with self._pending_lock:
                req = self._pending
                self._pending = None

            if req is None:
                continue

            try:
                buffers = self._render_fn(
                    req.c2w,
                    req.fx, req.fy, req.cx, req.cy,
                    req.width, req.height,
                    req.camera_model,
                )
                result = RenderResult(request_id=req.request_id, buffers=buffers)
            except BrokenPipeError as e:
                result = RenderResult(request_id=req.request_id, buffers={}, error=str(exc))
            except Exception as exc:
                import traceback
                traceback.print_exc()
                result = RenderResult(request_id=req.request_id, buffers={}, error=str(exc))

            # Non-blocking put; drop if sender is backed up
            try:
                self._result_queue.put_nowait(result)
            except queue.Full:
                # Drain one stale result and retry
                try:
                    self._result_queue.get_nowait()
                except queue.Empty:
                    pass
                try:
                    self._result_queue.put_nowait(result)
                except queue.Full:
                    pass


def encode_buffer_to_jpeg(
    tensor: Any,
    post_processor: Callable,
    quality: int,
    post_process_params: dict = {}
) -> bytes:
    import cv2

    tensor = post_processor(tensor, **post_process_params).cpu().numpy()

    success, buf = cv2.imencode(
        '.jpg', cv2.cvtColor(tensor, cv2.COLOR_RGB2BGR),
        [int(cv2.IMWRITE_JPEG_QUALITY), quality]
    )
    if success:
        return buf.tobytes()
