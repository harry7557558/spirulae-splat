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
    buffer_key: str = "rgb"       # which channel the viewer wants annotated
    show_training_cameras: bool = False
    show_grid: bool = False
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
        on_submit: Optional[Callable[[], None]] = None,
        on_idle: Optional[Callable[[], None]] = None,
    ) -> None:
        self._render_fn = render_fn
        # Called synchronously from the HTTP handler thread when submit() runs.
        # Used by the trainer to flip a "render desired" flag *immediately*
        # (before the worker thread has a chance to call render_fn) so the
        # very next training-iteration boundary can yield the lock to render.
        self._on_submit = on_submit
        # Called from the worker thread once it finishes a render AND has no
        # more pending work. Pairs with on_submit so the trainer's flag stays
        # set across back-to-back submissions (which can otherwise let
        # training sneak in between).
        self._on_idle = on_idle

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
        if self._on_submit is not None:
            self._on_submit()
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
                    req.buffer_key,
                    show_training_cameras=req.show_training_cameras,
                    show_grid=req.show_grid,
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

            # Clear the "render desired" flag only if no follow-up request
            # arrived while we were busy. Otherwise leave it set so training
            # keeps yielding while the worker chews through the backlog.
            if self._on_idle is not None:
                with self._pending_lock:
                    has_more = self._pending is not None
                if not has_more:
                    self._on_idle()


def encode_buffer_to_jpeg(
    tensor: Any,
    post_processor: Optional[Callable],
    quality: int,
    post_process_params: dict = {}
) -> bytes:
    import cv2

    # Trainer.render now invokes the post-processor + D->H copy itself under
    # the trainer lock, so by the time we get here `tensor` is already a CPU
    # numpy array. The legacy closure path (post_processor not None) is kept
    # as a fallback in case other render fns still return GPU tensors.
    if not isinstance(tensor, np.ndarray):
        if post_processor is not None:
            tensor = post_processor(tensor, **post_process_params).cpu().numpy()
        else:
            tensor = tensor.cpu().numpy()

    success, buf = cv2.imencode(
        '.jpg', cv2.cvtColor(tensor, cv2.COLOR_RGB2BGR),
        [int(cv2.IMWRITE_JPEG_QUALITY), quality]
    )
    if success:
        return buf.tobytes()
