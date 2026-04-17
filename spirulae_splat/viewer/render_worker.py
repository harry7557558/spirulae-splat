"""
render_worker.py
~~~~~~~~~~~~~~~~
Runs the user-supplied render function in a single dedicated thread so the
main thread is never blocked and thread-safety is preserved (the render
function is always called from this one thread, never concurrently).
"""

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
    """One render job submitted to the worker."""
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
    """Result produced by the worker."""
    request_id: int
    buffers: Dict[str, Any]       # dict[str, torch.Tensor] as returned by the render fn
    error: Optional[str] = None


class RenderWorker:
    """
    Wraps a render function and ensures it is always called from one thread.

    The worker keeps only the *latest* pending request; older ones are
    discarded so the viewer stays responsive under high load.
    """

    def __init__(
        self,
        render_fn: Callable,
        jpeg_quality: int = 85,
    ) -> None:
        self._render_fn = render_fn
        self.jpeg_quality = jpeg_quality

        # Latest pending request (None = idle)
        self._pending: Optional[RenderRequest] = None
        self._pending_lock = threading.Lock()
        self._work_event = threading.Event()

        # Completed results, consumed by the WebSocket sender
        self._result_queue: queue.Queue[RenderResult] = queue.Queue(maxsize=4)

        self._running = False
        self._thread: Optional[threading.Thread] = None
        self._request_counter = 0

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

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

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

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
                for key, value in buffers.items():
                    if len(value.shape) == 3:
                        buffers[key] = value[None]
                result = RenderResult(request_id=req.request_id, buffers=buffers)
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


# ------------------------------------------------------------------
# Image encoding helpers
# ------------------------------------------------------------------

def tensor_to_numpy(tensor: Any) -> np.ndarray:
    """Convert a (1, H, W, C) float32 torch tensor to uint8 numpy HWC."""
    import torch  # noqa: F401 — only needed at encode time
    if hasattr(tensor, "detach"):
        arr = tensor.detach().cpu().float().numpy()
    else:
        arr = np.asarray(tensor, dtype=np.float32)
    arr = arr[0]  # (H, W, C)
    return arr


def encode_buffer_to_jpeg(
    tensor: Any,
    quality: int = 85,
    colormap: str = "turbo",
) -> bytes:
    """
    Encode a rendered buffer tensor to JPEG bytes.

    * C==3 → clamp [0,1], convert to uint8 RGB
    * C==1 → apply matplotlib colormap, convert to uint8 RGB
    """
    import cv2

    arr = tensor_to_numpy(tensor)  # (H, W, C)
    c = arr.shape[2]

    if c == 3:
        rgb = np.clip(arr, 0.0, 1.0)
        rgb_u8 = (rgb * 255).astype(np.uint8)
    elif c == 1:
        mono = arr[:, :, 0]
        # Normalise ignoring NaN/inf
        valid = mono[np.isfinite(mono)]
        if valid.size > 0:
            lo, hi = valid.min(), valid.max()
            span = hi - lo
            if span > 1e-8:
                mono = (mono - lo) / span
            else:
                mono = np.zeros_like(mono)
        else:
            mono = np.zeros_like(mono)
        mono = np.clip(mono, 0.0, 1.0)
        try:
            import matplotlib.cm as cm
            cmap = cm.get_cmap(colormap)
            rgb_u8 = (cmap(mono)[:, :, :3] * 255).astype(np.uint8)
        except Exception:
            mono_u8 = (mono * 255).astype(np.uint8)
            rgb_u8 = np.stack([mono_u8] * 3, axis=2)
    else:
        raise ValueError(f"Unsupported channel count: {c}")

    success, buf = cv2.imencode(
        '.jpg', cv2.cvtColor(rgb_u8, cv2.COLOR_RGB2BGR),
        [int(cv2.IMWRITE_JPEG_QUALITY), quality]
    )
    if success:
        return buf.tobytes()
