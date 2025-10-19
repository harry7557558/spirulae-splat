import torch
import traceback
import math

def _color_text(text: str, color: str):
    if color == None or color == '':
        return text
    color = {
        'k': 30, 'r': 31, 'g': 32, 'y': 33, 'b': 34, 'p': 35, 'c': 36,
    }[color]
    return f"\033[{color}m{text}\033[m"

def check_close(name, a, b, atol=1e-5, rtol=1e-5):
    if a is None:
        print(name, _color_text("None None", 'r'), [*b.shape], str(b.dtype).lstrip("torch."))
        return
    if b is None:
        print(name, [*a.shape], str(a.dtype).lstrip("torch."), None, None, end=' | ')
        b = torch.zeros_like(a)
    else:
        print(name, [*a.shape], str(a.dtype).lstrip("torch."),
            [*b.shape], str(b.dtype).lstrip("torch."), end=' | ')
    if a.shape != b.shape:
        print(_color_text("Shape mismatch", 'r'))
        return
    if a.dtype != b.dtype:
        print(_color_text("dtype mismatch", 'r'))
        return
    abserr = torch.amax(torch.abs(a-b))
    relerr = abserr / torch.fmax(torch.abs(a), torch.abs(b)).clip(min=1e-20).mean()
    abserr_μ = torch.mean(torch.abs(a-b))
    relerr_μ = abserr_μ / torch.fmax(torch.abs(a), torch.abs(b)).clip(min=1e-20).mean()
    fmt = lambda err, tol: _color_text(
        f"{err:.2g}", 'g' if err < 0.1*tol else 'y' if err < tol else 'r')
    n = b.numel()
    xsc = math.log(n) + 0.5/n + 0.577  # E(max of n unit exponential random)
    print(f"diff: μa={fmt(abserr_μ, atol)}, μr={fmt(relerr_μ, rtol)}, xa={fmt(abserr, atol*xsc)}, xr={fmt(relerr, rtol*xsc)}")


def timeit(fun, name: str, repeat=20):
    from time import perf_counter

    for i in range(2):
        fun()
    torch.cuda.synchronize()

    time0 = perf_counter()
    for i in range(repeat):
        fun()
        torch.cuda.synchronize()
    time1 = perf_counter()

    dt = 1e3 * (time1-time0) / repeat
    print(f"{name}: {dt:.2f} ms")
