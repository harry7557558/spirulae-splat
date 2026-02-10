#!/usr/bin/env python3

import numpy as np
from spirulae_splat.viewer.model import SplatModel

from dataclasses import dataclass
from typing import Optional, Tuple

import matplotlib.pyplot as plt


@dataclass
class VignettingChannelParams:
    cx: np.ndarray
    cy: np.ndarray
    alpha0: np.ndarray
    alpha1: np.ndarray
    alpha2: np.ndarray

@dataclass
class ColorPPISPParams:
    b: np.ndarray
    r: np.ndarray
    g: np.ndarray
    n: np.ndarray

@dataclass
class CRFPPISPChannelParams:
    toe: np.ndarray
    shoulder: np.ndarray
    gamma: np.ndarray
    center: np.ndarray

@dataclass
class RQSCRFPPISPChannelParams:
    g0: np.ndarray
    g1: np.ndarray
    x0: np.ndarray
    y0: np.ndarray
    gc: np.ndarray

@dataclass
class PPISPParams:
    """Data class to hold unpacked PPISP parameters."""
    exposure: np.ndarray
    vignetting: Tuple[VignettingChannelParams, VignettingChannelParams, VignettingChannelParams]
    color: ColorPPISPParams
    crf: Optional[Tuple[CRFPPISPChannelParams, CRFPPISPChannelParams, CRFPPISPChannelParams]]
    rqs_crf: Optional[Tuple[RQSCRFPPISPChannelParams, RQSCRFPPISPChannelParams, RQSCRFPPISPChannelParams]]


def unpack_ppisp_params(ppisp_params: np.ndarray) -> PPISPParams:
    """Unpack PPISP parameters into their components."""
    exposure = ppisp_params[:, 0]
    vignetting = []
    for i in range(3):
        vignetting.append(VignettingChannelParams(
            cx=ppisp_params[:, 1 + i*5],
            cy=ppisp_params[:, 2 + i*5],
            alpha0=ppisp_params[:, 3 + i*5],
            alpha1=ppisp_params[:, 4 + i*5],
            alpha2=ppisp_params[:, 5 + i*5],
        ))
    color = ColorPPISPParams(
        b=ppisp_params[:, 16:18],
        r=ppisp_params[:, 18:20],
        g=ppisp_params[:, 20:22],
        n=ppisp_params[:, 22:24],
    )
    crf = None
    rqs_crf = None
    if ppisp_params.shape[1] == 36:
        crf = []
        for i in range(3):
            crf.append(CRFPPISPChannelParams(
                toe=ppisp_params[:, 24 + i*4],
                shoulder=ppisp_params[:, 25 + i*4],
                gamma=ppisp_params[:, 26 + i*4],
                center=ppisp_params[:, 27 + i*4],
            ))
    elif ppisp_params.shape[1] == 39:
        rqs_crf = []
        for i in range(3):
            rqs_crf.append(RQSCRFPPISPChannelParams(
                g0=ppisp_params[:, 24 + i*5],
                g1=ppisp_params[:, 25 + i*5],
                x0=ppisp_params[:, 26 + i*5],
                y0=ppisp_params[:, 27 + i*5],
                gc=ppisp_params[:, 28 + i*5],
            ))
    else:
        raise ValueError(f"Unexpected number of PPISP parameters: {ppisp_params.shape[1]}")
    return PPISPParams(
        exposure=exposure,
        vignetting=tuple(vignetting),
        color=color,
        crf=tuple(crf) if crf is not None else None,
        rqs_crf=tuple(rqs_crf) if rqs_crf is not None else None,
    )


def plot_ppisp_params(ppisp_params: PPISPParams):
    """Plot PPISP parameters
        - Exposure over time, as a line plot
        - Vignetting radial curves, as overlapped line plots, with RGB colors overlapped
        - Vignetting centers, as a scatter plot of cx vs cy, with RGB colors
        - Color correction parameters, as scatter plots on a color triangle
        - CRF parameters as curves, similar to vignetting curves
    """
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))

    # Exposure
    axs[0, 0].plot(ppisp_params.exposure)
    axs[0, 0].set_title("Exposure Offset")
    axs[0, 0].axhline(0.0, color='k', linestyle='--', alpha=0.5)

    # Vignetting curves
    r = np.linspace(0, np.sqrt(0.5), 100)  # Normalized radius from center to corner
    for i in range(len(ppisp_params.vignetting[0].alpha0)):
        for j, color in enumerate(['r', 'g', 'b']):
            vignetting = ppisp_params.vignetting[j]
            alpha0 = vignetting.alpha0[i]
            alpha1 = vignetting.alpha1[i]
            alpha2 = vignetting.alpha2[i]
            vignetting_curve = 1 + alpha0 * r**2 + alpha1 * r**4 + alpha2 * r**6
            axs[0, 1].plot(r, vignetting_curve, color=color, alpha=0.3, linewidth=1)
    axs[0, 1].set_ylim(0.0, 1.2)
    axs[0, 1].set_title("Vignetting Curves")

    # Vignetting centers
    for i in range(len(ppisp_params.vignetting[0].cx)):
        for j, color in enumerate(['r', 'g', 'b']):
            vignetting = ppisp_params.vignetting[j]
            axs[0, 2].scatter(vignetting.cx[i], vignetting.cy[i], color=color, marker='+', label=f'Channel {color.upper()}' if i == 0 else "")
    axs[0, 2].set_xlim(-0.5, 0.5)
    axs[0, 2].set_ylim(-0.5, 0.5)
    axs[0, 2].grid(True, linestyle='--', alpha=0.5)
    axs[0, 2].set_title("Vignetting Centers")
    axs[0, 2].legend()

    # Color correction parameters
    cos_60 = 0.5
    sin_60 = np.sqrt(3) / 2
    axs[1, 0].plot([-sin_60, sin_60, 0, -sin_60], [-cos_60, -cos_60, 1, -cos_60], color='k', label='Ideal Primaries')
    axs[1, 0].scatter(ppisp_params.color.r[:, 0] - sin_60, ppisp_params.color.r[:, 1] - cos_60, color='r', marker='.', label='Red Channel')
    axs[1, 0].scatter(ppisp_params.color.g[:, 0] + sin_60, ppisp_params.color.g[:, 1] - cos_60, color='g', marker='.', label='Green Channel')
    axs[1, 0].scatter(ppisp_params.color.b[:, 0], ppisp_params.color.b[:, 1] + 1, color='b', marker='.', label='Blue Channel')
    axs[1, 0].scatter(ppisp_params.color.n[:, 0], ppisp_params.color.n[:, 1], color='k', marker='.', label='Neutral')
    axs[1, 0].set_aspect('equal')
    axs[1, 0].set_title("Color Correction Parameters")
    axs[1, 0].legend()

    # CRF curves
    def bounded_positive_forward(raw, min_value):
        return min_value + np.log(1.0 + np.exp(raw))
    def clamped_forward(raw):
        return 1.0 / (1.0 + np.exp(-raw))
    if ppisp_params.crf is not None:
        for i in range(len(ppisp_params.crf[0].toe)):
            for j, color in enumerate(['r', 'g', 'b']):
                crf = ppisp_params.crf[j]
                x = np.linspace(0, 1, 100)
                toe = bounded_positive_forward(crf.toe[i], 0.3)
                shoulder = bounded_positive_forward(crf.shoulder[i], 0.3)
                gamma = bounded_positive_forward(crf.gamma[i], 0.1)
                center = clamped_forward(crf.center[i])
                lerp_val = toe * (1.0 - center) + shoulder * center
                a = (shoulder * center) / lerp_val
                b = 1.0 - a
                crf_curve = np.where(x < center, a * (x / center) ** toe, 1.0 - b * ((1.0 - x) / (1.0 - center)) ** shoulder)
                crf_curve = np.clip(crf_curve, 0.0, 1.0) ** gamma
                axs[1, 1].plot(x, crf_curve, color=color, alpha=0.3)
        axs[1, 1].set_title("CRF Curves")
    elif ppisp_params.rqs_crf is not None:
        for i in range(len(ppisp_params.rqs_crf[0].g0)):
            for j, color in enumerate(['r', 'g', 'b']):
                rqs_crf = ppisp_params.rqs_crf[j]
                x = np.linspace(0, 1, 100)
                g0 = np.exp(rqs_crf.g0[i])
                g1 = np.exp(rqs_crf.g1[i])
                x0 = clamped_forward(rqs_crf.x0[i])
                y0 = clamped_forward(rqs_crf.y0[i])
                gc = np.exp(rqs_crf.gc[i])
                rqs_curve = np.where(
                    x < x0,
                    y0 * (y0 / x0 * (x / x0) ** 2 + g0 * (x / x0) * (1 - x / x0)) / (y0 / x0 + (g0 + gc - 2 * y0 / x0) * (x / x0) * (1 - x / x0)),
                    y0 + (1 - y0) * ((1 - y0) / (1 - x0) * ((x - x0) / (1 - x0)) ** 2 + gc * ((x - x0) / (1 - x0)) * (1 - (x - x0) / (1 - x0))) / ((1 - y0) / (1 - x0) + (gc + g1 - 2 * (1 - y0) / (1 - x0)) * ((x - x0) / (1 - x0)) * (1 - (x - x0) / (1 - x0)))
                ) ** gc
                axs[1, 1].plot(x, rqs_curve, color=color, alpha=0.3, linewidth=1)
        axs[1, 1].set_title("CRF Curves (RQS)")

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("model_path", nargs=1, help="Path to the model to visualize")
    args = parser.parse_args()

    model = SplatModel(args.model_path[0])

    if model.ppisp_params is None:
        print("No PPISP parameters found in the model.")
        exit(0)

    ppisp_params = unpack_ppisp_params(model.ppisp_params.detach().cpu().numpy())

    # print(ppisp_params)
    plot_ppisp_params(ppisp_params)
