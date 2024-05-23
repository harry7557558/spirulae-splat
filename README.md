# spirulae-splat
My custom method for Nerfstudio. Modified from the `splatfacto` method in official Nerfstudio implementation.

Major changes / Novelties:

- Use 2D patches instead of 3D, add depth and normal regularization, as per SIGGRAPH 2024 paper [*2D Gaussian Splatting for Geometrically Accurate Radiance Fields*](https://arxiv.org/abs/2403.17888)
  - Rasterization is still based on traditional projection (rather than ray-surface intersection introduced in the paper); singularity is eliminated by replacing division by depth with multiplication, which perserves the unit normal

- Use polynomial spline kernel patches instead of Gaussians, as inspired by [kernels in computational physics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics)
  - This increases rendering speed by zeroing out opacities beyond a certain radius
  - A patch with higher overall opacity produces fewer Gaussians, especially when combined with removing patches with absolute gradient below a threshold
  - $\alpha=\max(1-r^2,0)$ empirically produces fewer Gaussians, faster convergence, and higher visual and geometric fidelity compared to kernels with higher continuities like commonly used SPH kernels

- Initialize patch scale and orientation based on SVD of neighbor SfM points

## Registering with Nerfstudio
Ensure that nerfstudio has been installed according to the [instructions](https://docs.nerf.studio/quickstart/installation.html). Clone or fork this repository and run the commands:

```
cd spirulae-splat/
git submodule update --init
conda activate nerfstudio
pip install -e .
ns-install-cli
```

## Running the new method
This repository creates a new Nerfstudio method named "spirulae". To train with it, run the command:
```
ns-train spirulae --data [PATH]
```