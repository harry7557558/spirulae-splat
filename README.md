# spirulae-splat
My custom method for Nerfstudio. Modified from the `splatfacto` method in official Nerfstudio implementation.

### Major changes / Novelties:

- Use 2D splats instead of 3D, add depth and normal regularization, use ray-splat intersection, as per SIGGRAPH 2024 paper [*2D Gaussian Splatting for Geometrically Accurate Radiance Fields*](https://arxiv.org/abs/2403.17888)
  - For normal regularization, singularity is eliminated by replacing division by depth with multiplication, which perserves the unit normal
  - For depth regularization, use L1 loss based on center of splat instead of L2 loss; This could be improved by figuring out fast per-pixel sorting

- Use polynomial spline kernel splats instead of Gaussians, as inspired by [kernels used in computational physics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics)
  - This increases rendering speed by zeroing out opacities beyond a certain radius
  - A splat with higher overall opacity produces fewer Gaussians, especially when combined with culling splats with absolute position gradient below a threshold
  - $\alpha=\max(1-r^2,0)$ empirically produces fewer Gaussians, faster convergence, and higher visual and geometric fidelity compared to kernels with higher continuities like commonly used SPH kernels

- Use [cylindrical harmonics](https://en.wikipedia.org/wiki/Cylindrical_harmonics) to model uv-dependent color
  - Base+SH color multiplied by the sigmoid of 2D "polar" harmonics
  - Able to better capture fine color details and potentially reduce export file size on complex scenes
  - Slows down training by about a half, due to resource spent on Bessel function evaluation and loading CH coefficients from global memory

- Initialize splat scale and orientation based on SVD of neighbor SfM points

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