# spirulae-splat
My custom method for Nerfstudio. Modified from the `splatfacto` method in official Nerfstudio implementation.

### Major changes / Novelties:

- Use 2D splats instead of 3D, add depth and normal regularization, use ray-splat intersection, as per SIGGRAPH 2024 paper [*2D Gaussian Splatting for Geometrically Accurate Radiance Fields*](https://arxiv.org/abs/2403.17888)
  - Use an interpolation function for median depth to reduce discontinuities and encourage denser gradient
  - For normal regularization, singularity is eliminated by replacing division by depth with multiplication, which perserves the unit normal

- Use polynomial spline kernel splats instead of Gaussians, as inspired by [kernels used in computational physics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics)
  - This increases rendering speed by zeroing out opacities beyond a certain radius
  - A splat with higher overall opacity produces fewer Gaussians, especially when combined with culling splats with absolute position gradient below a threshold
  - $\alpha=\max(1-r^2,0)$ empirically produces 1.5x-2x faster training and inference compared to Gaussians, with 0.2-0.5 drop in PSNR compared to Gaussian kernel for a similar number of splats

- Use [cylindrical harmonics](https://en.wikipedia.org/wiki/Cylindrical_harmonics) to model uv-dependent color
  - Base+SH color multiplied by the sigmoid of 2D "polar" harmonics
  - Able to better capture fine color details and potentially reduce export file size on scenes with fine detail
  - Use absolute gradient of CH coefficients as a criteria to split and duplicate splats
  - Slows down training by about a half, due to resource spent on Bessel function evaluation and loading CH coefficients from global memory

- Support [Markov Chain Monte Carlo (MCMC)](https://arxiv.org/abs/2404.09591) for adaptive control of splats, with relocation designed for polynomial splats

- Introduce splat anisotropy by multiplying opacity by a smoothstep of UV directional dot product $1-\mathrm{smoothstep}(\mathbf{a}\cdot(u,v))$, intended to better represent sharp edges

- Use spherical harmonics for direction-dependent background color

- Initialize splat scale and orientation based on SVD of neighbor SfM points


### To-do

- Fast per-pixel depth sorting
  - Has potential to stop "pops" and significantly reduce the number of splats, as per [this paper](https://arxiv.org/abs/2402.00525)
  - Has potential to improve depth regularization: So far L2 depth regularization performs worse than L1 regularization based on center of splats

- Better Gaussian splitting/duplication

- Anti-aliasing

- Per-image embedding to model lighting change and blur

- Try non-planar splats

- Support camera distortion


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