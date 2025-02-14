# spirulae-splat
My custom method for Nerfstudio. Based on the `splatfacto` method in official Nerfstudio implementation.

### Major changes / Novelties:

- Use 2D splats instead of 3D, add depth and normal regularization, use ray-splat intersection, as per SIGGRAPH 2024 paper [*2D Gaussian Splatting for Geometrically Accurate Radiance Fields*](https://arxiv.org/abs/2403.17888)
  - Use an interpolation function for median depth to reduce discontinuities and encourage denser gradient
  - For depth regularization, allow blending between two modes: pairwise L2 loss, and L1 loss to rendered depth
  <!-- - L1 loss to rendered mean depth requires two-pass rendering (i.e. one earlier pass for depth) that is 1.3 times as slow to train overall. The code automatically selects one-pass and two-pass rendering based on configuration. -->

- Use polynomial spline kernel splats instead of Gaussians, as inspired by [kernels used in computational physics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics)
  - This increases rendering speed by zeroing out opacities beyond a certain radius
  - A splat with higher overall opacity produces fewer Gaussians, especially when combined with culling splats with [absgrad](https://arxiv.org/abs/2404.10484) below a threshold
  - $\alpha=\max(1-r^2,0)$ empirically produces 1.5x-2x faster training and inference compared to Gaussians, with 0.2-0.5 drop in PSNR compared to Gaussian kernel for a similar number of splats

- Use [cylindrical harmonics](https://en.wikipedia.org/wiki/Cylindrical_harmonics) to model uv-dependent color
  - Base+SH color multiplied by the sigmoid of 2D "polar" harmonics
  - Able to better capture fine color details and potentially reduce export file size on scenes with fine detail
  - Optionally use absolute gradient of CH coefficients as a criteria to split and duplicate splats
  - Slows down training by about a half, due to resource spent on Bessel function evaluation and loading CH coefficients from global memory

<!-- - Support [Markov Chain Monte Carlo (MCMC)](https://arxiv.org/abs/2404.09591) for adaptive control of splats, with relocation designed for polynomial splats -->

- Introduce splat anisotropy by multiplying opacity by a smoothstep of UV directional dot product $1-\mathrm{smoothstep}(\mathbf{a}\cdot(u,v))$, intended to better represent sharp edges

- Handle exposure change by fitting output image to ground-truth image with a linear model before evaluating loss
  - Supported models: `gt ~ k * pred` with scalar or per-channel k, `log(gt) ~ k * log(pred) + b` with scalar or per-channel k and b, `gt ~ A * pred` with 3x3 matrix A, `log(gt) ~ A * log(pred) + b` with matrix A and vector b
  - Closed-form solution exists using linear least squares
  - To ensure predicted images look good during inference, the (weighted) mean and covariance of splat colors is fitted to mean and covariance of training image colors. This is alternatively achieved by adding a small loss component from pre-adjustment image when evaluating photometric loss with ground-truth image.

<!-- - Use an adaptive densification threshold based on absgrad, for consistency in number of splats across scenes and hyperparameter sets -->

- Use spherical harmonics for direction-dependent background color

- Initialize splat scale and orientation based on SVD of neighbor SfM points


### Ideas / To-do

- Fast per-pixel depth sorting
  - Has potential to stop "pops" and significantly reduce the number of splats, as per [this paper](https://arxiv.org/abs/2402.00525)
  - Has potential to improve depth regularization: So far L2 depth regularization on intersected depth performs worse than L1 regularization based on center of splats
  - Alternatively try ray tracing? Might not be fast for training, but seems to be an inference solution that's friendly with existing graphics pipelines

- Better model for view-dependent colors for glossy surfaces, e.g. floor where most camera views have low elevation angle, result doesn't generalize well to bird's eye views
  - Idea: low-degree SH with Fresnel reflectance?

- Depth supervision
   - Fit rendered depth to depth predicted using MVS or a foundational depth model
   - Special attention needed for biased depth, possibly using an approach similar to how exposure is handled

- Speed up uv-dependent color
  - Idea: use an [interpolated circular grid](https://www.desmos.com/calculator/0drgnclvod) that allow both low-latency access and accomodation for circular splat shape
  - Also [apply texture to opacity](https://arxiv.org/abs/2408.16982)?

- Anti-aliasing

- Try non-planar splats, like quadratic patches

- Support camera distortion, either in Gaussian projection or using ray tracing


## Installation
Install Nerfstudio (see [instructions](https://docs.nerf.studio/quickstart/installation.html)). Clone this repository and run the commands:

```
cd spirulae-splat/
git submodule update --init
pip install -e .
ns-install-cli
```

I've tested with nerfstudio 1.1.5, Python 3.8.10 and 3.10.12, Ubuntu 20.04 and 22.04, CUDA 11.8 and 12.5, both in Conda and in global environment. It may also work with other versions and platforms.

## Running the new method
This repository creates a new Nerfstudio method named "spirulae". To train with it, run the command:
```
ns-train spirulae --data [DATASET_PATH]
```

By default, spirulae-splat uses all available images for training. To support `ns-eval`, for nerfstudio dataset, use the following command for training:
```
ns-train spirulae --data [DATASET_PATH] nerfstudio-data --train-split-fraction 0.9
```

<!-- ### Dataset preparation
Some scripts I use to generate datasets from videos can be found [here](https://github.com/harry7557558/Graphics/tree/master/mapping/colmap_nerfstudio). See [here](https://github.com/harry7557558/Graphics/tree/master/mapping/video_imu_alignment) for a tool (under development) that recovers scene scale and orientation using IMU data. -->
