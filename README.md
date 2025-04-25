# spirulae-splat
My custom 2DGS method for Nerfstudio. Modified the `splatfacto` method in official Nerfstudio implementation.

<!-- This project is largely (although not entirely) a hobby project to test out ideas I have. I also use it to generate splats of some real-world scenes around me. -->

### Major changes / Novelties:

- Use 2D splats instead of 3D, add depth and normal regularization, use ray-splat intersection, as per SIGGRAPH 2024 paper [*2D Gaussian Splatting for Geometrically Accurate Radiance Fields*](https://arxiv.org/abs/2403.17888)
  - Use an interpolation function for median depth to reduce discontinuities and encourage denser gradient
  - For depth regularization, allow blending between two modes: pairwise L2 loss, and L1 loss to rendered depth
  <!-- - L1 loss to rendered mean depth requires two-pass rendering (i.e. one earlier pass for depth) that is 1.3 times as slow to train overall. The code automatically selects one-pass and two-pass rendering based on configuration. -->

- Per-pixel sorting: Sort intersections by depth for each intersected splat
  - A rasterization-like pass calculates intersected depths and splat indices, and indices are sorted by depth with utilization of shared memory for speed
  - During training, sorted indices are cached and reused in backward
  - Significantly improves generalizability to unseen views; for example:
    - Trained with camera views at small angles with surface, but rendered from a view perpendicular to surface (e.g. captured with a hand-held device but rendered from bird's eye view)
    - Reverse of above (e.g. customarily captured face or mural)
    - Rendered at focal lengths much lower than in training images (e.g. fisheye render for robot simulator)
  - Eliminates popping artifacts in animation (e.g. while navigating in a viewer)
  - Overall training is up to 20% slower as without per-pixel sorting but can be faster for some scenes; (for reference, [Stop-the-pop](https://arxiv.org/abs/2402.00525v3), a more complex per-pixel sorting method, reported consistent 4% slow down)
  - Makes inference a little tricky &ndash; A scene trained with per-pixel sorting doesn't look nice when rendered without per-pixel sorting, causes compatibility issues for existing viewers

- Use polynomial spline kernel splats instead of Gaussians, as inspired by [kernels used in computational physics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics)
  - This increases rendering speed by zeroing out opacities beyond a certain radius
  - A splat with higher overall opacity produces fewer Gaussians, especially when combined with culling splats with [absgrad](https://arxiv.org/abs/2404.10484) below a threshold
  - $\alpha=\max(1-r^2,0)$ empirically produces 1.5x-2x faster training and inference compared to Gaussians, with 0.2-0.5 drop in PSNR compared to Gaussian kernel for a similar number of splats

- Use [cylindrical harmonics](https://en.wikipedia.org/wiki/Cylindrical_harmonics) to model uv-dependent color
  - Base+SH color multiplied by the sigmoid of 2D "polar" harmonics
  - Able to better capture fine color details and potentially reduce export file size on scenes with fine detail
  - Optionally use absolute gradient of CH coefficients as a criteria to split and duplicate splats
  - Slows down training by about a half, due to resource spent on Bessel function evaluation and loading CH coefficients from global memory

- Set a maximum number of batches per epoch and dynamically set batch size based on number of training images, significantly improves performance for scenes with thousands of scattered images (e.g. large-scale indoor scenes)

- Support [Markov Chain Monte Carlo (MCMC)](https://arxiv.org/abs/2404.09591) for adaptive control of splats, with a relocation formula specifically designed for polynomial splats by matching two moments

<!-- - Introduce splat anisotropy by multiplying opacity by a smoothstep of UV directional dot product $1-\mathrm{smoothstep}(\mathbf{a}\cdot(u,v))$, intended to better represent sharp edges

DEPRECATED: although this seems to improve metrics (e.g. PSNR) for a similar number of splats, visualization shows anisotropy isn't actually high near sharp edges; deprecate to save memory footprint, and considering similar but more configurable methods (like https://arxiv.org/abs/2408.16982) already exist -->

- Handle exposure change by fitting output image to ground-truth image with a linear model before evaluating loss
  - Supported models: `gt ~ k * pred` with scalar or per-channel k, `log(gt) ~ k * log(pred) + b` with scalar or per-channel k and b, `gt ~ A * pred` with 3x3 matrix A, `log(gt) ~ A * log(pred) + b` with matrix A and vector b
  - Closed-form solution exists using linear least squares
  - To ensure predicted images look good during inference, the (weighted) mean and covariance of splat colors is fitted to mean and covariance of training image colors. This is alternatively achieved by adding a small loss component from pre-adjustment image when evaluating photometric loss with ground-truth image.

- Depth supervision
  - Directly supervising depth, normal, and alpha with depth predicted by [MoGe](https://arxiv.org/abs/2410.19115) significantly improves performance on low-texture and reflective surfaces.
  - Before MoGe, an attempt was made to utilize [Depth Anything v2](https://arxiv.org/abs/2406.09414) predictions, a generalized linear model fitting was applied to correct distortion before evaluating loss, although result was poor for most outdoor scenes.

<!-- - Use an adaptive densification threshold based on absgrad, for consistency in number of splats across scenes and hyperparameter sets -->

- Support fisheye camera
  - Exact fisheye render is achieved by using a pre-computed undistortion map in rasterization, with distorted bounding box from projection.
  - This is different from [Fisheye-GS](https://arxiv.org/abs/2409.04751) that distorts splats using first-order approximation during projection, which has issues that our method don't have:
    - Large approximation error near the 180° circle &ndash; a distorted splat should be a crescent, not an ellipse
    - Culls large splats with parts behind camera that could otherwise contribute radiance, which is much more common for 180° fisheye cameras than for conventional cameras

- Use spherical harmonics for direction-dependent background color
  - Intended to better represent sky in outdoor scenes
  - Automatically removes floaters from sky, although sometimes cause splat opacity to vanish
  - Improves metrics (like PSNR) with fewer splats, but reduces generalizability to unseen views

- Initialize splat scale and orientation based on SVD of neighbor SfM points


### Ideas / To-do

- Train on low-quality images, as well as general captures from general users
  - Mobile video captures often have mixed high-quality images and blurry images.
    - List of other common degradations: shot noise, ringing overshoot, JPEG/AVC/HEVC artifacts, glass reflection, bloom, flare
  - Ideas
    - Physically model motion/defocus blur? (Google "gaussian splatting deblur" and check first few results)
    - Process input images using an image restoration model like [this](https://github.com/swz30/Restormer)? (Need to train one for general degradation. I've tried to train one with [albumentations](https://albumentations.ai/) degradation, it worked for shot noise and JPEG artifacts but not real-world motion blur. Possibly need good synthetic data generated in a way similar to [this](https://arxiv.org/abs/1711.07064) to handle motion blur.)
    - Use some loss or regularization term to encourage splats to learn from high-quality parts and ignore low-quality parts, when there's mixed low-quality and high-quality training data
  - May also be useful to handle uncalibrated distortion, as in following cases:
    - Captured through car window or curved museum glass
    - Rolling shutter cameras under motion
    - Cameras with intrinsics that don't fit perfectly to COLMAP models (e.g. cheap Arduino cameras, clipped external lens for smartphones)

- Speed up bilateral grid from splatfacto with a fully-fused CUDA implementation

- Fit a splat in 8 scalars for memory speed
  - Current implementation represents geometry of projected splat with position and two axes, which appears to be faster than 2DGS implemented in gsplat, but it still takes high memory footprint.
  - Quick benchmark shows loading two `float4` from global memory is 2 times faster than loading three non-contiguous `float3`.
  - Idea: 5 scalars for projected ellipse, plus $z=a(x-x_0)+b(y-y_0)+c$ for screen-space depth, which doesn't introduce additional singularity; Pack them into one contiguous tensor.

- Better model for view-dependent colors for glossy surfaces, e.g. floor where most camera views have low elevation angle, result doesn't generalize well to bird's eye views
  - Idea: low-degree SH with Fresnel reflectance?

- Speed up uv-dependent color
  - Idea: use an [interpolated circular grid](https://www.desmos.com/calculator/0drgnclvod) that allow both low-latency access and accomodation for circular splat shape
  - Also [apply texture to opacity](https://arxiv.org/abs/2408.16982)?

- Anti-aliasing, like in [Mip-Splatting](https://arxiv.org/abs/2311.16493)

- Try non-planar splats, like quadratic patches

- *\*\*Add some pics to this README*

## Installation
Install Nerfstudio (see [instructions](https://docs.nerf.studio/quickstart/installation.html)). Clone this repository and run the commands:

```
cd spirulae-splat/
git checkout stable
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

## Dataset preparation

See `scripts` folder for some scripts I use to prepare datasets.
