# spirulae-splat
My custom method for Nerfstudio. Based on the `splatfacto` method in official Nerfstudio implementation.

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
  - Makes inference a little tricky - A scene trained with per-pixel sorting doesn't look nice when rendered without per-pixel sorting, causes compatibility issues for existing viewers

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

<!-- - Support [Markov Chain Monte Carlo (MCMC)](https://arxiv.org/abs/2404.09591) for adaptive control of splats, with relocation designed for polynomial splats -->

<!-- - Introduce splat anisotropy by multiplying opacity by a smoothstep of UV directional dot product $1-\mathrm{smoothstep}(\mathbf{a}\cdot(u,v))$, intended to better represent sharp edges

DEPRECATED: although this seems to improve metrics (e.g. PSNR) for a similar number of splats, visualization shows anisotropy isn't actually high near sharp edges; deprecate to save memory footprint, and considering similar but more configurable methods (like https://arxiv.org/abs/2408.16982) already exist -->

- Handle exposure change by fitting output image to ground-truth image with a linear model before evaluating loss
  - Supported models: `gt ~ k * pred` with scalar or per-channel k, `log(gt) ~ k * log(pred) + b` with scalar or per-channel k and b, `gt ~ A * pred` with 3x3 matrix A, `log(gt) ~ A * log(pred) + b` with matrix A and vector b
  - Closed-form solution exists using linear least squares
  - To ensure predicted images look good during inference, the (weighted) mean and covariance of splat colors is fitted to mean and covariance of training image colors. This is alternatively achieved by adding a small loss component from pre-adjustment image when evaluating photometric loss with ground-truth image.

<!-- - Use an adaptive densification threshold based on absgrad, for consistency in number of splats across scenes and hyperparameter sets -->

- Use spherical harmonics for direction-dependent background color
  - Intended to better represent sky in outdoor scenes
  - Automatically removes floaters from sky, although sometimes cause splat opacity to vanish
  - Improves metrics (like PSNR) with fewer splats, but reduces generalizability to unseen views

- Initialize splat scale and orientation based on SVD of neighbor SfM points


### Ideas / To-do

- Support camera distortion in rendering and training
  - Render large image + distort is inefficient, undistorting fisheye image for training crops boundary and loses information
  - Idea: distort using Jacobian during projection, or calculate distorted bounding box in projection and undistort pixel in rasterization
    - Expect the latter to work better for large distortion; ([visualization](https://www.desmos.com/calculator/eqo66jd3qa))

- Better model for view-dependent colors for glossy surfaces, e.g. floor where most camera views have low elevation angle, result doesn't generalize well to bird's eye views
  - Idea: low-degree SH with Fresnel reflectance?

- Depth supervision
   - Fit rendered depth to depth predicted using MVS or a foundational depth model
   - Special attention needed for biased depth, possibly using an approach similar to how exposure is handled, or bilateral grid as in splatfacto
   - \*\*UPDATE: a prototype was made, while I haven't got to make depth better, supervising with alpha predicted by the depth model does seem to effectively remove floaters from sky (with spherical harmonics background)

- Speed up uv-dependent color
  - Idea: use an [interpolated circular grid](https://www.desmos.com/calculator/0drgnclvod) that allow both low-latency access and accomodation for circular splat shape
  - Also [apply texture to opacity](https://arxiv.org/abs/2408.16982)?

- Anti-aliasing, like in [Mip-Splatting](https://arxiv.org/abs/2311.16493)

- Fit a splat in 8 scalars for memory speed

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
