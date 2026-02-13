# spirulae-splat
My custom 3D Gaussian Splatting method for Nerfstudio. Modified the `splatfacto` method in official Nerfstudio implementation.

### Currently supports
- MCMC densification
- Bilateral grid for exposure/WB correction
- Masking (sky mode and people/car mode)
- Use SH for background color, with regularization to balance sky removal and discouraging transparency
- Depth and normal supervision using monocular geometry models
- Camera models: perspective and equidistant fisheye (supports >180° fov), fully supports radial, tangential, and thin prism distortion coefficients
- Primitives: 3DGUT (default), 3DGS (original and anti-aliased); 3DGUT provides options for compatibility with viewers that assume vanilla 3DGS
- Caching training images on disk for extremely large scenes
- Select a subset of images for validation, early stop training when validation loss starts to increase
- [PhyGaussian](https://arxiv.org/abs/2311.12198) and [erank](https://arxiv.org/abs/2406.11672) regularization to reduce spiky Gaussians
- Batching for very large scenes

### Partially supports
- PPISP for exposure/WB correction
- Training on images in linear and various wide-gamut color spaces
- 2DGS-like depth regularization to discourage floaters
- Batch tiles instead of full images (via `spirulae-patched` method)
- Opaque triangle splatting (via `spirulae-triangle` method)
- Vanilla and Abs-GS densification
- Mesh export

### To-Do
- Better camera pose/intrinsics optimizer
- Scale agnostic optimizer
- Faster training (Stop-the-pop tile culling, Morton sort Gaussians, etc.)
- Multi-GPU training

### Additional features
- Display number of Gaussians and training loss/PSNR/SSIM in terminal during training
- Reduced system RAM usage for data loader when cached on CPU (up to 2x)
- Fast backward implementation based on Taming-3DGS

<!--### WebGL viewer features (see `webgl`)
- Fisheye distortion (supports >180deg fisheye)
- Compression
- Collision detection (WIP)-->

### Scripts (see `scripts`)
- Extract frames from video, auto skip blurry frames
- Auto generate depth and normal maps
- Segmentation, with scripts for point (with GUI) and natural language prompts
- Downscale/Undistort datasets
- And more, etc.

## Installation
Install [Nerfstudio](https://docs.nerf.studio/) using `pip install nerfstudio`. We recommend a more recent CUDA and PyTorch version combination. Clone this repository and run the commands:

```
cd spirulae-splat/
git submodule update --init
pip install -e . --no-build-isolation -v
ns-install-cli
```

If you are running out of system resources in the `pip install` step, set environment variable `MAX_JOBS` to a lower number (default is max number of concurrent CPU threads).

## Running the method
This repository creates a new Nerfstudio method named "spirulae". To train with it, run the command:
```
ns-train spirulae --data [DATASET_PATH]
```

Other experimental methods include `spirulae-patched`, `spirulae-triangle`, etc.

By default, spirulae-splat uses all available images for training. To support `ns-eval`, for nerfstudio dataset, use the following command for training:
```
ns-train spirulae --data [DATASET_PATH] nerfstudio-data --eval-mode interval --eval-interval 8
```
