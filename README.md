# spirulae-splat
My custom 3D Gaussian Splatting method for Nerfstudio. Modified the `splatfacto` method in official Nerfstudio implementation.

### Currently supports
- Vanilla, Abs-GS, and MCMC densification
- Bilateral grid (fully fused CUDA implementation), as well as another exposure correction method based on linear least squares
- [PhyGaussian](https://arxiv.org/abs/2311.12198) and [erank](https://arxiv.org/abs/2406.11672) regularization to reduce spiky Gaussians
- 3DGUT for fisheye, with an MCMC-friendly way to eliminate large/spiky Gaussians for compatibility with vanilla 3DGS viewers
- Use SH for background color, removes floaters from sky when combined with opacity regularization
- Batching for very large scenes

### Partially supports
- Masking (sky mode and people/car mode)
- Depth and normal supervision using monocular geometry models
- Regularization to balance sky removal and discouraging transparency

### TODO
- Multi resolution loss
- Better camera pose optimizer
- Better mesh export
- Faster training (Taming-GS backward, fused/lazy optimizer, Stop-the-pop tile culling, etc.)
- Multi-GPU training

### Additional features
- Display number of Gaussians and training loss/PSNR/SSIM in terminal during training
- Reduced system RAM usage for data loader when cached on CPU (up to 2x)

### WebGL viewer features (see `webgl`)
- Fisheye distortion (supports >180deg fisheye)
- Compression
- Collision detection (WIP)

### Scripts (see `scripts`)
- Extract frames from video, auto skip blurry frames
- Segmentation using SAM-2
- Process a folder of camera raw images to specific color space

## Installation
Install Nerfstudio (see [instructions](https://docs.nerf.studio/quickstart/installation.html)). Clone this repository and run the commands:

```
cd spirulae-splat/
git checkout stable
git submodule update --init
MAX_JOBS=8 pip install -e . --no-build-isolation
ns-install-cli
```

## Running the method
This repository creates a new Nerfstudio method named "spirulae". To train with it, run the command:
```
ns-train spirulae --data [DATASET_PATH]
```

By default, spirulae-splat uses all available images for training. To support `ns-eval`, for nerfstudio dataset, use the following command for training:
```
ns-train spirulae --data [DATASET_PATH] nerfstudio-data --train-split-fraction 0.9
```
