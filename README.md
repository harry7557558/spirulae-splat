# spirulae-splat
My custom 3D Gaussian Splatting method for Nerfstudio. Modified the `splatfacto` method in official Nerfstudio implementation.

### Currently supports
- MCMC densification
- Bilateral grid (fully fused CUDA implementation)
- Masking (sky mode and people/car mode)
- Use SH for background color, with regularization to balance sky removal and discouraging transparency
- Depth and normal supervision using monocular geometry models
- Camera models: perspective and equidistant fisheye (supports >180° fov), fully supports radial, tangential, and thin prism distortion coefficients
- 3DGUT, with option for improving compatibility with vanilla 3DGS viewers
- [PhyGaussian](https://arxiv.org/abs/2311.12198) and [erank](https://arxiv.org/abs/2406.11672) regularization to reduce spiky Gaussians
- Batching for very large scenes

### Partially supports
- Vanilla and Abs-GS densification
- Mesh export
- Batch tiles instead of full images (via `spirulae-patched` method)
- Opaque triangle splatting (via `spirulae-triangle` method)

### TODO
- Better results for HDR and various non-sRGB color spaces
- Multi resolution loss
- Better camera pose optimizer
- Faster training (fused/lazy optimizer, Stop-the-pop tile culling, Morton sort Gaussians, etc.)
- Multi-GPU training

### Additional features
- Display number of Gaussians and training loss/PSNR/SSIM in terminal during training
- Reduced system RAM usage for data loader when cached on CPU (up to 2x)
- Fast backward implementation based on Taming-3DGS

### WebGL viewer features (see `webgl`)
- Fisheye distortion (supports >180deg fisheye)
- Compression
- Collision detection (WIP)

### Scripts (see `scripts`)
- Extract frames from video, auto skip blurry frames
- Auto generate depth and normal maps
- Segmentation, with scripts for point (with GUI) and natural language prompts
- And more, etc.

## Installation
Install Nerfstudio (see [instructions](https://docs.nerf.studio/quickstart/installation.html)). Clone this repository and run the commands:

```
cd spirulae-splat/
git submodule update --init
pip install -e . --no-build-isolation
ns-install-cli
```

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
