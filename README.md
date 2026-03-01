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
Install [Nerfstudio](https://docs.nerf.studio/) using `pip install nerfstudio`. We recommend a more recent CUDA and PyTorch version combination. Clone the repository and run the commands:

```
cd spirulae-splat/
git submodule update --init
pip install -e . --no-build-isolation -v
ns-install-cli
```

If you are running out of system resources in the `pip install` step, set environment variable `MAX_JOBS` to a lower number (default is max number of concurrent CPU threads).

Use default (master) branch for a stable version. Use `dev` branch if you want to try some more recent features.

## Running the method
This repository creates a new Nerfstudio method named "spirulae". To train with it, run the command:
```
ns-train spirulae --data [DATASET_PATH]
```

By default, spirulae-splat uses all available images for training. To support `ns-eval`, for nerfstudio dataset, use the following command for training:
```
ns-train spirulae --data [DATASET_PATH] nerfstudio-data --eval-mode interval --eval-interval 8
```

## Quick start

Gaussian representation
- Change number of Gaussians: `--pipeline.model.mcmc_cap_max 6000000` (default 1000000)
- Change SH degree: `--pipeline.model.sh_degree 1` (default 3)
- Set primitive using `--pipeline.model.primitive` (default `3dgut`, change to `3dgs` or `mip` for potentially better compatibility and faster training)

Exposure/WB correction
- Bilateral grid is enabled by default, disable using `--pipeline.model.use_bilateral_grid False`
- Change shape from default `(16, 16, 8)` to `(8, 8, 4)` using `--pipeline.model.bilagrid_shape 8 8 4` (sometimes give less color shift)
- `--pipeline.model.use_ppisp True` to enable PPISP, for less color shift but more floaters when there's environment lighting change (Note: this feature is more stable in `dev` branch, which also implements predicting PPISP parameters with bilateral grid)

Fisheye images
- Fisheye images are automatically handled using 3DGUT
- `--pipeline.model.mcmc_max_screen_size 0.15` is enabled by default for compatibility conventional viewers. Increase it for potentially better quality in built-in viewer, decrease it for better compatibility with other viewers (e.g. SuperSplat viewer, especially if you notice spikes or large floaters)
- To fall back to a fisheye-GS style method: set `--pipeline.model.primitive` to `3dgs` (not anti-aliased), or `mip` (anti-aliased).

Background control
- By default, background is trained with a skybox, which effectively removes floaters from sky for outdoor scenes.
- To set to conventional black background, use `--pipeline.model.background_color black --pipeline.model.train_background_color False`
- To discourage transparency in background for indoor scenes, use `--pipeline.model.randomize_background True`
- If mask is provided, set `--pipeline.model.apply_loss_for_mask` to True to mask e.g. sky, background, and False to mask e.g. people and cars.

Training very large-scale scenes
- Cache images on disk for large datasets (instead of loading everything into RAM): `--pipeline.model.cache_images disk` (default: `cpu-pageable`)
- If you notice "splat blobs", increase `--pipeline.model.relative_scale` aggressively (default 1.0)
- Batching is enabled by default for scenes with large number of images, can be configured with `--pipeline.datamanager.max_batch_per_epoch` (default 768)

For consistency with academic benchmarks (if you want an implementation that's as faithful to original 3DGS MCMC as possible)
- `ns-train spirulae --data [DATASET_PATH] --pipeline.model.primitive 3dgs --pipeline.model.background_color black --pipeline.model.train_background_color False --pipeline.model.alpha_reg_weight 0.0 --pipeline.model.use_bilateral_grid False --pipeline.model.use_ppisp False --pipeline.model.ssim_lambda 0.2 --pipeline.model.mcmc_max_screen_size inf --pipeline.model.erank_reg 0.0 --pipeline.model.erank_reg_s3 0.0 --pipeline.model.alpha_loss_weight 0.0 --pipeline.model.alpha_loss_weight_under 0.0 --pipeline.model.stop_refine_at 25000 --pipeline.datamanager.max_batch_per_epoch 999999 nerfstudio-data --eval-mode interval --eval-interval 8`
- This may give you higher PSNR, SSIM, etc. on standard benchmark datasets (e.g. Mip-NeRF 360). Not recommended for large-scale scenes, or scenes captured with large lighting change.
- Note that this is still not 100% identical to official 3DGS MCMC. Notably differences: sorting using ray depth instead of linear depth, no SH degree warmup, regularization for quaternion norm, learning rate for Gaussian means, etc.

Unstable features
- Batch many tiny tiles instead of whole images: `ns-train spirulae-patched ...` instead of `ns-train spirulae`
- Validation (early stop training if loss on validation images start to increase): append `nerfstudio-data --validation_fraction 0.1` to the end of training command
- Second-order optimizer using Jacobian-residual product and Hessian diagnonal (in `dev` branch): `ns-train spirulae^2-pos` or `spirulae^2` instead of `spirulae`

Scripts
- Use `scripts/export_ply_3dgs.py` to export PLY
- To process data, use `scripts/process_data_colmap.py` and `scripts/process_data_metashape.py`, will bypass ns-process-data limitations (e.g. `THIN_PRISM_FISHEYE`)
- Use `scripts/mask.py` to generate masks (Example usage: `python3 scripts/mask.py path/to/dataset --prompt "person; car; fisheye border`)
- Use `scripts/predict_geometry.py` to generate depth and normal maps, as well as optionally sky segmentation maps
