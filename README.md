# spirulae-splat
My custom 3D Gaussian Splatting method for Nerfstudio. Modified the `splatfacto` method in official Nerfstudio implementation.

### Currently supports
- MCMC densification, with option to enhance background details
- Bilateral grid and PPISP for exposure/WB correction
- Masking (sky mode and people/car mode)
- Use SH for background color, with regularization to balance sky removal and discouraging transparency
- Depth and normal supervision using monocular geometry models
- Camera models: perspective and equidistant fisheye (supports >180° fov), fully supports radial, tangential, and thin prism distortion coefficients
- Primitives: 3DGUT (default), 3DGS (original and anti-aliased); 3DGUT provides options for compatibility with viewers that assume vanilla 3DGS
- Caching training images on disk for extremely large scenes
- Select a subset of images for validation, early stop training when validation loss starts to increase
- [PhyGaussian](https://arxiv.org/abs/2311.12198) and [erank](https://arxiv.org/abs/2406.11672) regularization to reduce spiky Gaussians
- Batching for very large scenes
- Display number of Gaussians and training loss/PSNR/SSIM in terminal during training

### Partially supports
- Training on images in linear and various wide-gamut color spaces, with EXR or 16-bit PNG images as input
- Second-order optimizer (loosely based on [this paper](https://arxiv.org/abs/2602.00395))
- 2DGS-like depth regularization to discourage floaters
- Batch tiles instead of full images (via `spirulae-patched` method)
- Opaque triangle splatting (via `spirulae-triangle` method)
- Vanilla and Abs-GS densification
- Mesh export

### To-Do
- Better camera pose/intrinsics optimizer
- Faster training (Stop-the-pop tile culling, Morton sort Gaussians, etc.)
- Multi-GPU training

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

The `pip install` step may take a few minutes. If you are running out of system resources during installation, set environment variable `MAX_JOBS` to a lower number (default is max number of concurrent CPU threads).

Use default (master) branch for a stable version. Use `dev` branch if you want to try some more recent features.

## Running the method
This repository creates a new Nerfstudio method named "spirulae". To train with it, run the command (replace `[DATASET_PATH]` and `<additional args>` with actual dataset path and additional arguments):
```
ns-train spirulae --data [DATASET_PATH] <additional args>
```

By default, spirulae-splat uses all available images for training. To support `ns-eval`, for nerfstudio dataset, use the following command for training:
```
ns-train spirulae --data [DATASET_PATH] <additional args> nerfstudio-data --eval-mode interval --eval-interval 8
```

## Quick start

Presets
- Spirulae-splat provides presets. Run `ns-train <preset name> --data [DATASET_PATH] <additional args>` to use a preset.
- List of presets:
    - `spirulae-preset-confined`: For confined environments (indoor scenes or outdoor scenes without much sky visible)
    - `spirulae-preset-open`: For open environments (outdoor scenes with sky visible)
    - `spirulae-preset-confined-low-texture`: For confined environments with large textureless surfaces
    - `spirulae-preset-open-low-texture`: For open environments with large textureless surfaces
    - `spirulae-preset-centered-object`: For small centered objects (e.g. small item on turntable, human avatar)
    - `spirulae-preset-academic-baseline`: Use this to replicate academic baseline (3DGS MCMC) as faithful as it can be
- Notes:
    - `spirulae-preset-academic-baseline` is generally the fastest and most memory efficient option and the one that achieves the highest PSNR and SSIM on standard benchmarks (e.g. Mip-NeRF 360). For better visual quality on real-world datasets and compatibility across different viewers, choose a different preset.
        - Note that this does not 100% replicate academic baseline. Known mismatches: no SH degree warm up, quaternion optimizer and initialization. When benchmarked on Mip-NeRF 360, it generally achieves better SSIM and LPIPS for outdoor scenes but worse PSNR for indoor scenes.
    - `low-texture` presets are for large surfaces with nearly no texture (e.g. full white wall). For scenes with moderate texture, you can likely get better visual results without `low-texture`.
    - `open` presets will train a sky box, and the PLY export script will export it to an equirectangular map. Choose a `confined` preset if you wish to keep sky as splats.

Gaussian representation
- Change number of Gaussians: `--pipeline.model.mcmc_cap_max 6000000` (default 1000000)
- Change SH degree: `--pipeline.model.sh_degree 1` (default 3)
- Set primitive using `--pipeline.model.primitive` (default `3dgut`, change to `3dgs` or `mip` for potentially better compatibility across viewers and faster training)

Exposure/WB correction
- Bilateral grid is enabled by default, disable using `--pipeline.model.use_bilateral_grid False`
- Change shape from default `(16, 16, 8)` to `(8, 8, 4)` using `--pipeline.model.bilagrid_shape 8 8 4` (sometimes gives less color shift)
- `--pipeline.model.use_ppisp True` to enable PPISP, for less color shift but more floaters when there's environment lighting change
- Enable bilateral grid and set `--pipeline.model.bilagrid_type ppisp` to make bilateral grid predict PPISP parameters (exposure and color). Generally achieves less color shift for low-texture surfaces.

Distorted/Fisheye images
- Images with fisheye camera model (as well as general distorted images) are automatically handled using 3DGUT
- `--pipeline.model.mcmc_max_screen_size 0.15` is enabled by default for compatibility conventional viewers. Increase it for potentially better quality in built-in viewer, decrease it for better compatibility with other viewers (e.g. SuperSplat viewer, especially if you notice spikes or large floaters)
- To fall back to a Fisheye-GS style method: set `--pipeline.model.primitive` to `3dgs` (not anti-aliased), or `mip` (anti-aliased).
- Supported camera models: perspective, equidistant fisheye (supports >180° fov); Supported distortion parameters: k1-k4, p1, p2, sx1, sy1, b1, b2. For better reliability, use `scripts/process_data_(colmap|metashape).py` (instead of `ns-process-data`) to process data.

Background control
- By default, background is trained with a skybox, which removes floaters from sky for outdoor scenes, but can introduce unwanted transparency for some datasets.
- To set to conventional black background, use `--pipeline.model.background_color black --pipeline.model.train_background_color False`
- To discourage transparency in background for indoor scenes, use `--pipeline.model.randomize_background True`
- If mask is provided, set `--pipeline.model.apply_loss_for_mask` to True to mask e.g. sky, background, and False to mask e.g. people and cars.

Training very large-scale scenes
- Cache images on disk for large datasets (instead of loading everything into RAM that cause OOM): `--pipeline.model.cache_images disk` (default: `cpu-pageable`)
- If you notice "splat blobs" with a `low-texture` preset, increase `--pipeline.model.relative_scale` aggressively (default 10.0 for open and 1.0 otherwise)
- Experimental multi-resolution loss that helps with convergence with high-resolution images: `--pipeline.model.num_loss_scales 2` (default 0)
- Batching for scenes with large number of images can be configured with `--pipeline.datamanager.max_batch_per_epoch` (default 800), which automatically enables batching when number of input images is above this number.
- Tips for training with limited VRAM
    - In batching mode, all images are processed in a single batch by default. Set `--pipeline.datamanager.split_batch True` to process one image at once, which can save VRAM a lot if you have a large number of high-resolution images.
    - Use `--pipeline.model.optimizer_offload all` to offload Adam optimizer momentum to CPU, can save VRAM significantly at cost of slower training.
    - On some platforms, setting environment variable `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True` can reduce VRAM usage by 10%-20% with minimum speed overhead. Setting `PYTORCH_NO_CUDA_MEMORY_CACHING=1` can reduce peak VRAM even more but with slight speed overhead.
    - If you are not using depth and normal supervision, setting `--pipeline.model.use_bilateral_grid_for_geometry False` may save VRAM slightly.
- To skip viewer thumbnail loading (if it takes too long in the beginning of training), append `nerfstudio-data --load_thumbnails False` to the end of training command.

Unstable features
- Training on images in linear color spaces: `--pipeline.model.use_linear_color_space True`; Wide-gamut color spaces: `--pipeline.model.image_color_gamut ACEScg` (supports `ACES2065-1`, `ACEScg`, `Rec.2020`, `AdobeRGB`); Specify `--pipeline.model.convert_initial_point_cloud_color True` if colors in initial point cloud is in sRGB.
- Batch many tiny tiles instead of whole images: `ns-train spirulae-patched ...` instead of `ns-train spirulae`
- Validation (early stop training if loss on validation images start to increase): append `nerfstudio-data --validation_fraction 0.1` to the end of training command
- Second-order optimizer using Jacobian-residual product and Hessian diagonal: `ns-train spirulae^2-pos` (more stable) or `spirulae^2` (less stable) instead of `spirulae`. We also provide presets `spirulae^2-preset-confined` and `spirulae^2-preset-open` for the corresponding presets with `spirulae^2` methods, which otherwise run on `spirulae^2-pos`.
- 2DGS-like depth regularization to discourage floaters: `--pipeline.model.depth_reg_weight 0.01`. Similar regularization can also be applied to RGB by setting `--pipeline.model.rgb_distortion_reg_weight` to a positive value.

Scripts
- Use `scripts/export_ply_3dgs.py` to export PLY
- To process data, use `scripts/process_data_colmap.py` and `scripts/process_data_metashape.py`, will bypass `ns-process-data` limitations (e.g. `THIN_PRISM_FISHEYE`)
- Use `scripts/mask.py` to generate masks (Example usage: `python3 scripts/mask.py path/to/dataset --prompt "person; car; fisheye border"`). By default, this runs on [lang-sam](https://github.com/luca-medeiros/lang-segment-anything) model. Use `--model sam3` to switch to [SAM 3](https://github.com/facebookresearch/sam3) model for often better results (may require applying for access and logging in to Huggingface).
- Use `scripts/predict_geometry.py` to generate depth and normal maps using [Metric3D v2](https://github.com/YvanYin/Metric3D) model, and optionally sky segmentation maps with various model options.
- Use `scripts/extract_frames.py` to extract frames from a video, while skipping blurry frames. Supports various video formats, including most `.mp4`, `.mov`, and `.insv` videos.
- `scripts/downscale_dataset.py`, `scripts/undistort_dataset.py`: self-explanatory

<!-- ## Trivia

spirulae-splat is named after now-inactive project [spirulae](https://github.com/harry7557558/spirulae), which is named after the [deep-ocean cephalopod mollusk](https://en.wikipedia.org/wiki/Spirula).
 -->
