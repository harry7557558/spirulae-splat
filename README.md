# spirulae-splat

> **Note to existing users:** Since mid-April 2026, Nerfstudio and GSplat dependencies are no longer needed. The interface has gone through major change (see "Quick Start" section below). This branch (`dev`) will act as a staging branch for next few weeks to expose any issue, before being merged into the `master` branch. If you see error or notice significant quality degrade compared to the `master` branch, please send me a message on Discord (@spirulae) or email (the one in most of my git commits).

My custom method used to train 3D Gaussian Splatting (3DGS) models.

If you find spirulae-splat helpful for your research, please cite corresponding works (see "Acknowledgement" section below).

If you share 3DGS models trained with spirulae-splat, or incorporate any feature, idea, or the software itself into your code, product, or service, a mention of spirulae-splat with a link to this page is much appreciated.

I'm also considering adding a few visuals to enhance this README. If you have cool splats made with spirulae-splat and are willing to share a few images publicly, feel free to reach out.

### Currently supports
- Unified densification strategy combining elements of MCMC and [IGS+](https://arxiv.org/abs/2603.08661)/[MRNF](https://lichtfeld.io/)
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
- Batch tiles instead of full images (via `3dgs-patched` method)
- Opaque triangle splatting (via `triangle` method)
- Mesh export

### Scripts (see `scripts`)
- Extract frames from video, auto skip blurry frames
- Auto generate depth and normal maps
- Segmentation, with scripts for point (with GUI) and natural language prompts
- Downscale/Undistort datasets
- And more, etc.

## Installation
Make sure you have a recent version of CUDA and PyTorch installed. Clone the repository and run the commands:

```
cd spirulae-splat/
git submodule update --init
pip install -e . --no-build-isolation -v
```

The `pip install` step may take a few minutes. If you are running out of system resources during installation, set environment variable `MAX_JOBS` to a lower number (default is max number of concurrent CPU threads).

Use default (master) branch for a stable version. Use `dev` branch if you want to try some more recent features.

If you see a help message when running `spirulae-train --help`, you likely have installed it successfully.

## Quick start

Presets
- Spirulae-splat provides presets. Run `spirulae-train <preset name> --data [DATASET_PATH] <additional args>` to use a preset.
- List of presets:
    - `3dgs-confined`: For confined environments (indoor scenes or outdoor scenes without much sky visible)
    - `3dgs-open`: For open environments (outdoor scenes with sky visible)
    - `3dgs-confined-low-texture`: For confined environments with large textureless surfaces
    - `3dgs-open-low-texture`: For open environments with large textureless surfaces
    - `3dgs-centered-object`: For small centered objects (e.g. small item on turntable, human avatar)
    - `academic-baseline`: Use this to replicate academic baseline (3DGS MCMC) as faithful as it can be
- Notes:
    - `academic-baseline` is generally the fastest and most memory efficient option and the one that achieves the highest PSNR and SSIM on standard benchmarks (e.g. Mip-NeRF 360). For better visual quality on real-world datasets and compatibility across different viewers, choose a different preset.
        - Note that this does not 100% replicate academic baseline. Known mismatches: no SH degree warm up, quaternion optimizer and initialization. When benchmarked on Mip-NeRF 360, it generally achieves better SSIM and LPIPS for outdoor scenes but worse PSNR for indoor scenes.
    - `low-texture` presets are for large surfaces with nearly no texture (e.g. full white wall). For scenes with moderate texture, you can likely get better visual results without `low-texture`.
    - `open` presets will train a sky box, and the PLY export script will export it to an equirectangular map. Choose a `confined` preset if you wish to keep sky as splats.

Viewer
- Similar to Nerfstudio and GSplat, you can open the link `http://localhost:7007/` in a web browser to view training progress.

Gaussian representation
- Change number of Gaussians: `--model.mcmc_cap_max 6000000` (default 1000000)
- Change SH degree: `--model.sh_degree 1` (default 3)
- Set primitive using `--model.primitive` (default `3dgut`, change to `3dgs` or `mip` for potentially better compatibility across viewers and faster training)

Exposure/WB correction
- Bilateral grid is enabled by default, disable using `--model.no_use_bilateral_grid`
- Change shape from default `(16, 16, 8)` to `(8, 8, 4)` using `--model.bilagrid_shape 8 8 4` (sometimes gives less color shift)
- `--model.use_ppisp` to enable PPISP, for less color shift but more floaters when there's environment lighting change
- Enable bilateral grid and set `--model.bilagrid_type ppisp` to make bilateral grid predict PPISP parameters (exposure and color). Generally achieves less color shift for low-texture surfaces.

Distorted/Fisheye images
- Images with fisheye camera model (as well as general distorted images) are automatically handled using 3DGUT
- `--model.mcmc_max_screen_size 0.15` is enabled by default for compatibility conventional viewers. Increase it for potentially better quality in built-in viewer, decrease it for better compatibility with other viewers (e.g. SuperSplat viewer, especially if you notice spikes or large floaters)
- To fall back to a Fisheye-GS style method: set `--model.primitive` to `3dgs` (not anti-aliased), or `mip` (anti-aliased).
- Supported camera models: perspective, equidistant fisheye (supports >180° fov); Supported distortion parameters: k1-k4, p1, p2, sx1, sy1, b1, b2. For better reliability, use `scripts/process_data_(colmap|metashape).py` (instead of `ns-process-data`) to process data.

Background control
- By default, background varies depends on preset, which can be black, random noise, or a sky box. To configure custom background:
- To disable random noise (which intends to discourage transparency), use `--model.no_randomize_background`; (and otherwise `--model.randomize_background` to enable it)
- To set to conventional black background, use `--model.background_color black --model.no_train_background_color`
- To train a skybox, use `--model.background_color gray --model.background_sh_degree 4`
- If mask is provided, set `--model.apply_loss_for_mask` to True to mask e.g. sky, background, and False to mask e.g. people and cars.

Training very large-scale scenes
- Cache images on disk for large datasets (instead of loading everything into RAM that cause OOM): `--datamanager.cache_images disk` (default: `cpu-pageable`)
- If you notice "splat blobs" with a `low-texture` preset, increase `--model.relative_scale` aggressively (default 10.0 for open and 1.0 otherwise)
- Experimental multi-resolution loss that helps with convergence with high-resolution images: `--model.num_loss_scales 2` (default 0)
- Batching for scenes with large number of images can be configured with `--datamanager.max_batch_per_epoch` (default 800), which automatically enables batching when number of input images is above this number.
- Tips for training with limited VRAM
    - In batching mode, all images are processed in a single batch by default. Set `--datamanager.split_batch True` to process one image at once, which can save VRAM a lot if you have a large number of high-resolution images.
    - Use `--model.optimizer_offload all` to offload Adam optimizer momentum to CPU, can save VRAM significantly at cost of slower training.
    - On some platforms, setting environment variable `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True` can reduce VRAM usage by 10%-20% with minimum speed overhead. Setting `PYTORCH_NO_CUDA_MEMORY_CACHING=1` can reduce peak VRAM even more but with slight speed overhead.
    - If you are not using depth and normal supervision, setting `--model.use_bilateral_grid_for_geometry False` may save VRAM slightly.
- To skip viewer thumbnail loading (if it takes too long in the beginning of training), append `nerfstudio-data --load_thumbnails False` to the end of training command.

Unstable features
- Training on images in linear color spaces: `--model.image_color_is_linear True`; Wide-gamut color spaces: `--model.image_color_gamut ACEScg` (supports `ACES2065-1`, `ACEScg`, `Rec.2020`, `AdobeRGB`, `DCI-3`)
    - If you want splat and image to be in different color spaces, specify `--model.splat_color_is_linear` and `--model.splat_color_gamut` (supports all options in `--model.image_color_gamut` plus `Rec.709`).
    - Specify `--model.convert_initial_point_cloud_color True` if colors in initial point cloud is in sRGB, and color in initial point cloud will be converted to splat's color space. If you don't specify True or False, it will auto decide based on arguments.
- Batch many tiny tiles instead of whole images: `ns-train spirulae-patched ...` instead of `ns-train spirulae`
- Validation (early stop training if loss on validation images start to increase): append `nerfstudio-data --validation_fraction 0.1` to the end of training command
- Second-order optimizer using Jacobian-residual product and Hessian diagonal: `ns-train spirulae^2-pos` (more stable) or `spirulae^2` (less stable) instead of `spirulae`. We also provide presets `spirulae^2-preset-confined` and `spirulae^2-preset-open` for the corresponding presets with `spirulae^2` methods, which otherwise run on `spirulae^2-pos`.
- 2DGS-like depth regularization to discourage floaters: `--model.depth_reg_weight 0.01`. Similar regularization can also be applied to RGB by setting `--model.rgb_distortion_reg_weight` to a positive value.

Scripts
- Use `scripts/export_ply_3dgs.py` to export PLY
- To process data, use `scripts/process_data_colmap.py` and `scripts/process_data_metashape.py`, will bypass `ns-process-data` limitations (e.g. `THIN_PRISM_FISHEYE`)
- Use `scripts/mask.py` to generate masks (Example usage: `python3 scripts/mask.py path/to/dataset --prompt "person; car; fisheye border"`). By default, this runs on [lang-sam](https://github.com/luca-medeiros/lang-segment-anything) model. Use `--model sam3` to switch to [SAM 3](https://github.com/facebookresearch/sam3) model for often better results (may require applying for access and logging in to Huggingface).
- Use `scripts/predict_geometry.py` to generate depth and normal maps using [Metric3D v2](https://github.com/YvanYin/Metric3D) model, and optionally sky segmentation maps with various model options.
- Use `scripts/extract_frames.py` to extract frames from a video, while skipping blurry frames. Supports various video formats, including most `.mp4`, `.mov`, and `.insv` videos.
- `scripts/downscale_dataset.py`, `scripts/undistort_dataset.py`: self-explanatory

## Acknowledgement

Spirulae-splat begins as a fork of:
- Nerfstudio: https://github.com/nerfstudio-project/nerfstudio/
- GSplat: https://github.com/nerfstudio-project/gsplat

In addition, spirulae-splat has been inspired by ideas from the following works:

Representation:
- *Fisheye-GS: Lightweight and Extensible Gaussian Splatting Module for Fisheye Cameras*, by Liao et al. - https://arxiv.org/abs/2409.04751
- *3DGUT: Enabling Distorted Cameras and Secondary Rays in Gaussian Splatting*, by Wu et al. - https://arxiv.org/abs/2412.12507
- *Efficient Perspective-Correct 3D Gaussian Splatting Using Hybrid Transparency*, by Hahlbohm et al. - https://fhahlbohm.github.io/htgs/
- *Triangle Splatting+: Differentiable Rendering with Opaque Triangles*, by Held et al. - https://arxiv.org/abs/2509.25122
- *Sparse Voxels Rasterization: Real-time High-fidelity Radiance Field Rendering*, by Sun et al. - https://arxiv.org/abs/2412.04459

Densification:
- *3D Gaussian Splatting as Markov Chain Monte Carlo*, by Kheradmand et al. - https://arxiv.org/abs/2404.09591
- *Taming 3DGS: High-Quality Radiance Fields with Limited Resources*, by Mallick et al. - https://arxiv.org/abs/2406.15643
- *Improving Densification in 3D Gaussian Splatting for High-Fidelity Rendering*, by Deng et al. - https://arxiv.org/abs/2508.12313
- *ImprovedGS+: A High-Performance C++/CUDA Re-Implementation Strategy for 3D Gaussian Splatting*, by Jordi Muñoz Vicente - https://arxiv.org/abs/2603.08661
- LichtFeld Studio, by MrNeRF and other contributors - https://lichtfeld.io/
- *AbsGS: Recovering Fine Details for 3D Gaussian Splatting*, by Ye et al. - https://arxiv.org/abs/2404.10484

Other quality improvement:
- *3DGS^2-TR: Scalable Second-Order Trust-Region Method for 3D Gaussian Splatting*, by Hsiao et al. - https://arxiv.org/abs/2602.00395
- *Effective Rank Analysis and Regularization for Enhanced 3D Gaussian Splatting*, by Hyung et al. - https://arxiv.org/abs/2406.11672
- *PhysGaussian: Physics-Integrated 3D Gaussians for Generative Dynamics*, by Xie et al. - https://arxiv.org/abs/2311.12198
- *Tile-wise vs. Image-wise: Random-Tile Loss and Training Paradigm for Gaussian Splatting*, by Zhang et al. - [openaccess.thecvf.com](https://openaccess.thecvf.com/content/ICCV2025/html/Zhang_Tile-wise_vs._Image-wise_Random-Tile_Loss_and_Training_Paradigm_for_Gaussian_ICCV_2025_paper.html)

Optimization:
- *Taming 3DGS: High-Quality Radiance Fields with Limited Resources*, by Mallick et al. - https://arxiv.org/abs/2406.15643
- *Faster-GS: Analyzing and Improving Gaussian Splatting Optimization*, by Hahlbohm et al. - https://arxiv.org/abs/2602.09999 (originally LichtFeld Studio bounty 001)
- *CLM: Removing the GPU Memory Barrier for 3D Gaussian Splatting*, by Zhao et al. - https://arxiv.org/abs/2511.04951
- *Fast BVH Construction on GPUs*, by Lauterbach et al. - https://luebke.us/publications/eg09.pdf

Foundation models:
- *Metric3Dv2: A Versatile Monocular Geometric Foundation Model for Zero-shot Metric Depth and Surface Normal Estimation*, by Hu et al. - https://arxiv.org/abs/2404.15506
- Language Segment-Anything - https://github.com/luca-medeiros/lang-segment-anything
- *SAM 3: Segment Anything with Concepts*, by Meta Research - https://github.com/facebookresearch/sam3
- *Depth Anything 3: Recovering the Visual Space from Any Views*, by Lin et al. - https://arxiv.org/abs/2511.10647

Additionally, thanks Kerbl et al. for introducing [3DGS](https://repo-sam.inria.fr/fungraph/3d-gaussian-splatting/) in the first place, as well as various members from MrNeRF & Brush (and previously Nerfstudio) Discord communities for providing feedback and support.

<!-- ## Trivia

Spirulae-splat is named after the now-inactive project [spirulae](https://github.com/harry7557558/spirulae), which is named after the [deep-ocean cephalopod mollusk](https://en.wikipedia.org/wiki/Spirula).
 -->
