# Dataset Preparation Scripts

Scripts I use to prepare datasets for spirulae-splat as well as other nerfstudio methods.

Mostly tested on Ubuntu 20.04 and 22.04.


## Image extraction from video

Extract frames from a video, pick ones with minimum motion blur. Run `python3 extract_frames.py --help` for details.

After extracting frames, you can use `run_colmap.bash` to process datasets into nerfstudio format. You may edit the bash script to customize it for your needs.

### Generating segmentation masks using SAM2

I experimented with scanning small items where all-direction views can be useful, which unavoidably involves rotating the object while capturing and introduce distractors (like hands and moving background). See [here](https://www.youtube.com/watch?v=ugqZpQzKix8) for an example video capture. This feature may also be helpful in the case of turntable capture, or when one needs to remove background.

To build reconstruction from this type of capture:
- Use `--sam 1` while extracting frames from video with the above script. It should generate `run_colmap.bash` and `run_sam2.bash` in the dataset directory.
- Install [Segment Anything 2](https://github.com/facebookresearch/sam2) following instruction and download checkpoints.
  - In my case, I modified `setup.py` to disable auto PyTorch update, and I ran it in the same environment as nerfstudio without issue.
- From the `sam2` directory, run `path/to/dataset/run_sam2.bash`. It should launch a Gradio site that lets you generate and save segmentation masks.
  - See [this issue](https://github.com/harry7557558/spirulae-splat/issues/1) for steps about using the Gradio site.
- From the dataset directory, run `run_colmap.bash` to process dataset into nerfstudio format.
  - If reconstruction fails, you may retry with GLOMAP instead of COLMAP.


## Coloring GLOMAP reconstructions

GLOMAP outperforms COLMAP in many datasets, but by default, the generated point clouds don't have color, which can be a minor issue for some Gaussian splatting methods. Use `glomap_color.py` to color the point cloud after running reconstruction.


