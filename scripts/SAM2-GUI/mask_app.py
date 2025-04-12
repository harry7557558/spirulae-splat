import torch
# use bfloat16 for the entire notebook
torch.autocast(device_type="cuda", dtype=torch.bfloat16).__enter__()

if torch.cuda.get_device_properties(0).major >= 8:
    # turn on tfloat32 for Ampere GPUs (https://pytorch.org/docs/stable/notes/cuda.html#tensorfloat-32-tf32-on-ampere-devices)
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.backends.cudnn.allow_tf32 = True

# torch.multiprocessing.set_start_method("spawn")

import colorsys
import datetime
import os
import subprocess

import cv2
import gradio as gr
import imageio.v2 as iio
import numpy as np

from loguru import logger as guru

from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

from sam2.build_sam import build_sam2_video_predictor



class PromptGUI(object):
    def __init__(self, checkpoint_dir, model_cfg):
        self.checkpoint_dir = checkpoint_dir
        self.model_cfg = model_cfg
        self.sam_model = None
        self.tracker = None

        self.selected_points = []
        self.selected_labels = []
        self.cur_label_val = 1.0

        self.frame_index = 0
        self.image = None
        self.cur_mask_idx = 0
        # can store multiple object masks
        # saves the masks and logits for each mask index
        self.cur_masks = {}
        self.cur_logits = {}
        self.index_masks_all = []
        self.color_masks_all = []

        self.img_dir = ""
        self.img_paths = []
        self.init_sam_model()

    def init_sam_model(self):
        if self.sam_model is None:
            self.sam_model = build_sam2_video_predictor(self.model_cfg, self.checkpoint_dir)
            guru.info(f"loaded model checkpoint {self.checkpoint_dir}")


    def clear_points(self) -> tuple[None, None, str]:
        self.selected_points.clear()
        self.selected_labels.clear()
        message = "Cleared points, select new points to update mask"
        return None, None, message

    def add_new_mask(self):
        self.cur_mask_idx += 1
        self.clear_points()
        message = f"Creating new mask with index {self.cur_mask_idx}"
        return None, message

    def make_index_mask(self, masks):
        assert len(masks) > 0
        idcs = list(masks.keys())
        idx_mask = masks[idcs[0]].astype("uint8")
        for i in idcs:
            mask = masks[i]
            idx_mask[mask] = i + 1
        return idx_mask

    def _clear_image(self):
        """
        clears image and all masks/logits for that image
        """
        self.image = None
        self.cur_mask_idx = 0
        self.frame_index = 0
        self.cur_masks = {}
        self.cur_logits = {}
        self.index_masks_all = []
        self.color_masks_all = []

    def reset(self):
        self._clear_image()
        self.sam_model.reset_state(self.inference_state)

    def set_img_dir(self, img_dir: str) -> int:
        self._clear_image()
        self.img_dir = img_dir
        self.img_paths = [
            f"{img_dir}/{p}" for p in sorted(os.listdir(img_dir)) if isimage(p)
        ]
        
        return len(self.img_paths)

    def set_input_image(self, i: int = 0) -> np.ndarray | None:
        guru.debug(f"Setting frame {i} / {len(self.img_paths)}")
        if i < 0 or i >= len(self.img_paths):
            return self.image
        self.clear_points()
        self.frame_index = i
        image = iio.imread(self.img_paths[i])
        self.image = image

        return image

    def get_sam_features(self) -> tuple[str, np.ndarray | None]:
        self.inference_state = self.sam_model.init_state(
            video_path=self.img_dir,
            offload_video_to_cpu=True,
            offload_state_to_cpu=True,
            async_loading_frames=True
        )
        self.sam_model.reset_state(self.inference_state)
        msg = (
            "SAM features extracted. "
            "Click points to update mask, and submit when ready to start tracking"
        )
        return msg, self.image

    def set_positive(self) -> str:
        self.cur_label_val = 1.0
        return "Selecting positive points. Submit the mask to start tracking"

    def set_negative(self) -> str:
        self.cur_label_val = 0.0
        return "Selecting negative points. Submit the mask to start tracking"

    def add_point(self, frame_idx, i, j):
        """
        get the index mask of the objects
        """
        self.selected_points.append([j, i])
        self.selected_labels.append(self.cur_label_val)
        # masks, scores, logits if we want to update the mask
        masks = self.get_sam_mask(
            frame_idx, np.array(self.selected_points, dtype=np.float32), np.array(self.selected_labels, dtype=np.int32)
        )
        mask = self.make_index_mask(masks)

        return mask
    

    def get_sam_mask(self, frame_idx, input_points, input_labels):
        """
        :param frame_idx int
        :param input_points (np array) (N, 2)
        :param input_labels (np array) (N,)
        return (H, W) mask, (H, W) logits
        """
        assert self.sam_model is not None
        

        if torch.cuda.get_device_properties(0).major >= 8:
            # turn on tfloat32 for Ampere GPUs (https://pytorch.org/docs/stable/notes/cuda.html#tensorfloat-32-tf32-on-ampere-devices)
            torch.backends.cuda.matmul.allow_tf32 = True
            torch.backends.cudnn.allow_tf32 = True
        
        with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
            _, out_obj_ids, out_mask_logits = self.sam_model.add_new_points_or_box(
                inference_state=self.inference_state,
                frame_idx=frame_idx,
                obj_id=self.cur_mask_idx,
                points=input_points,
                labels=input_labels,
            )

        return  {
                out_obj_id: (out_mask_logits[i] > 0.0).squeeze().cpu().numpy()
                for i, out_obj_id in enumerate(out_obj_ids)
            }


    def run_tracker(self) -> tuple[str, str]:

        # read images and drop the alpha channel
        images = [iio.imread(p)[:, :, :3] for p in self.img_paths]
        
        video_segments = {}  # video_segments contains the per-frame segmentation results
        if torch.cuda.get_device_properties(0).major >= 8:
            # turn on tfloat32 for Ampere GPUs (https://pytorch.org/docs/stable/notes/cuda.html#tensorfloat-32-tf32-on-ampere-devices)
            torch.backends.cuda.matmul.allow_tf32 = True
            torch.backends.cudnn.allow_tf32 = True
        
        with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
            for out_frame_idx, out_obj_ids, out_mask_logits in self.sam_model.propagate_in_video(self.inference_state, start_frame_idx=0):
                masks = {
                    out_obj_id: (out_mask_logits[i] > 0.0).squeeze().cpu().numpy()
                    for i, out_obj_id in enumerate(out_obj_ids)
                }
                video_segments[out_frame_idx] = masks
                # don't run out of RAM
                del out_mask_logits
                torch.cuda.empty_cache()
            # index_masks_all.append(self.make_index_mask(masks))

        self.index_masks_all = [self.make_index_mask(v) for k, v in video_segments.items()]

        max_num_images = 100
        skip = max(len(images) // max_num_images, 1)
        out_frames, self.color_masks_all = colorize_masks(images, self.index_masks_all)
        out_vidpath = "tracked_colors.mp4"
        print("iio.mimwrite", out_vidpath)
        iio.mimwrite(out_vidpath, out_frames[::skip])
        message = f"Wrote current tracked video to {out_vidpath}."
        instruct = "Save the masks to an output directory if it looks good!"
        return out_vidpath, f"{message} {instruct}"

    def save_masks_to_dir(self, output_dir: str, max_workers=None) -> str:
        assert self.color_masks_all is not None
        os.makedirs(output_dir, exist_ok=True)
        
        def process_and_save_mask(args):
            img_path, clr_mask = args
            name = os.path.basename(img_path) + ".png"
            out_path = f"{output_dir}/{name}"
            
            clr_mask = (clr_mask[:,:,0:1] != 0) | (clr_mask[:,:,1:2] != 0) | (clr_mask[:,:,2:3] != 0)
            clr_mask = (255*clr_mask).astype(np.uint8)
            
            # iio.imwrite(out_path, clr_mask)
            cv2.imwrite(out_path, clr_mask)
            return out_path
        
        args_list = list(zip(self.img_paths, self.color_masks_all))
        saved_paths = []
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Use tqdm to track progress
            for result in tqdm(
                executor.map(process_and_save_mask, args_list),
                total=len(args_list),
                desc="Saving masks"
            ):
                saved_paths.append(result)
        
        message = f"Saved {len(saved_paths)} masks to {output_dir}!"
        guru.debug(message)
        return message

def isimage(p):
    ext = os.path.splitext(p.lower())[-1]
    return ext in [".png", ".jpg", ".jpeg"]


def draw_points(img, points, labels):
    out = img.copy()
    for p, label in zip(points, labels):
        x, y = int(p[0]), int(p[1])
        color = (0, 255, 0) if label == 1.0 else (255, 0, 0)
        out = cv2.circle(out, (x, y), 10, color, -1)
    return out


def get_hls_palette(
    n_colors: int,
    lightness: float = 0.5,
    saturation: float = 0.7,
) -> np.ndarray:
    """
    returns (n_colors, 3) tensor of colors,
        first is black and the rest are evenly spaced in HLS space
    """
    hues = np.linspace(0, 1, int(n_colors) + 1)[1:-1]  # (n_colors - 1)
    # hues = (hues + first_hue) % 1
    palette = [(0.0, 0.0, 0.0)] + [
        colorsys.hls_to_rgb(h_i, lightness, saturation) for h_i in hues
    ]
    return (255 * np.asarray(palette)).astype("uint8")


def colorize_masks(images, index_masks, fac: float = 0.5, max_workers=None):
    max_idx = max([m.max() for m in index_masks])
    guru.debug(f"{max_idx=}")

    palette = get_hls_palette(max_idx + 1)
    color_masks = [None] * len(images)
    out_frames = [None] * len(images)
    
    def process_single_image(args):
        idx, img, mask = args
        clr_mask = palette[mask.astype("int")]
        out_u = compose_img_mask(img, clr_mask, fac)
        return idx, clr_mask, out_u
    
    args_list = [(i, img, mask) for i, (img, mask) in enumerate(zip(images, index_masks))]
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for idx, clr_mask, out_u in tqdm(
            executor.map(process_single_image, args_list),
            total=len(images),
            desc="Colorizing masks"
        ):
            color_masks[idx] = clr_mask
            out_frames[idx] = out_u
    
    return out_frames, color_masks


def compose_img_mask(img, color_mask, fac: float = 0.5):
    out_f = fac * img / 255 + (1 - fac) * color_mask / 255
    out_u = (255 * out_f).astype("uint8")
    return out_u


def listdir(vid_dir):
    if vid_dir is not None and os.path.isdir(vid_dir):
        return sorted(os.listdir(vid_dir))
    return []


def make_demo(
    checkpoint_dir,
    model_cfg,
    root_dir,
    vid_name: str = "videos",
    img_name: str = "images",
    mask_name: str = "masks",
):
    prompts = PromptGUI(checkpoint_dir, model_cfg)

    start_instructions = (
        "Select a video file to extract frames from, "
        "or select an image directory with frames already extracted."
    )
    vid_root, img_root = (f"{root_dir}/{vid_name}", f"{root_dir}/{img_name}")
    with gr.Blocks() as demo:
        instruction = gr.Textbox(
            start_instructions, label="Instruction", interactive=False
        )
        with gr.Row():
            root_dir_field = gr.Text(root_dir, label="Dataset root directory")
            vid_name_field = gr.Text(vid_name, label="Video subdirectory name")
            img_name_field = gr.Text(img_name, label="Image subdirectory name")
            mask_name_field = gr.Text(mask_name, label="Mask subdirectory name")
            seq_name_field = gr.Text(None, label="Sequence name", interactive=False)

        with gr.Row():
            with gr.Column():
                vid_files = listdir(vid_root)
                vid_files_field = gr.Dropdown(label="Video files", choices=vid_files)
                input_video_field = gr.Video(label="Input Video")

                with gr.Row():
                    start_time = gr.Number(0, label="Start time (s)")
                    end_time = gr.Number(0, label="End time (s)")
                    sel_fps = gr.Number(30, label="FPS")
                    sel_height = gr.Number(540, label="Height")
                    extract_button = gr.Button("Extract frames")

            with gr.Column():
                img_dirs = listdir(img_root)
                img_dirs = [img_root]
                img_dirs_field = gr.Dropdown(
                    label="Image directories", choices=img_dirs
                )
                img_dir_field = gr.Text(
                    None, label="Input directory", interactive=False
                )
                frame_index = gr.Slider(
                    label="Frame index",
                    minimum=0,
                    maximum=len(prompts.img_paths) - 1,
                    value=0,
                    step=1,
                )
                sam_button = gr.Button("Get SAM features")
                reset_button = gr.Button("Reset")
                input_image = gr.Image(
                    prompts.set_input_image(0),
                    label="Input Frame",
                    every=1,
                )
                with gr.Row():
                    pos_button = gr.Button("Toggle positive")
                    neg_button = gr.Button("Toggle negative")
                clear_button = gr.Button("Clear points")

            with gr.Column():
                output_img = gr.Image(label="Current selection")
                add_button = gr.Button("Add new mask")
                submit_button = gr.Button("Submit mask for tracking")
                final_video = gr.Video(label="Masked video")
                mask_dir_field = gr.Text(
                    None, label="Path to save masks", interactive=False
                )
                save_button = gr.Button("Save masks")

        def update_vid_root(root_dir, vid_name):
            vid_root = f"{root_dir}/{vid_name}"
            vid_paths = listdir(vid_root)
            guru.debug(f"Updating video paths: {vid_paths=}")
            return vid_paths

        def update_img_root(root_dir, img_name):
            img_root = f"{root_dir}/{img_name}"
            img_dirs = listdir(img_root)
            guru.debug(f"Updating img dirs: {img_dirs=}")
            return img_root, img_dirs

        def update_mask_dir(root_dir, mask_name, seq_name):
            # return f"{root_dir}/{mask_name}/{seq_name}"
            return seq_name.replace("images", "masks")

        def update_root_paths(root_dir, vid_name, img_name, mask_name, seq_name):
            return (
                update_vid_root(root_dir, vid_name),
                update_img_root(root_dir, img_name),
                update_mask_dir(root_dir, mask_name, seq_name),
            )

        def select_video(root_dir, vid_name, seq_file):
            seq_name = os.path.splitext(seq_file)[0]
            guru.debug(f"Selected video: {seq_file=}")
            vid_path = f"{root_dir}/{vid_name}/{seq_file}"
            return seq_name, vid_path

        def extract_frames(
            root_dir, vid_name, img_name, vid_file, start, end, fps, height, ext="png"
        ):
            seq_name = os.path.splitext(vid_file)[0]
            vid_path = f"{root_dir}/{vid_name}/{vid_file}"
            # out_dir = f"{root_dir}/{img_name}/{seq_name}"
            out_dir = seq_name
            guru.debug(f"Extracting frames to {out_dir}")
            os.makedirs(out_dir, exist_ok=True)

            def make_time(seconds):
                return datetime.time(
                    seconds // 3600, (seconds % 3600) // 60, seconds % 60
                )

            start_time = make_time(start).strftime("%H:%M:%S")
            end_time = make_time(end).strftime("%H:%M:%S")
            cmd = (
                f"ffmpeg -ss {start_time} -to {end_time} -i {vid_path} "
                f"-vf 'scale=-1:{height},fps={fps}' {out_dir}/%05d.{ext}"
            )
            print(cmd)
            subprocess.call(cmd, shell=True)
            img_root = f"{root_dir}/{img_name}"
            img_dirs = listdir(img_root)
            return out_dir, img_dirs

        def select_image_dir(root_dir, img_name, seq_name):
            # img_dir = f"{root_dir}/{img_name}/{seq_name}"
            img_dir = seq_name
            guru.debug(f"Selected image dir: {img_dir}")
            return seq_name, img_dir

        def update_image_dir(root_dir, img_name, seq_name):
            # img_dir = f"{root_dir}/{img_name}/{seq_name}"
            img_dir = seq_name
            num_imgs = prompts.set_img_dir(img_dir)
            slider = gr.Slider(minimum=0, maximum=num_imgs - 1, value=0, step=1)
            message = (
                f"Loaded {num_imgs} images from {img_dir}. Choose a frame to run SAM!"
            )
            return slider, message

        def get_select_coords(frame_idx, img, evt: gr.SelectData):
            i = evt.index[1]  # type: ignore
            j = evt.index[0]  # type: ignore
            index_mask = prompts.add_point(frame_idx, i, j)
            guru.debug(f"{index_mask.shape=}")
            palette = get_hls_palette(index_mask.max() + 1)
            color_mask = palette[index_mask]
            out_u = compose_img_mask(img, color_mask)
            out = draw_points(out_u, prompts.selected_points, prompts.selected_labels)
            return out

        # update the root directory
        # and associated video, image, and mask root directories
        root_dir_field.submit(
            update_root_paths,
            [
                root_dir_field,
                vid_name_field,
                img_name_field,
                mask_name_field,
                seq_name_field,
            ],
            outputs=[vid_files_field, img_dirs_field, mask_dir_field],
        )
        vid_name_field.submit(
            update_vid_root,
            [root_dir_field, vid_name_field],
            outputs=[vid_files_field],
        )
        img_name_field.submit(
            update_img_root,
            [root_dir_field, img_name_field],
            outputs=[img_dirs_field],
        )
        mask_name_field.submit(
            update_mask_dir,
            [root_dir_field, mask_name_field, seq_name_field],
            outputs=[mask_dir_field],
        )

        # selecting a video file
        vid_files_field.select(
            select_video,
            [root_dir_field, vid_name_field, vid_files_field],
            outputs=[seq_name_field, input_video_field],
        )

        # when the img_dir_field changes
        img_dir_field.change(
            update_image_dir,
            [root_dir_field, img_name_field, seq_name_field],
            [frame_index, instruction],
        )
        seq_name_field.change(
            update_mask_dir,
            [root_dir_field, mask_name_field, seq_name_field],
            outputs=[mask_dir_field],
        )

        # selecting an image directory
        img_dirs_field.select(
            select_image_dir,
            [root_dir_field, img_name_field, img_dirs_field],
            [seq_name_field, img_dir_field],
        )

        # extracting frames from video
        extract_button.click(
            extract_frames,
            [
                root_dir_field,
                vid_name_field,
                img_name_field,
                vid_files_field,
                start_time,
                end_time,
                sel_fps,
                sel_height,
            ],
            outputs=[img_dir_field, img_dirs_field],
        )

        frame_index.change(prompts.set_input_image, [frame_index], [input_image])
        input_image.select(get_select_coords, [frame_index, input_image], [output_img])

        sam_button.click(prompts.get_sam_features, outputs=[instruction, input_image])
        reset_button.click(prompts.reset)
        clear_button.click(
            prompts.clear_points, outputs=[output_img, final_video, instruction]
        )
        pos_button.click(prompts.set_positive, outputs=[instruction])
        neg_button.click(prompts.set_negative, outputs=[instruction])

        add_button.click(prompts.add_new_mask, outputs=[output_img, instruction])
        submit_button.click(prompts.run_tracker, outputs=[final_video, instruction])
        save_button.click(
            prompts.save_masks_to_dir, [mask_dir_field], outputs=[instruction]
        )

    return demo


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=8890)
    parser.add_argument("--checkpoint_dir", type=str, default="checkpoints/sam2.1_hiera_large.pt")
    parser.add_argument("--model_cfg", type=str, default="configs/sam2.1/sam2.1_hiera_l.yaml")
    parser.add_argument("--root_dir", type=str, required=True)
    parser.add_argument("--vid_name", type=str, default="videos")
    parser.add_argument("--img_name", type=str, default="images")
    parser.add_argument("--mask_name", type=str, default="masks")
    args = parser.parse_args()

    # device = "cuda" if torch.cuda.is_available() else "cpu"

    demo = make_demo(
        args.checkpoint_dir,
        args.model_cfg,
        args.root_dir,
        args.vid_name,
        args.img_name,
        args.mask_name
    )
    demo.launch(server_port=args.port)
