#!/usr/bin/env python3

# extract frames from video or rosbag file
# pick good frames that minimize motion blur


import cv2
import numpy as np
import os

import subprocess

quality = 95
rotate = 0
scale = 1.0

ENHANCE = False  # unsupported for now


def write_image(image, filename, filename_enh=None, _enhance_model=[]):
    if rotate != 0:
        rotate_param = {
            90: cv2.ROTATE_90_CLOCKWISE,
            180: cv2.ROTATE_180,
            -90: cv2.ROTATE_90_COUNTERCLOCKWISE,
            270: cv2.ROTATE_90_COUNTERCLOCKWISE
        }[rotate]
        image = cv2.rotate(image, rotate_param)
    if image is None:
        return

    # resize image
    tile = 8 if ENHANCE else 1
    if scale < 1.0:
        h, w = image.shape[:2]
        image = cv2.resize(image, (int(scale*w/tile+0.5)*tile, int(scale*h/tile+0.5)*tile),
                           interpolation=cv2.INTER_AREA)
    elif scale == int(scale) and scale != 1.0:
        iscale = int(scale)
        h, w, c = image.shape
        h, w = (h//(iscale*tile))*tile, (w//(iscale*tile))*tile
        image = image[:h*iscale, :w*iscale]
        image = np.mean(np.mean(
            image.reshape(h, iscale, w, iscale, c),
            axis=1), axis=2).astype(np.uint8)

    if image is None:
        return

    def write_im(image, filename):
        if 0 <= quality <= 100:
            filename += '.jpg'
            encode_param = [int(cv2.IMWRITE_JPEG_QUALITY), quality]
            cv2.imwrite(filename, image, encode_param)
        else:
            filename += '.png'
            cv2.imwrite(filename, image)

    write_im(image, filename)

    if filename_enh is not None:
        raise NotImplementedError()
        if len(_enhance_model) == 0:
            import colmap.image_enhance as enhance
            image_enhance_model = enhance.load_image_enhance_model_02()
            image_enhance_model = image_enhance_model.float()
            _enhance_model.append(enhance.torch)
            _enhance_model.append(image_enhance_model)
            torch = enhance.torch
        else:
            torch, image_enhance_model = _enhance_model

        image = torch.from_numpy(image).cuda().permute(2,0,1).unsqueeze(0).float()
        image = image_enhance_model(image)
        image = torch.clip(image, 0.0, 255.0)[0].byte().permute(1,2,0).cpu().numpy()
        write_im(image, filename_enh)

    return filename


class FrameSelector:
    def __init__(self, image_dir, enhance_dir, max_frames, skip, keep):
        os.makedirs(image_dir, exist_ok=True)
        if enhance_dir is not None:
            os.makedirs(enhance_dir, exist_ok=True)

        self.image_dir = image_dir
        self.enhance_dir = enhance_dir
        self.max_frames = max_frames
        self.skip = skip
        self.keep = keep

        self.frame_count = 0
        self.written_frame_count = 0
        self.iqms = []
        self.frames_buf = []

    def add_frame(self, frame, comment=""):
        if self.written_frame_count >= self.max_frames:
            print("Maximum number of frames reached.")
            return False
        
        # check motion blur
        s = 256
        gray = cv2.cvtColor(cv2.resize(frame, (2*s, 2*s)), cv2.COLOR_BGR2GRAY)
        mu, sigma = np.mean(gray), np.std(gray)
        gray = (gray-mu)#/sigma
        
        iqm = cv2.Laplacian(gray, cv2.CV_64F).var()

        self.iqms.append((mu, sigma, iqm))

        self.frames_buf.append((self.frame_count, iqm, frame, comment))
        if len(self.frames_buf) > max(self.keep, 1):
            del self.frames_buf[0]

        if self.keep != 0:
            self.frame_count += 1
        if self.frame_count % self.skip == 0:
            fid, iqm, frame, comment = sorted(self.frames_buf, key=lambda x: x[1], reverse=True)[0]
            filename = os.path.join(self.image_dir, f"{fid:05d}")
            filename_enh = None if self.enhance_dir is None else os.path.join(self.enhance_dir, f"{fid:05d}")
            if False:  # gamma adjustment
                gamma = 0.6
                lut = (255*np.linspace(0, 1, 256)**gamma).astype(np.uint8)
                frame = lut[frame]
            filename = write_image(frame, filename, filename_enh)
            print(filename, comment)
            self.written_frame_count += 1
        if self.keep == 0:
            self.frame_count += 1

        return True

    def conclude(self):
        print(f"Extracted {self.written_frame_count} frames.")



def extract_rosbag_frames(bag_path, image_dir, enhance_dir, max_frames, skip, keep, required_topic=None):

    def _imgmsg_to_cv2(img_msg):
        
        dtype = np.dtype("uint8")
        dtype = dtype.newbyteorder('>' if img_msg.is_bigendian else '<')
        image_opencv = np.ndarray(shape=(img_msg.height, img_msg.width, 3),
                                dtype=dtype, buffer=img_msg.data)
        
        if img_msg.encoding == "rgb8":
            image_opencv = cv2.cvtColor(image_opencv, cv2.COLOR_BGR2RGB)
        elif img_msg.encoding != "bgr8":
            raise ValueError(f"Unsupported encoding: {img_msg.encoding}")
        
        return image_opencv

    def _compressed_imgmsg_to_cv2(compressed_img_msg):
        np_arr = np.frombuffer(compressed_img_msg.data, np.uint8)
        return cv2.imdecode(np_arr, cv2.IMREAD_COLOR)

    def get_image(msg):
        if msg._type == 'sensor_msgs/Image':
            cv_image = _imgmsg_to_cv2(msg)
        else:  # CompressedImage
            cv_image = _compressed_imgmsg_to_cv2(msg)
        return cv_image

    frame_selector = FrameSelector(image_dir, enhance_dir, max_frames, skip, keep)

    all_image_topics = set()
    frame_count = 0
    with __import__('rosbag').Bag(bag_path, 'r') as bag:
        for topic, msg, t in bag.read_messages():
            if required_topic is not None and topic != required_topic:
                continue
            if msg._type in ['sensor_msgs/Image', 'sensor_msgs/CompressedImage']:
                all_image_topics.add(topic)
                if frame_count % skip != 0:
                    frame_count += 1
                    continue
                image = get_image(msg)
                
                if not frame_selector.add_frame(image, 'topic='+topic):
                    break

    frame_selector.conclude()

    if required_topic is None:
        print("Detected image topics:", all_image_topics)


def extract_video_frames(video_path, image_dir, enhance_dir, max_frames, skip, keep):

    def get_video_streams(video_path):
        cmd = [
            "ffprobe", "-v", "error", "-select_streams", "v",
            "-show_entries", "stream=index", "-of", "csv=p=0", video_path
        ]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
        except FileNotFoundError:
            print("ffprobe not found. Please install ffmpeg.")
            exit(0)
        if result.returncode != 0:
            raise RuntimeError(f"ffprobe failed: {result.stderr}")
        streams = [int(x.strip()) for x in result.stdout.splitlines() if x.strip()]
        return streams

    streams = get_video_streams(video_path) if video_path.endswith(".insv") else [0]

    if len(streams) == 0:
        print(f"Error: no video streams found in {video_path}")
        return

    for idx, stream_idx in enumerate(streams):
        # Handle multi-track video by splitting track via ffmpeg if needed
        if len(streams) > 1:
            print(f"Extracting track {stream_idx} to cam{idx}")
            temp_path = os.path.join(image_dir, f"temp_cam{idx}.mp4")
            subprocess.run([
                "ffmpeg", "-y", "-i", video_path,
                "-map", f"0:v:{idx}", "-c", "copy", temp_path
            ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            track_video_path = temp_path
            out_image_dir = os.path.join(image_dir, f"cam{idx}")
            out_enhance_dir = os.path.join(enhance_dir, f"cam{idx}") if enhance_dir is not None else None
        else:
            track_video_path = video_path
            out_image_dir = image_dir
            out_enhance_dir = enhance_dir

        video = cv2.VideoCapture(track_video_path)
        if not video.isOpened():
            print(f"Error: could not open {track_video_path}")
            continue

        frame_selector = FrameSelector(out_image_dir, out_enhance_dir, max_frames, skip, keep)

        while True:
            success, frame = video.read()
            if not success:
                break
            if not frame_selector.add_frame(frame):
                break

        frame_selector.conclude()
        video.release()

        if len(streams) > 1:
            os.remove(track_video_path)

    print("Frame extraction completed.")


def extract_frames(filename, dirname, max_frames=100000, skip=1, keep=-1, ros_topic=""):
    if ENHANCE:
        image_dir = os.path.join(dirname, 'images_raw')
        enhance_dir = os.path.join(dirname, 'images')
    else:
        image_dir = os.path.join(dirname, 'images')
        enhance_dir = None

    if not os.path.exists(image_dir):
        os.makedirs(image_dir)
    if enhance_dir is not None and not os.path.exists(enhance_dir):
        os.makedirs(enhance_dir)

    if keep == -1:
        keep = int(0.5*skip+0.5)

    if filename.split('.')[-1].lower() == 'bag':
        if ros_topic == "":
            ros_topic = None
        extract_rosbag_frames(filename, image_dir, enhance_dir, max_frames, skip, keep, ros_topic)

    else:  # video
        extract_video_frames(filename, image_dir, enhance_dir, max_frames, skip, keep)



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Extract frames from a rosbag or a video.")
    parser.add_argument("input_file", nargs=1, help="The input rosbag or video file.")
    parser.add_argument("--max_frames", "-n", default="100000", help="Maximum number of frames to output.")
    parser.add_argument("--skip", "-s", default="1", help="Take one frame every this number of frames.")
    parser.add_argument("--keep", "-k", default="-1", help="A fraction of --skip, nonzero for selecting images with minimum motion blur")
    parser.add_argument("--topic", default="", help="rosbag topic. If not set, it will exact all detected frames.")
    parser.add_argument("--quality", "-q", default="95", help="Quality to save JPEG images. Save lossless PNG if this is not between 0 and 100.")
    parser.add_argument("--rotate", "-r", default="0", help="Rotate clockwise this degrees, must be one of {0, 90, 180, 270}.")
    parser.add_argument("--scale", default="1", help="Relative image size scale, less than 1.0 for downscale.")
    parser.add_argument("--mask", default="0", help="Whether to use export scripts to generate masks.")
    args = parser.parse_args()

    quality = int(args.quality)
    rotate = int(args.rotate)
    scale = float(args.scale)
    filename = args.input_file[0]
    dirname = filename[:filename.rfind('.')]
    extract_frames(
        filename, dirname,
        int(args.max_frames), int(args.skip), int(args.keep), args.topic
    )

    cur_dir = os.path.dirname(__file__)

    # generate masks
    if int(args.mask) != 0:
        # run_colmap.bash
        open(os.path.join(dirname, 'run_colmap.bash'), 'w').write(f"""
# {' '.join(__import__('sys').argv)}

{open(os.path.join(cur_dir, "run_colmap_mask.bash")).read()}""".lstrip())
        # run_sam2.bash
        open(os.path.join(dirname, 'run_sam2.bash'), 'w').write(f"""
# Run this from SAM2 install directory, after downloading checkpoints
python3 {os.path.join(cur_dir, "SAM2-GUI", "mask_app.py")} --root_dir {os.path.abspath(dirname)}
""".lstrip())
        # run_lang_sam.bash

    # no mask
    else:
        # run_colmap.bash
        open(os.path.join(dirname, 'run_colmap.bash'), 'w').write(f"""
# {' '.join(__import__('sys').argv)}

{open(os.path.join(cur_dir, "run_colmap.bash")).read()}""".lstrip())
