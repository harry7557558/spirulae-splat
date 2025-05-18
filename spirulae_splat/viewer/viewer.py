import sys
import os
import numpy as np
from scipy.spatial.transform import Rotation as R

import tkinter as tk
from tkinter import ttk
from PIL import Image, ImageTk
import matplotlib

import torch
import torch.nn.functional as F
from spirulae_splat.viewer import Camera, SplatModel
from spirulae_splat.splat._torch_impl import depth_inv_map
from spirulae_splat.splat import depth_to_normal

import time


class RenderViewer(tk.Tk):
    colormap = matplotlib.colormaps['magma']
    vis_modes = ["rgb", "depth", "normal", "shaded"]
    move_mode = ["pinch", "navigate"][0]

    def __init__(self):
        super().__init__()
        self.title("spirulae-splat Viewer")
        self.geometry(f"{camera.w}x{camera.h}")
        
        # Create image label
        self.image_label = ttk.Label(self)
        self.image_label.pack(expand=True, fill=tk.BOTH)
        
        self.last_mouse_pos = None
        self.click_pos_3d = None
        self.active_keys = set()
        self.vis_mode = 0

        self.c2w = np.array([
            [1, 0, 0, 0],
            [0, 0, 1, 0],
            [0, -1, 0, 0],
            [0, 0, 0, 1]
        ], dtype=np.float32)
        
        # Event bindings
        self.image_label.bind("<ButtonPress-1>", self.on_mouse_press)
        self.image_label.bind("<B1-Motion>", self.on_mouse_move)
        self.bind("<MouseWheel>", self.on_mouse_wheel)  # Windows and macOS
        self.bind("<Button-4>", self.on_mouse_wheel_linux)  # Linux scroll up
        self.bind("<Button-5>", self.on_mouse_wheel_linux)  # Linux scroll down
        self.bind("<KeyPress>", self.on_key_press)
        self.bind("<KeyRelease>", self.on_key_release)
        
        # Make label focusable for keyboard events
        self.image_label.config(takefocus=1)
        self.image_label.focus_set()
        
        self.fps = 60
        self.last_image_update = -1.0
        self.update_image()
        
        # Animation timer
        self.after(1000//self.fps, self.animation_frame)  # 60 fps
    
    def color_depth(self, depth):
        depth = depth[..., 0] / depth.max()
        colored_data = self.colormap(depth)
        return (colored_data[:, :, :3] * 255).astype(np.uint8)

    def update_image(self):
        # Reduce lag while scrolling
        if time.perf_counter() - self.last_image_update < 0.5/self.fps:
            return

        # Render image
        image, depth = ssplat_model.render(camera, self.c2w, True)
        self.last_depth = depth_inv_map(depth)
        
        if self.vis_modes[self.vis_mode] == "rgb":
            image = (255*image.clip(0, 1)).byte().cpu().numpy()
        elif self.vis_modes[self.vis_mode] == "depth":
            depth = depth.cpu().numpy()
            image = self.color_depth(depth)
        elif self.vis_modes[self.vis_mode] == "normal":
            normal = depth_to_normal(depth_inv_map(depth), camera._to_ssplat_camera())
            image = (255*(0.5-0.5*normal)).byte().cpu().numpy()
        elif self.vis_modes[self.vis_mode] == "shaded":
            normal, point = depth_to_normal(depth_inv_map(depth), camera._to_ssplat_camera(), return_points=True)
            ldir = F.normalize(torch.tensor([[[1.0, -1.0, -1.0]]]), dim=-1).to(image)
            rdir = F.normalize(point, dim=-1)
            ambient = 0.3 + 0.02 * normal
            diffuse = torch.relu((normal*ldir).sum(-1, True)) * torch.tensor([[[0.5, 0.55, 0.6]]]).to(image)
            backlit = torch.relu(-(normal*ldir).sum(-1, True)) * torch.tensor([[[0.15, 0.1, 0.1]]]).to(image)
            specular = ((rdir - 2.0*(rdir*normal).sum(-1, True) * normal)*ldir).sum(-1, True).clip(0, 1) ** 5 \
                * torch.tensor([[[0.2, 0.15, 0.1]]]).to(image)
            image = ambient + diffuse - backlit + specular
            image = (255*image.clip(0, 1)).byte().cpu().numpy()
        
        # Convert numpy array to PIL Image and then to ImageTk
        pil_img = Image.fromarray(image)
        if not hasattr(self, 'tk_img'):
            self.tk_img = ImageTk.PhotoImage(image=pil_img)
            self.image_label.config(image=self.tk_img)
        else:
            self.tk_img.paste(pil_img)

        self.last_image_update = time.perf_counter()
        
    def on_mouse_press(self, event):
        self.last_mouse_pos = (event.x, event.y)
        
        # Get position for 3D click
        x, y = event.x, event.y
        ssplat_camera = camera._to_ssplat_camera()
        if ssplat_camera.is_distorted():
            pass
        else:
            fx, fy, cx, cy = ssplat_camera.intrins
            # Use fixed depth for now
            depth = 0.0
            self.click_pos_3d = np.array([0, 0, 1]) * depth
    
    def on_mouse_move(self, event):
        if self.last_mouse_pos:
            dx = (event.x - self.last_mouse_pos[0]) * 0.01
            dy = (event.y - self.last_mouse_pos[1]) * 0.01
            matR = self._mat_rotate(-dy, -dx, 0.0)
            if self.move_mode == "pinch":
                w2c = np.linalg.inv(self.c2w)
                matT = self._mat_translate(*w2c[:3, 3])
            else:
                matT = self._mat_translate(*self.click_pos_3d)
            mat = matT @ matR @ np.linalg.inv(matT)
            self.c2w = self.c2w @ mat
            self.last_mouse_pos = (event.x, event.y)
            self.update_image()
    
    def on_mouse_wheel(self, event):
        # Windows and macOS wheel event
        zoom_factor = event.delta / 2000
        self.c2w[:3, 3] += zoom_factor * self.c2w[:3, 2]
        self.update_image()
    
    def on_mouse_wheel_linux(self, event):
        # Linux scroll events
        if event.num == 4:  # scroll up
            zoom_factor = 0.05
        else:  # scroll down (event.num == 5)
            zoom_factor = -0.05
        self.c2w[:3, 3] += zoom_factor * self.c2w[:3, 2]
        self.update_image()
    
    def on_key_press(self, event):
        self.active_keys.add(event.keysym)
        
        # Space to switch between RGB and depth
        if event.keysym == "space":
            self.vis_mode = (self.vis_mode + 1) % len(self.vis_modes)
            self.update_image()
        
        # Ctrl+C or Escape to exit
        if (event.keysym == "c" and "Control_L" in self.active_keys) or event.keysym == "Escape":
            self.quit()
    
    def on_key_release(self, event):
        if event.keysym in self.active_keys:
            self.active_keys.remove(event.keysym)
    
    def animation_frame(self):
        # Handle WASD rotation and movement
        mat = np.eye(4)
        rotate_speed = 0.025
        move_speed = 0.02
        
        if self.move_mode == "pinch":
            if "w" in self.active_keys:
                mat = self._mat_rotate(rotate_speed, 0, 0)
            elif "s" in self.active_keys:
                mat = self._mat_rotate(-rotate_speed, 0, 0)
            elif "a" in self.active_keys:
                mat = self._mat_rotate(0, -rotate_speed, 0)
            elif "d" in self.active_keys:
                mat = self._mat_rotate(0, rotate_speed, 0)
        elif self.move_mode == "navigate":
            if "w" in self.active_keys:
                mat = self._mat_translate(0, 0, move_speed)
            elif "s" in self.active_keys:
                mat = self._mat_translate(0, 0, -move_speed)
            elif "a" in self.active_keys:
                mat = self._mat_translate(-move_speed, 0, 0)
            elif "d" in self.active_keys:
                mat = self._mat_translate(move_speed, 0, 0)
                
        if "q" in self.active_keys:
            mat = self._mat_rotate(0, 0, rotate_speed)
        elif "e" in self.active_keys:
            mat = self._mat_rotate(0, 0, -rotate_speed)

        if (mat != np.eye(4)).any():
            self.c2w = self.c2w @ mat
            self.update_image()
        
        # Schedule next frame
        self.after(1000//self.fps, self.animation_frame)
        
    @staticmethod
    def _mat3_to_mat4(mat3):
        mat4 = np.eye(4, dtype=mat3.dtype)
        mat4[:3, :3] = mat3
        return mat4

    @staticmethod
    def _mat_translate(dx, dy, dz):
        return np.array([
            [1, 0, 0, dx],
            [0, 1, 0, dy],
            [0, 0, 1, dz],
            [0, 0, 0, 1]
        ])

    @staticmethod
    def _mat_rotate(dx, dy, dz):
        rot = R.from_rotvec([dx, 0, 0]) \
            * R.from_rotvec([0, dy, 0]) \
            * R.from_rotvec([0, 0, dz])
        return RenderViewer._mat3_to_mat4(rot.as_matrix())


if __name__ == "__main__":
    # Parse command line arguments
    if len(sys.argv) > 1:
        file_path = sys.argv[1]
    else:
        print("Usage: python3 path/to/viewer.py path/to/config.yml")
        sys.exit(-1)

    # Load camera
    camera_path = os.path.join(os.path.dirname(__file__), "cameras/s21.yaml")
    camera = Camera(camera_path)

    # Initialize model
    ssplat_model = SplatModel(file_path)
    ssplat_model.bgr = False
    ssplat_model.return_torch = True
    ssplat_model.convert_to_input_frame()
    print(ssplat_model.num_splats(), "splats")
    print("sort_per_pixel:", ssplat_model.sort_per_pixel)

    # Center the model
    if False:
        ssplat_model.gauss_params["means"] -= torch.mean(
            torch.nan_to_num(ssplat_model.means, 0.0, 0.0, 0.0),
            dim=0, keepdim=True)

    # Start tkinter app
    viewer = RenderViewer()
    viewer.mainloop()
