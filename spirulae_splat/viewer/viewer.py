import sys
import os
import numpy as np
from scipy.spatial.transform import Rotation as R

from PyQt6.QtWidgets import QApplication, QMainWindow, QLabel, QVBoxLayout, QWidget
from PyQt6.QtGui import QImage, QPixmap
from PyQt6.QtCore import Qt, QTimer
import matplotlib

import torch
import torch.nn.functional as F
from spirulae_splat.viewer import Camera, SplatModel
from spirulae_splat.splat._torch_impl import depth_inv_map
from spirulae_splat.splat import depth_to_normal


class RenderViewer(QMainWindow):
    colormap = matplotlib.colormaps['magma']
    vis_modes = ["rgb", "depth", "normal", "shaded"]
    move_mode = ["pinch", "navigate"][0]

    def __init__(self):
        super().__init__()
        self.setWindowTitle("spirulae-splat Viewer")
        self.setGeometry(100, 100, 512, 512)
        
        self.image_label = QLabel(self)
        self.image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        layout = QVBoxLayout()
        layout.addWidget(self.image_label)
        
        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)
        
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
        self.update_image()

        self.timer = QTimer(self)
        self.timer.timeout.connect(self.animation_frame)
        self.timer.start(16)  # 60 fps

    def color_depth(self, depth):
        depth = depth[..., 0] / depth.max()
        colored_data = self.colormap(depth)
        return (colored_data[:, :, :3] * 255).astype(np.uint8)

    def update_image(self):
        image, depth = ssplat_model.render(camera, self.c2w, True)
        self.last_depth = depth_inv_map(depth)
        if self.vis_modes[self.vis_mode] == "rgb":
            image = (255*image.clip(0, 1)).byte().cpu().numpy()
        if self.vis_modes[self.vis_mode] == "depth":
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
        h, w, c = image.shape
        qimage = QImage(image.data, w, h, 3 * w, QImage.Format.Format_RGB888)
        self.image_label.setPixmap(QPixmap.fromImage(qimage))

    def resizeEvent(self, event):
        return
        w = event.size().width()
        h = event.size().height()
        w0 = event.oldSize().width()
        h0 = event.oldSize().height()
        if (w != w0 or h != h0) and (w0 > 0 and h0 > 0):
            ssplat_camera.resize(
                w-w0+ssplat_camera.w,
                h-h0+ssplat_camera.h
            )
            self.update_image()

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.last_mouse_pos = event.pos()

            # get position
            x, y = event.pos().x(), event.pos().y()
            ssplat_camera = camera._to_ssplat_camera()
            if ssplat_camera.is_distorted():
                pass
            else:
                fx, fy, cx, cy = ssplat_camera.intrins
                if False:
                    depth = self.last_depth[int(y+0.5), int(x+0.5), 0].item()
                    self.click_pos_3d = np.array([(x-cx)/fx, (y-cy)/fy, 1]) * depth
                elif True:
                    # depth = torch.median(self.last_depth).item()
                    # depth = 2.0
                    depth = 0.0
                    self.click_pos_3d = np.array([0, 0, 1]) * depth

    def mouseMoveEvent(self, event):
        if self.last_mouse_pos is not None:
            dx = (event.pos().x() - self.last_mouse_pos.x()) * 0.01
            dy = (event.pos().y() - self.last_mouse_pos.y()) * 0.01
            matR = self._mat_rotate(-dy, -dx, 0.0)
            if self.move_mode == "pinch":
                w2c = np.linalg.inv(self.c2w)
                matT = self._mat_translate(*w2c[:3, 3])
            else:
                matT = self._mat_translate(*self.click_pos_3d)
            mat = matT @ matR @ np.linalg.inv(matT)
            self.c2w = self.c2w @ mat
            self.last_mouse_pos = event.pos()
            self.update_image()

    def wheelEvent(self, event):
        zoom_factor = event.angleDelta().y() / 2000
        self.c2w[:3, 3] += zoom_factor * self.c2w[:3, 2]
        self.update_image()

    def keyPressEvent(self, event):
        self.active_keys.add(event.key())

        # print(event.key(), self.active_keys)
        # print(Qt.Key.Key_Control.value)

        # space to switch between RGB and depth
        if event.key() == ord(' '):
            self.vis_mode = (self.vis_mode + 1) % len(self.vis_modes)
            self.update_image()

        # Ctrl+C or Esc to exit
        if (
            (event.key() == ord('C') and Qt.Key.Key_Control.value in self.active_keys)
            or event.key() == Qt.Key.Key_Escape.value
        ):
            app.exit()

    def keyReleaseEvent(self, event):
        self.active_keys.remove(event.key())
        pass

    def animation_frame(self):
        # WSAD for rotation
        mat = np.eye(4)
        rotate_speed = 0.025
        move_speed = 0.02
        if self.move_mode == "pinch":
            if ord('W') in self.active_keys:
                mat = self._mat_rotate(rotate_speed, 0, 0)
            elif ord('S') in self.active_keys:
                mat = self._mat_rotate(-rotate_speed, 0, 0)
            elif ord('A') in self.active_keys:
                mat = self._mat_rotate(0, -rotate_speed, 0)
            elif ord('D') in self.active_keys:
                mat = self._mat_rotate(0, rotate_speed, 0)
        elif self.move_mode == "navigate":
            if ord('W') in self.active_keys:
                mat = self._mat_translate(0, 0, move_speed)
            elif ord('S') in self.active_keys:
                mat = self._mat_translate(0, 0, -move_speed)
            elif ord('A') in self.active_keys:
                mat = self._mat_translate(-move_speed, 0, 0)
            elif ord('D') in self.active_keys:
                mat = self._mat_translate(move_speed, 0, 0)
        if ord('Q') in self.active_keys:
            mat = self._mat_rotate(0, 0, rotate_speed)
        elif ord('E') in self.active_keys:
            mat = self._mat_rotate(0, 0, -rotate_speed)

        if (mat != np.eye(4)).any():
            self.c2w = self.c2w @ mat
            self.update_image()

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

    camera_path = os.path.join(os.path.dirname(__file__), "cameras/s21.yaml")
    camera = Camera(camera_path)

    if len(sys.argv) > 1:
        file_path = sys.argv[1]
    else:
        print("Usage: python3 path/to/viewer.py path/to/config.yml")
        exit(-1)

    ssplat_model = SplatModel(file_path)
    ssplat_model.bgr = False
    ssplat_model.return_torch = True
    print(ssplat_model.num_splats(), "splats")
    print("sort_per_pixel:", ssplat_model.sort_per_pixel)

    ssplat_model.gauss_params["means"] -= torch.mean(
        torch.nan_to_num(ssplat_model.means, 0.0, 0.0, 0.0),
        dim=0, keepdim=True)

    app = QApplication(sys.argv)
    viewer = RenderViewer()
    viewer.show()
    sys.exit(app.exec())
