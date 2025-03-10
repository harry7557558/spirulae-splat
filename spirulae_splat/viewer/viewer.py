import sys
import os
import numpy as np
from scipy.spatial.transform import Rotation as R

from PyQt6.QtWidgets import QApplication, QMainWindow, QLabel, QVBoxLayout, QWidget
from PyQt6.QtGui import QImage, QPixmap
from PyQt6.QtCore import Qt, QTimer

from spirulae_splat.viewer import Camera, SplatModel


class RenderViewer(QMainWindow):
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
        
        self.c2w = np.array([
            [1, 0, 0, 0],
            [0, 0, 1, 0],
            [0, -1, 0, 0],
            [0, 0, 0, 1]
        ], dtype=np.float32)
        self.update_image()
        
        self.last_mouse_pos = None
        self.active_keys = set()

        self.timer = QTimer(self)
        self.timer.timeout.connect(self.animation_frame)
        self.timer.start(16)  # 60 fps

    def update_image(self):
        image = ssplat_model.render(ssplat_camera, self.c2w)
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

    def mouseMoveEvent(self, event):
        if self.last_mouse_pos is not None:
            dx = (event.pos().x() - self.last_mouse_pos.x()) * 0.01
            dy = (event.pos().y() - self.last_mouse_pos.y()) * 0.01
            matR = self._mat_rotate(-dy, -dx, 0.0)
            w2c = np.linalg.inv(self.c2w)
            matT = self._mat_translate(*w2c[:3, 3])
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
        if ord('W') in self.active_keys:
            mat = self._mat_rotate(rotate_speed, 0, 0)
        elif ord('S') in self.active_keys:
            mat = self._mat_rotate(-rotate_speed, 0, 0)
        elif ord('A') in self.active_keys:
            mat = self._mat_rotate(0, -rotate_speed, 0)
        elif ord('D') in self.active_keys:
            mat = self._mat_rotate(0, rotate_speed, 0)
        elif ord('Q') in self.active_keys:
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

    camera_path = os.path.join(os.path.dirname(__file__), "camera.yaml")
    # camera_path = os.path.join(os.path.dirname(__file__), "camera_t265.yaml")
    ssplat_camera = Camera(camera_path)

    file_path = "/home/harry7557558/nerfstudio/data/sfm_test/mosque/outputs/8/spirulae/2024-12-24_141831/nerfstudio_models/step-000029999.ckpt"
    if len(sys.argv) > 1:
        file_path = sys.argv[1]

    ssplat_model = SplatModel(file_path)
    ssplat_model.bgr = False

    app = QApplication(sys.argv)
    viewer = RenderViewer()
    viewer.show()
    sys.exit(app.exec())
