import sys
import numpy as np
from scipy.spatial.transform import Rotation as R

from PyQt6.QtWidgets import QApplication, QMainWindow, QLabel, QVBoxLayout, QWidget
from PyQt6.QtGui import QImage, QPixmap
from PyQt6.QtCore import Qt

from model import SplatModel
from camera import Camera

ssplat_camera = Camera("viewer/camera.yaml")
ssplat_model = SplatModel("/home/harry7557558/nerfstudio/data/sfm_test/mosque/outputs/8/spirulae/2024-12-24_141831/nerfstudio_models/step-000029999.ckpt")

def render(c2w: np.ndarray):
    ssplat_model.bgr = False
    return ssplat_model.render(ssplat_camera, c2w)


class RenderViewer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Render Viewer")
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

    def update_image(self):
        image = render(self.c2w)
        h, w, c = image.shape
        qimage = QImage(image.data, w, h, 3 * w, QImage.Format.Format_RGB888)
        self.image_label.setPixmap(QPixmap.fromImage(qimage))

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.last_mouse_pos = event.pos()

    def mouseMoveEvent(self, event):
        if self.last_mouse_pos is not None:
            dx = (event.pos().x() - self.last_mouse_pos.x()) * 0.01
            dy = (event.pos().y() - self.last_mouse_pos.y()) * 0.01
            rot = R.from_rotvec(-dy * self.c2w.T[0,:3]) * R.from_rotvec(-dx * self.c2w.T[1,:3])
            self.c2w[:3] = rot.as_matrix() @ self.c2w[:3]
            self.last_mouse_pos = event.pos()
            self.update_image()

    def wheelEvent(self, event):
        zoom_factor = event.angleDelta().y() / 2000
        self.c2w[:3, 3] += zoom_factor * self.c2w[:3, 2]
        self.update_image()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    viewer = RenderViewer()
    viewer.show()
    sys.exit(app.exec())
