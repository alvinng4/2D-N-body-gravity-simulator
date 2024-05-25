from pathlib import Path
import sys

path = Path(__file__).parent.parent.parent.absolute() / "gravity_sim"
sys.path.insert(0, str(path))

from camera import Camera


def test_camera_movement():
    camera = Camera()
    camera.moving_right = True
    camera.moving_up = True
    assert camera.pos[0] == 0
    assert camera.pos[1] == 0
    camera.update_movement()
    assert camera.pos[0] == 10
    assert camera.pos[1] == -10
    camera.moving_left = True
    camera.moving_down = True
    camera.update_movement()
    assert camera.pos[0] == 10
    assert camera.pos[1] == -10
    camera.moving_right = False
    camera.moving_up = False
    camera.update_movement()
    camera.update_movement()
    assert camera.pos[0] == -10
    assert camera.pos[1] == 10
