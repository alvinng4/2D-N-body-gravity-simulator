import sys

import pygame
import pygame.gfxdraw
from pygame.sprite import Sprite


class Grav_obj(Sprite):
    def __init__(self, grav_sim, x0, y0, img_path: str = None):
        super().__init__()
        self.screen = grav_sim.screen
        self.settings = grav_sim.settings
        self.color = self.settings.GRAV_OBJ_COLOR

        if img_path:
            try:
                self.image = pygame.image.load(img_path)
                self.rect = self.image.get_rect()
            except FileNotFoundError:
                sys.exit(
                    "Error: Image not found. Check that you are in the correct directory, and the image path provided for Grav_obj is correct."
                )

        self.rect.centerx = x0
        self.rect.centery = y0

    def draw(self):
        """Draw the object at its current location."""
        self.screen.blit(self.image, self.rect)


def solar_system(screen):
    """Initialize the solar system"""
    sun = Grav_obj(screen, "images/sun.png")
    mercury = Grav_obj(screen, "images/mercury.png")
    venus = Grav_obj(screen, "images/venus.png")
    earth = Grav_obj(screen, "images/earth.png")
    mars = Grav_obj(screen, "images/mars.png")
    jupiter = Grav_obj(screen, "images/jupiter.png")
    saturn = Grav_obj(screen, "images/saturn.png")
    uranus = Grav_obj(screen, "images/uranus.png")
    neptune = Grav_obj(screen, "images/neptune.png")
