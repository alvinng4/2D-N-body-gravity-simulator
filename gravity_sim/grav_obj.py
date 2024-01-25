import sys
import os

import pygame
from pygame.sprite import Sprite

from settings import Settings


class Grav_obj(Sprite):
    # Gravitational constant (AU^3 / d^2 / M_sun):
    G = 0.00029591220828559

    # Solar radius in AU
    SOLAR_RADIUS = 0.004650467261

    def __init__(
        self,
        grav_sim,
        params: dict,
        img_path: str = None,
        name: str = None,
    ):
        super().__init__()
        self.screen = grav_sim.screen
        self.screen_rect = self.screen.get_rect()
        self.camera = grav_sim.camera
        self.settings = grav_sim.settings
        self.params = params
        self.diameter = 2 * self.params["R"]
        if name == "Sun":
            self.img_diameter = self.diameter * self.settings.star_img_scale
        else:
            self.img_diameter = self.diameter * self.settings.planet_img_scale

        if img_path:
            try:
                load_image = pygame.image.load(img_path).convert_alpha()
                self.image = pygame.transform.scale(
                    load_image, (self.img_diameter, self.img_diameter)
                )
                self.rect = self.image.get_rect()
            except FileNotFoundError:
                sys.exit(
                    "Error: Image not found. Make sure the image path provided for Grav_obj is correct."
                )

    def update(self):
        """Update the apparent position of all grav_objs with camera"""
        self.rect.center = (
            self.params["r1"] * self.settings.distance_scale
            + self.screen_rect.centerx
            - self.camera.pos[0],
            -self.params["r2"] * self.settings.distance_scale
            + self.screen_rect.centery
            - self.camera.pos[1],
        )

    @staticmethod
    def create_star(grav_sim, mouse_pos, camera_pos, drag_mouse_pos, drag_camera_pos):
        main_dir_path = os.path.dirname(__file__)
        path_sun = os.path.join(main_dir_path, "assets/images/sun.png")
        m = (
            1
            * 0.5
            * grav_sim.stats.holding_rclick_time
            * grav_sim.settings.new_star_mass_scale
        )
        R = Grav_obj.SOLAR_RADIUS * (m ** (1.0 / 3.0))
        grav_obj = Grav_obj(
            grav_sim,
            {
                "r1": (
                    mouse_pos[0] - grav_sim.screen.get_rect().centerx + camera_pos[0]
                )
                / grav_sim.settings.distance_scale,
                "r2": -(
                    mouse_pos[1] - grav_sim.screen.get_rect().centery + camera_pos[1]
                )
                / grav_sim.settings.distance_scale,
                "r3": 0.0,
                "v1": -(
                    (drag_mouse_pos[0] - mouse_pos[0])
                    + (drag_camera_pos[0] - camera_pos[0])
                )
                * Settings.DEFAULT_NEW_STAR_VELOCITY_SCALE
                / (grav_sim.settings.distance_scale / Settings.DEFAULT_DISTANCE_SCALE),
                "v2": (
                    (drag_mouse_pos[1] - mouse_pos[1])
                    + (drag_camera_pos[1] - camera_pos[1])
                )
                * Settings.DEFAULT_NEW_STAR_VELOCITY_SCALE
                / (grav_sim.settings.distance_scale / Settings.DEFAULT_DISTANCE_SCALE),
                "v3": 0.0,
                "m": m,
                "R": R,
            },
            path_sun,
            name="Sun",
        )
        grav_sim.grav_objs.add(grav_obj)
        grav_sim.simulator.is_initialize = True
        grav_sim.simulator.is_initialize_integrator = grav_sim.simulator.current_integrator

    @staticmethod
    def create_solor_system(grav_sim):
        """
        Create the solar system
        Data dated on A.D. 2024-Jan-01 00:00:00.0000 TDB
        Computational data generated by NASA JPL Horizons System https://ssd.jpl.nasa.gov/horizons/
        """
        main_dir_path = os.path.dirname(__file__)
        path_sun = os.path.join(main_dir_path, "assets/images/sun.png")
        path_mercury = os.path.join(main_dir_path, "assets/images/mercury.png")
        path_venus = os.path.join(main_dir_path, "assets/images/venus.png")
        path_earth = os.path.join(main_dir_path, "assets/images/earth.png")
        path_mars = os.path.join(main_dir_path, "assets/images/mars.png")
        path_jupiter = os.path.join(main_dir_path, "assets/images/jupiter.png")
        path_saturn = os.path.join(main_dir_path, "assets/images/saturn.png")
        path_uranus = os.path.join(main_dir_path, "assets/images/uranus.png")
        path_neptune = os.path.join(main_dir_path, "assets/images/neptune.png")
        # r1 - r3: Positions (AU), v1 - v3: Velocities (AU/d), m: Mass (Solar masses)
        sun = Grav_obj(
            grav_sim,
            {
                "r1": -0.007967955691534,
                "r2": -0.002906227441573,
                "r3": 0.000210305430155,
                "v1": 0.000004875094764,
                "v2": -0.000007057133214,
                "v3": -0.000000045734537,
                "m": 1.0,
                "R": Grav_obj.SOLAR_RADIUS,
            },
            path_sun,
            name="Sun",
        )
        mercury = Grav_obj(
            grav_sim,
            {
                "r1": -0.282598326953863,
                "r2": 0.197455979595808,
                "r3": 0.041774335580637,
                "v1": -0.022321659001897,
                "v2": -0.021572071031763,
                "v3": 0.000285519341050,
                "m": 1.66051140935277e-07,
                "R": 1.63083872e-05,
            },
            path_mercury,
        )
        venus = Grav_obj(
            grav_sim,
            {
                "r1": -0.723210370166638,
                "r2": -0.079483020263124,
                "r3": 0.040428714281743,
                "v1": 0.002034068201002,
                "v2": -0.020208286265930,
                "v3": -0.000394563984386,
                "m": 2.44827371182131e-06,
                "R": 4.04537843e-05,
            },
            path_venus,
        )
        earth = Grav_obj(
            grav_sim,
            {
                "r1": -0.173819201725705,
                "r2": 0.966324555023514,
                "r3": 0.000155390185490,
                "v1": -0.017230012325382,
                "v2": -0.002967721342619,
                "v3": 0.000000638212538,
                "m": 3.00329789031573e-06,
                "R": 4.25875046e-05,
            },
            path_earth,
        )
        mars = Grav_obj(
            grav_sim,
            {
                "r1": -0.301326239258265,
                "r2": -1.454029331393295,
                "r3": -0.023005314339914,
                "v1": 0.014248322593453,
                "v2": -0.001579236181581,
                "v3": -0.000382372279616,
                "m": 3.22773848604808e-07,
                "R": 2.26574081e-05,
            },
            path_mars,
        )
        jupiter = Grav_obj(
            grav_sim,
            {
                "r1": 3.485202469657675,
                "r2": 3.552136904413157,
                "r3": -0.092710354427984,
                "v1": -0.005470970658852,
                "v2": 0.005642487338479,
                "v3": 0.000098961906021,
                "m": 0.000954532562518104,
                "R": 4.6732617e-04,
            },
            path_jupiter,
        )
        saturn = Grav_obj(
            grav_sim,
            {
                "r1": 8.988104223143450,
                "r2": -3.719064854634689,
                "r3": -0.293193777732359,
                "v1": 0.001822013845554,
                "v2": 0.005143470425888,
                "v3": -0.000161723590489,
                "m": 0.00028579654259599,
                "R": 3.89256877e-04 * (157 / 57),  # Img scale for Saturn's ring
            },
            path_saturn,
        )
        uranus = Grav_obj(
            grav_sim,
            {
                "r1": 12.263024178975050,
                "r2": 15.297387924805450,
                "r3": -0.102054902688356,
                "v1": -0.003097615358317,
                "v2": 0.002276781932346,
                "v3": 0.000048604332222,
                "m": 4.3655207025844e-05,
                "R": 1.69534499e-04,
            },
            path_uranus,
        )
        neptune = Grav_obj(
            grav_sim,
            {
                "r1": 29.835014609847410,
                "r2": -1.793812957956852,
                "r3": -0.650640113225459,
                "v1": 0.000167653661182,
                "v2": 0.003152098732862,
                "v3": -0.000068775010957,
                "m": 5.1499991953912e-05,
                "R": 1.64587904e-04,
            },
            path_neptune,
        )
        grav_sim.grav_objs.add(sun)
        grav_sim.grav_objs.add(mercury)
        grav_sim.grav_objs.add(venus)
        grav_sim.grav_objs.add(earth)
        grav_sim.grav_objs.add(mars)
        grav_sim.grav_objs.add(jupiter)
        grav_sim.grav_objs.add(saturn)
        grav_sim.grav_objs.add(uranus)
        grav_sim.grav_objs.add(neptune)

    @staticmethod
    def create_figure_8(grav_sim):
        """
        Create a figure-8 orbit
        Data from the book Moving Planets Around: An Introduction to
        N-Body Simulations Applied to Exoplanetary Systems, Ch.7, Page 109
        As the data given use G = 1, the mass is converted by m / G, since a = GM/r^2.
        """
        # Currently use sun as object. May or may not change later.
        main_dir_path = os.path.dirname(__file__)
        path_sun = os.path.join(main_dir_path, "assets/images/sun.png")
        object_1 = Grav_obj(
            grav_sim,
            {
                "r1": 0.970043,
                "r2": -0.24308753,
                "r3": 0.0,
                "v1": 0.466203685,
                "v2": 0.43236573,
                "v3": 0.0,
                "m": 1.0 / Grav_obj.G,
                "R": Grav_obj.SOLAR_RADIUS,  # The radius is arbitrary. Here we give it the solar radii as it have 1 solar mass
            },
            path_sun,
            name="Sun",
        )
        object_2 = Grav_obj(
            grav_sim,
            {
                "r1": -0.970043,
                "r2": 0.24308753,
                "r3": 0.0,
                "v1": 0.466203685,
                "v2": 0.43236573,
                "v3": 0.0,
                "m": 1.0 / Grav_obj.G,
                "R": Grav_obj.SOLAR_RADIUS,  # The radius is arbitrary. Here we give it the solar radii as it have 1 solar mass
            },
            path_sun,
            name="Sun",
        )
        object_3 = Grav_obj(
            grav_sim,
            {
                "r1": 0.0,
                "r2": 0.0,
                "r3": 0.0,
                "v1": -0.93240737,
                "v2": -0.86473146,
                "v3": 0.0,
                "m": 1.0 / Grav_obj.G,
                "R": Grav_obj.SOLAR_RADIUS,  # The radius is arbitrary. Here we give it the solar radii as it have 1 solar mass
            },
            path_sun,
            name="Sun",
        )
        grav_sim.grav_objs.add(object_1)
        grav_sim.grav_objs.add(object_2)
        grav_sim.grav_objs.add(object_3)

    @staticmethod
    def create_pyth_3_body(grav_sim):
        """
        Create a Pythagorean three-body orbit
        Data from the book Moving Planets Around: An Introduction to
        N-Body Simulations Applied to Exoplanetary Systems, Ch.7, Page 109
        As the data given use G = 1, the mass is converted by m / G, since a = GM / r^2.
        """
        # Currently use sun as object. May or may not change later.
        main_dir_path = os.path.dirname(__file__)
        path_sun = os.path.join(main_dir_path, "assets/images/sun.png")
        object_1 = Grav_obj(
            grav_sim,
            {
                "r1": 1.0,
                "r2": 3.0,
                "r3": 0.0,
                "v1": 0.0,
                "v2": 0.0,
                "v3": 0.0,
                "m": 3.0 / Grav_obj.G,
                "R": Grav_obj.SOLAR_RADIUS,  # The radius is arbitrary. Here we give it the solar radii as it have 1 solar mass
            },
            path_sun,
            name="Sun",
        )
        object_2 = Grav_obj(
            grav_sim,
            {
                "r1": -2.0,
                "r2": -1.0,
                "r3": 0.0,
                "v1": 0.0,
                "v2": 0.0,
                "v3": 0.0,
                "m": 4.0 / Grav_obj.G,
                "R": Grav_obj.SOLAR_RADIUS,  # The radius is arbitrary. Here we give it the solar radii as it have 1 solar mass
            },
            path_sun,
            name="Sun",
        )
        object_3 = Grav_obj(
            grav_sim,
            {
                "r1": 1.0,
                "r2": -1.0,
                "r3": 0.0,
                "v1": 0.0,
                "v2": 0.0,
                "v3": 0.0,
                "m": 5.0 / Grav_obj.G,
                "R": Grav_obj.SOLAR_RADIUS,  # The radius is arbitrary. Here we give it the solar radii as it have 1 solar mass
            },
            path_sun,
            name="Sun",
        )
        grav_sim.grav_objs.add(object_1)
        grav_sim.grav_objs.add(object_2)
        grav_sim.grav_objs.add(object_3)
