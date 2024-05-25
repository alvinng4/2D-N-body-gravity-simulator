"""
test file for gravity_sim/settings.py
This file mainly test the setter and getter values. 
Note that 
1. some variables are int 
2. some variables round to ndigits=15 to avoid round off errors
3. Max min rk iteration are coupled.
"""
from pathlib import Path
import random
import sys

path = Path(__file__).parent.parent.parent.absolute() / "gravity_sim"
sys.path.insert(0, str(path))

from settings import Settings


def test_max_range():
    settings = Settings(screen_width=1920, screen_height=1080)
    assert settings.MAX_RANGE > 0
    # MAX_RANGE * self.settings.distance_scale + self.screen_rect.center - self.camera.pos must be smaller than pygame int value
    assert abs(settings.MAX_RANGE) * settings.MAX_DISTANCE_SCALE + 100000 < 2147483647


def test_settings_screen_size():
    settings1 = Settings(screen_width=1920, screen_height=1080)
    assert settings1._screen_width == 1920
    assert settings1._screen_height == 1080
    settings2 = Settings(screen_width=200.5, screen_height=201.5)
    assert settings2._screen_width == 200.5
    assert settings2._screen_height == 201.5
    settings3 = Settings(screen_width=-100, screen_height=-101)
    assert settings3._screen_width == 0
    assert settings3._screen_height == 0
    settings4 = Settings(screen_width=0, screen_height=0)
    assert settings4._screen_width == 0
    assert settings4._screen_height == 0


def test_settings_star_img_scale():
    settings = Settings(screen_width=1920, screen_height=1080)
    assert settings.MAX_STAR_IMG_SCALE > settings.MIN_STAR_IMG_SCALE
    assert settings.star_img_scale == settings.DEFAULT_STAR_IMG_SCALE
    random_value = random.randint(
        settings.MIN_STAR_IMG_SCALE, settings.MAX_STAR_IMG_SCALE
    )
    settings.star_img_scale = random_value
    assert settings.star_img_scale == random_value
    test_value = int((settings.MAX_STAR_IMG_SCALE - settings.MIN_STAR_IMG_SCALE) / 2)
    settings.star_img_scale = settings.MAX_STAR_IMG_SCALE - test_value
    assert settings.star_img_scale == settings.MAX_STAR_IMG_SCALE - test_value
    settings.star_img_scale = settings.MAX_STAR_IMG_SCALE + test_value
    assert settings.star_img_scale == settings.MAX_STAR_IMG_SCALE
    settings.star_img_scale = settings.MIN_STAR_IMG_SCALE + test_value
    assert settings.star_img_scale == settings.MIN_STAR_IMG_SCALE + test_value
    settings.star_img_scale = settings.MIN_STAR_IMG_SCALE - test_value
    assert settings.star_img_scale == settings.MIN_STAR_IMG_SCALE


def test_settings_planet_img_scale():
    settings = Settings(screen_width=1920, screen_height=1080)
    assert settings.MAX_PLANET_IMG_SCALE > settings.MIN_PLANET_IMG_SCALE
    assert settings.planet_img_scale == settings.DEFAULT_PLANET_IMG_SCALE
    random_value = random.randint(
        settings.MIN_PLANET_IMG_SCALE, settings.MAX_PLANET_IMG_SCALE
    )
    settings.planet_img_scale = random_value
    assert settings.planet_img_scale == random_value
    test_value = int(
        (settings.MAX_PLANET_IMG_SCALE - settings.MIN_PLANET_IMG_SCALE) / 2
    )
    settings.planet_img_scale = settings.MAX_PLANET_IMG_SCALE - test_value
    assert settings.planet_img_scale == settings.MAX_PLANET_IMG_SCALE - test_value
    settings.planet_img_scale = settings.MAX_PLANET_IMG_SCALE + test_value
    assert settings.planet_img_scale == settings.MAX_PLANET_IMG_SCALE
    settings.planet_img_scale = settings.MIN_PLANET_IMG_SCALE + test_value
    assert settings.planet_img_scale == settings.MIN_PLANET_IMG_SCALE + test_value
    settings.planet_img_scale = settings.MIN_PLANET_IMG_SCALE - test_value
    assert settings.planet_img_scale == settings.MIN_PLANET_IMG_SCALE


def test_settings_distance_scale():
    settings = Settings(screen_width=1920, screen_height=1080)
    assert settings.MAX_DISTANCE_SCALE > settings.MIN_DISTANCE_SCALE
    assert settings.distance_scale == settings.DEFAULT_DISTANCE_SCALE
    random_value = random.randint(
        settings.MIN_DISTANCE_SCALE, settings.MAX_DISTANCE_SCALE
    )
    settings.distance_scale = random_value
    assert settings.distance_scale == random_value
    test_value = int((settings.MAX_DISTANCE_SCALE - settings.MIN_DISTANCE_SCALE) / 2)
    settings.distance_scale = settings.MAX_DISTANCE_SCALE - test_value
    assert settings.distance_scale == settings.MAX_DISTANCE_SCALE - test_value
    settings.distance_scale = settings.MAX_DISTANCE_SCALE + test_value
    assert settings.distance_scale == settings.MAX_DISTANCE_SCALE
    settings.distance_scale = settings.MIN_DISTANCE_SCALE + test_value
    assert settings.distance_scale == settings.MIN_DISTANCE_SCALE + test_value
    settings.distance_scale = settings.MIN_DISTANCE_SCALE - test_value
    assert settings.distance_scale == settings.MIN_DISTANCE_SCALE


def test_new_star_mass_scale():
    settings = Settings(screen_width=1920, screen_height=1080)
    assert settings.MAX_NEW_STAR_MASS_SCALE > settings.MIN_NEW_STAR_MASS_SCALE
    assert settings.new_star_mass_scale == settings.DEFAULT_NEW_STAR_MASS_SCALE
    random_value = round(
        random.uniform(
            settings.MIN_NEW_STAR_MASS_SCALE, settings.MAX_NEW_STAR_MASS_SCALE
        ),
        ndigits=15,
    )
    settings.new_star_mass_scale = random_value
    assert settings.new_star_mass_scale == random_value
    test_value = (
        settings.MAX_NEW_STAR_MASS_SCALE - settings.MIN_NEW_STAR_MASS_SCALE
    ) / 2.0
    settings.new_star_mass_scale = round(
        settings.MAX_NEW_STAR_MASS_SCALE - test_value, ndigits=15
    )
    assert settings.new_star_mass_scale == round(
        settings.MAX_NEW_STAR_MASS_SCALE - test_value, ndigits=15
    )
    settings.new_star_mass_scale = round(
        settings.MAX_NEW_STAR_MASS_SCALE + test_value, ndigits=15
    )
    assert settings.new_star_mass_scale == round(
        settings.MAX_NEW_STAR_MASS_SCALE, ndigits=15
    )
    settings.new_star_mass_scale = round(
        settings.MIN_NEW_STAR_MASS_SCALE + test_value, ndigits=15
    )
    assert settings.new_star_mass_scale == round(
        settings.MIN_NEW_STAR_MASS_SCALE + test_value, ndigits=15
    )
    settings.new_star_mass_scale = round(
        settings.MIN_NEW_STAR_MASS_SCALE - test_value, ndigits=15
    )
    assert settings.new_star_mass_scale == round(
        settings.MIN_NEW_STAR_MASS_SCALE, ndigits=15
    )


def test_new_star_speed_scale():
    settings = Settings(screen_width=1920, screen_height=1080)
    assert settings.MAX_NEW_STAR_SPEED_SCALE > settings.MIN_NEW_STAR_SPEED_SCALE
    assert settings.new_star_speed_scale == settings.DEFAULT_NEW_STAR_SPEED_SCALE
    random_value = random.randint(
        settings.MIN_NEW_STAR_SPEED_SCALE, settings.MAX_NEW_STAR_SPEED_SCALE
    )
    settings.new_star_speed_scale = random_value
    assert settings.new_star_speed_scale == random_value
    test_value = int(
        (settings.MAX_NEW_STAR_SPEED_SCALE - settings.MIN_NEW_STAR_SPEED_SCALE) / 2
    )
    settings.new_star_speed_scale = settings.MAX_NEW_STAR_SPEED_SCALE - test_value
    assert (
        settings.new_star_speed_scale == settings.MAX_NEW_STAR_SPEED_SCALE - test_value
    )
    settings.new_star_speed_scale = settings.MAX_NEW_STAR_SPEED_SCALE + test_value
    assert settings.new_star_speed_scale == settings.MAX_NEW_STAR_SPEED_SCALE
    settings.new_star_speed_scale = settings.MIN_NEW_STAR_SPEED_SCALE + test_value
    assert (
        settings.new_star_speed_scale == settings.MIN_NEW_STAR_SPEED_SCALE + test_value
    )
    settings.new_star_speed_scale = settings.MIN_NEW_STAR_SPEED_SCALE - test_value
    assert settings.new_star_speed_scale == settings.MIN_NEW_STAR_SPEED_SCALE


def test_settings_dt():
    settings = Settings(screen_width=1920, screen_height=1080)
    assert settings.MAX_DT > settings.MIN_DT
    assert settings.dt == settings.DEFAULT_DT
    random_value = round(random.uniform(settings.MIN_DT, settings.MAX_DT), ndigits=15)
    settings.dt = random_value
    assert settings.dt == random_value
    test_value = (settings.MAX_DT - settings.MIN_DT) / 2.0
    settings.dt = round(settings.MAX_DT - test_value, ndigits=15)
    assert settings.dt == round(settings.MAX_DT - test_value, ndigits=15)
    settings.dt = round(settings.MAX_DT + test_value, ndigits=15)
    assert settings.dt == round(settings.MAX_DT, ndigits=15)
    settings.dt = round(settings.MIN_DT + test_value, ndigits=15)
    assert settings.dt == round(settings.MIN_DT + test_value, ndigits=15)
    settings.dt = round(settings.MIN_DT - test_value, ndigits=15)
    assert settings.dt == round(settings.MIN_DT, ndigits=15)


def test_settings_time_speed():
    settings = Settings(screen_width=1920, screen_height=1080)
    assert settings.MAX_TIME_SPEED > settings.MIN_TIME_SPEED
    assert settings.time_speed == settings.DEFAULT_TIME_SPEED
    random_value = random.randint(settings.MIN_TIME_SPEED, settings.MAX_TIME_SPEED)
    settings.time_speed = random_value
    assert settings.time_speed == random_value
    test_value = int((settings.MAX_TIME_SPEED - settings.MIN_TIME_SPEED) / 2)
    settings.time_speed = settings.MAX_TIME_SPEED - test_value
    assert settings.time_speed == settings.MAX_TIME_SPEED - test_value
    settings.time_speed = settings.MAX_TIME_SPEED + test_value
    assert settings.time_speed == settings.MAX_TIME_SPEED
    settings.time_speed = settings.MIN_TIME_SPEED + test_value
    assert settings.time_speed == settings.MIN_TIME_SPEED + test_value
    settings.time_speed = settings.MIN_TIME_SPEED - test_value
    assert settings.time_speed == settings.MIN_TIME_SPEED


def test_settings_max_iteration():
    settings = Settings(screen_width=1920, screen_height=1080)
    assert settings.MAX_MAX_ITERATION > settings.MIN_MAX_ITERATION
    assert settings.max_iteration == settings.DEFAULT_MAX_ITERATION
    settings.min_iteration = (
        settings.MIN_MAX_ITERATION
    )  # This is neccessary as max and min are coupled
    random_value = random.randint(
        settings.MIN_MAX_ITERATION, settings.MAX_MAX_ITERATION
    )
    settings.max_iteration = random_value
    assert settings.max_iteration == random_value
    test_value = int((settings.MAX_MAX_ITERATION - settings.MIN_MAX_ITERATION) / 2)
    settings.max_iteration = settings.MAX_MAX_ITERATION - test_value
    assert settings.max_iteration == settings.MAX_MAX_ITERATION - test_value
    settings.max_iteration = settings.MAX_MAX_ITERATION + test_value
    assert settings.max_iteration == settings.MAX_MAX_ITERATION
    settings.max_iteration = settings.MIN_MAX_ITERATION + test_value
    assert settings.max_iteration == settings.MIN_MAX_ITERATION + test_value
    settings.max_iteration = settings.MIN_MAX_ITERATION - test_value
    assert settings.max_iteration == settings.MIN_MAX_ITERATION


def test_settings_min_iteration():
    settings = Settings(screen_width=1920, screen_height=1080)
    assert settings.MAX_MIN_ITERATION > settings.MIN_MIN_ITERATION
    assert settings.min_iteration == settings.DEFAULT_MIN_ITERATION
    settings.max_iteration = (
        settings.MAX_MAX_ITERATION
    )  # This is neccessary as max and min are coupled
    random_value = random.randint(
        settings.MIN_MIN_ITERATION, settings.MIN_MIN_ITERATION
    )
    settings.min_iteration = random_value
    assert settings.min_iteration == random_value
    test_value = int((settings.MIN_MIN_ITERATION - settings.MIN_MIN_ITERATION) / 2)
    settings.min_iteration = settings.MIN_MIN_ITERATION - test_value
    assert settings.min_iteration == settings.MIN_MIN_ITERATION - test_value
    settings.min_iteration = settings.MIN_MIN_ITERATION + test_value
    assert settings.min_iteration == settings.MIN_MIN_ITERATION
    settings.min_iteration = settings.MIN_MIN_ITERATION + test_value
    assert settings.min_iteration == settings.MIN_MIN_ITERATION + test_value
    settings.min_iteration = settings.MIN_MIN_ITERATION - test_value
    assert settings.min_iteration == settings.MIN_MIN_ITERATION


def test_settings_tolerance():
    settings = Settings(screen_width=1920, screen_height=1080)
    assert settings.MAX_TOLERANCE > settings.MIN_TOLERANCE
    assert settings.tolerance == settings.DEFAULT_TOLERANCE
    random_value = round(
        random.uniform(settings.MIN_TOLERANCE, settings.MAX_TOLERANCE), ndigits=15
    )
    settings.tolerance = random_value
    assert settings.tolerance == random_value
    test_value = (settings.MAX_TOLERANCE - settings.MIN_TOLERANCE) / 2.0
    settings.tolerance = round(settings.MAX_TOLERANCE - test_value, ndigits=15)
    assert settings.tolerance == round(settings.MAX_TOLERANCE - test_value, ndigits=15)
    settings.tolerance = round(settings.MAX_TOLERANCE + test_value, ndigits=15)
    assert settings.tolerance == round(settings.MAX_TOLERANCE, ndigits=15)
    settings.tolerance = round(settings.MIN_TOLERANCE + test_value, ndigits=15)
    assert settings.tolerance == round(settings.MIN_TOLERANCE + test_value, ndigits=15)
    settings.tolerance = round(settings.MIN_TOLERANCE - test_value, ndigits=15)
    assert settings.tolerance == round(settings.MIN_TOLERANCE, ndigits=15)


def test_settings_expected_time_scale():
    settings = Settings(screen_width=1920, screen_height=1080)
    assert settings.MAX_EXPECTED_TIME_SCALE > settings.MIN_EXPECTED_TIME_SCALE
    assert settings.expected_time_scale == settings.DEFAULT_EXPECTED_TIME_SCALE
    random_value = random.uniform(
        settings.MIN_EXPECTED_TIME_SCALE, settings.MAX_EXPECTED_TIME_SCALE
    )
    settings.expected_time_scale = random_value
    assert settings.expected_time_scale == random_value
    test_value = (
        settings.MAX_EXPECTED_TIME_SCALE - settings.MIN_EXPECTED_TIME_SCALE
    ) / 2.0
    settings.expected_time_scale = settings.MAX_EXPECTED_TIME_SCALE - test_value
    assert settings.expected_time_scale == settings.MAX_EXPECTED_TIME_SCALE - test_value
    settings.expected_time_scale = settings.MAX_EXPECTED_TIME_SCALE + test_value
    assert settings.expected_time_scale == settings.MAX_EXPECTED_TIME_SCALE
    settings.expected_time_scale = settings.MIN_EXPECTED_TIME_SCALE + test_value
    assert settings.expected_time_scale == settings.MIN_EXPECTED_TIME_SCALE + test_value
    settings.expected_time_scale = settings.MIN_EXPECTED_TIME_SCALE - test_value
    assert settings.expected_time_scale == settings.MIN_EXPECTED_TIME_SCALE
