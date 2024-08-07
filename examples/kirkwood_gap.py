"""
Demonstration on using the gravity simulator to simulate the Kirkwood gap.
You will need to install the `Pillow` library for this script.

Note: N = 1000 is enough to observe some gaps, but it may not be very clear.
Since the simulation is O(N), seting N = 50000 would take 24 hours to a few days to 
finish, and reducing N to 25000 will reduce the runtime by half.

Warning: This script will take a lot of storage space on your computer (probably a few GBs).
         When combining the individual frames, the pillow library will take a lot of memories.
         I have reduced number of asteroids, dpi and frame size to reduce the runtime.
         You may increase them if you want better quality.

         The script will ask user's permission to erase the data after the video is generated.

TODO: Calculations for the 2D scatter plot is not vectorized and is extremely slow. 
"""

import csv
from pathlib import Path
import sys
sys.path.append(str(Path(__file__).parent.parent))

import numpy as np
import PIL
import matplotlib.pyplot as plt

from gravity_sim import GravitySimulator
from gravity_sim.common import get_bool
from gravity_sim.common import Progress_bar

N = 10000
FPS = 30
DPI = 150

G = 0.00029591220828411956
M = 1.0


def main():
    # ---------- Initialization ---------- #

    grav_sim = GravitySimulator()
    system = grav_sim.create_system()
    grav_sim.set_current_system(system)

    system.load("solar_system")
    system.remove(name="Mercury")
    system.remove(name="Venus")
    system.remove(name="Earth")
    system.remove(name="Neptune")
    system.remove(name="Uranus")
    colors = [
        "orange",
        "red",
        "darkgoldenrod",
        "gold",
    ]
    labels = system.objects_names.copy()
    marker_sizes = [6.0, 1.5, 4.0, 3.5]

    massive_objects_count = system.objects_count
    inner_objects_count = 2
    outer_objects_count = 2

    rng = np.random.default_rng()
    a = rng.uniform(2.0, 3.35, size=N)
    ecc = np.abs(rng.normal(loc=0.0, scale=0.12, size=N))
    inc = np.abs(rng.normal(loc=0.0, scale=0.3, size=N))
    argument_of_periapsis = rng.uniform(0, 2 * np.pi, size=N)
    long_asc_node = rng.uniform(0, 2 * np.pi, size=N)
    true_anomaly = rng.uniform(0, 2 * np.pi, size=N)

    for i in range(N):
        system.add_keplerian(
            semi_major_axis=a[i],
            eccentricity=ecc[i],
            inclination=inc[i],
            argument_of_periapsis=argument_of_periapsis[i],
            longitude_of_ascending_node=long_asc_node[i],
            true_anomaly=true_anomaly[i],
            m=0.0,
            primary_object_name="Sun",
        )
    system.sort_by_distance(primary_object_name="Sun")
    system.center_of_mass_correction()
    system.name = f"kirkwood_gap_N{N}"

    # ---------- Simulation ---------- #
    file_path = Path(__file__).parent.parent / "gravity_sim" / "results"
    file_path.mkdir(parents=True, exist_ok=True)
    data_path = file_path / "kirkwood_gap_sim.csv"

    print("Simulating kirkwood gap...")
    # Store about 2000 points in total
    tf = grav_sim.years_to_days(5000000)
    dt = 180.0
    store_every_n = int((grav_sim.years_to_days(5000000) // dt) // 500)
    grav_sim.launch_simulation(
        "whfast",
        tf=tf,
        dt=dt,
        store_every_n=store_every_n,
        acceleration_method="massless",
        storing_method="flush",
        flush_results_path=str(data_path),
        no_print=True,
        kepler_tol=1e-12,
        kepler_max_iter=500,
        kepler_auto_remove=True,
        kepler_auto_remove_tol=1e-8,
    )

    # ---------- Data Analysis and drawing frames ---------- #

    # In the API, we use PillowWriter to generate animations.
    # However, for some reason, the PillowWriter run out of memory
    # in this case. Therefore, we save each frames as images and
    # combine them as gif instead.
    save_count_semi_major_axes = 0
    save_count_visualization = 0
    new_field_lim = sys.maxsize
    while True:
        try:
            csv.field_size_limit(new_field_lim)
            break
        except OverflowError:
            new_field_lim = new_field_lim // 10

    with open(data_path, "r") as file:
        reader = csv.reader(file)

        # Read data size
        data_size = None
        for row in reader:
            if row[0].startswith("#"):
                if row[0].startswith("# Data size: "):
                    data_size = int(row[0].replace("# Data size: ", ""))
                else:
                    pass
            else:
                break

        file.seek(0)

        progress_bar = Progress_bar()

        print()
        print("Drawing frames for semi-major-axes plots...")
        fig1, ax1 = plt.subplots()
        ax1_xlim_min = 1.8
        ax1_xlim_max = 3.5
        ax1_ylim_min = 0.0
        ax1_ylim_max = 1.0
        with progress_bar:
            if data_size is not None:
                task = progress_bar.add_task("", total=data_size)

            for row in reader:
                if len(row) == 0 or row[0].startswith("#"):
                    continue

                year = grav_sim.days_to_years(float(row[0])) / 1e6
                row = row[3:]
                row = list(map(float, row))

                # fmt: off

                # Due to removal of asteroids by Kepler auto clear,
                # we need to count the number of asteroids every row
                asteroids_count = int(len(row) // 6) - massive_objects_count
                sun_x = np.array(row[:3])
                sun_v = np.array(
                    row[(inner_objects_count + asteroids_count + outer_objects_count) * 3 
                        : ( inner_objects_count + asteroids_count + outer_objects_count) * 3 + 3
                    ]
                )
                asteroids_x = (
                    np.array(
                        row[(inner_objects_count * 3) 
                            : (inner_objects_count + asteroids_count) * 3
                        ]
                    ).reshape(asteroids_count, 3)
                )
                asteroids_v = (
                    np.array(
                        row[inner_objects_count * 6 + (asteroids_count + outer_objects_count) * 3
                            : (inner_objects_count + asteroids_count) * 6 + outer_objects_count * 3
                        ]
                    ).reshape(asteroids_count, 3)
                )

                eccentricity = calculate_eccentricity(asteroids_x - sun_x, asteroids_v - sun_v, 0.0, G, M)
                semi_major_axes = calculate_semi_major_axis(asteroids_x - sun_x, asteroids_v - sun_v, 0.0, G, M)
                # fmt: on

                # Scatter plot
                ax1.scatter(
                    semi_major_axes,
                    eccentricity,
                    s=2,
                    marker=".",
                    color="black",
                    alpha=0.5,
                )

                # Annotate time
                ax1.annotate(
                    f"{year:.2f} Myr",
                    xy=(0.95, 0.95),
                    xycoords="axes fraction",
                    fontsize=12,
                    ha="right",
                    va="top",
                )

                # Draw dashed lines indicating the gaps
                # fmt: off
                ax1.axvline(x=2.065, color="black", linestyle="--", linewidth=0.5, alpha=0.2, zorder=0)
                ax1.axvline(x=2.502, color="black", linestyle="--", linewidth=0.5, alpha=0.2, zorder=0)
                ax1.axvline(x=2.825, color="black", linestyle="--", linewidth=0.5, alpha=0.2, zorder=0)
                ax1.axvline(x=2.958, color="black", linestyle="--", linewidth=0.5, alpha=0.2, zorder=0)
                ax1.axvline(x=3.279, color="black", linestyle="--", linewidth=0.5, alpha=0.2, zorder=0)

                # Annotate resonance ratios
                ax1.text(2.065, ax1_ylim_max, '4:1', color='black', ha='center', va='bottom', fontsize=10)
                ax1.text(2.502, ax1_ylim_max, '3:1', color='black', ha='center', va='bottom', fontsize=10)
                ax1.text(2.825, ax1_ylim_max, '5:2', color='black', ha='center', va='bottom', fontsize=10)
                ax1.text(2.958, ax1_ylim_max, '7:3', color='black', ha='center', va='bottom', fontsize=10)
                ax1.text(3.279, ax1_ylim_max, '2:1', color='black', ha='center', va='bottom', fontsize=10)

                # Set labels
                ax1.set_xlabel("Semi-major axis (AU)")
                ax1.set_ylabel("Eccentricity")

                # Set axes
                ax1.set_xlim(ax1_xlim_min, ax1_xlim_max)
                ax1.set_ylim(ax1_ylim_min, ax1_ylim_max)

                fig1.tight_layout()

                # Capture the frame
                plt.savefig(file_path / f"semi_major_axes_frames_{save_count_semi_major_axes:04d}.png", dpi=DPI)
                save_count_semi_major_axes += 1

                # Clear the plot to prepare for the next frame
                ax1.clear()

                if data_size is not None:
                    progress_bar.update(task, completed=save_count_semi_major_axes)
                # fmt: on

        progress_bar.stop_task(task)
        plt.close("all")
        file.seek(0)

        print()
        print("Drawing frames for visualization plots...")
        fig2, ax2 = plt.subplots(figsize=(4.8, 4.8))
        ax2.set_facecolor("black")
        ax2_xlim_min = -3.5
        ax2_xlim_max = 3.5
        ax2_ylim_min = -3.5
        ax2_ylim_max = 3.5
        with progress_bar:
            if data_size is not None:
                task = progress_bar.add_task("", total=data_size)

            for row in reader:
                if len(row) == 0 or row[0].startswith("#"):
                    continue

                year = grav_sim.days_to_years(float(row[0])) / 1e6
                row = row[3:]
                row = list(map(float, row))

                # fmt: off

                # Due to removal of asteroids by Kepler auto clear,
                # we need to count the number of asteroids every row
                asteroids_count = int(len(row) // 6) - massive_objects_count

                sun_x = np.array(row[:3])
                sun_v = np.array(
                    row[(inner_objects_count + asteroids_count + outer_objects_count) * 3 
                        : ( inner_objects_count + asteroids_count + outer_objects_count) * 3 + 3
                    ]
                )
                asteroids_x = (
                    np.array(
                        row[(inner_objects_count * 3) 
                            : (inner_objects_count + asteroids_count) * 3
                        ]
                    ).reshape(asteroids_count, 3)
                )
                asteroids_v = (
                    np.array(
                        row[inner_objects_count * 6 + (asteroids_count + outer_objects_count) * 3
                            : (inner_objects_count + asteroids_count) * 6 + outer_objects_count * 3
                        ]
                    ).reshape(asteroids_count, 3)
                )

                x = np.zeros((asteroids_count, 3))
                for i in range(asteroids_count):
                    semi_major_axis, _, true_anomaly, _, arg_per, long_asc_nodes = cartesian_to_orbital_elements(
                        0.0, M, asteroids_x[i] - sun_x, asteroids_v[i] - sun_v, G=G
                    )
                    x[i] = keplerian_to_cartesian(
                        semi_major_axis=semi_major_axis,
                        eccentricity=0.0,
                        inclination=0.0,
                        argument_of_periapsis=arg_per,
                        longitude_of_ascending_node=long_asc_nodes,
                        true_anomaly=true_anomaly,
                    )

                # Plotting the sun
                ax2.plot(
                    sun_x[0],
                    sun_x[1],
                    "o",
                    label=labels[0],
                    color=colors[0],
                    markersize=marker_sizes[0],
                )

                # Plotting the asteroids
                ax2.scatter(
                    x[:, 0],
                    x[:, 1],
                    color="white",
                    marker=".",
                    s=0.1,
                    alpha=0.4,
                )

                # Annotate time
                ax2.annotate(
                    f"{year:.2f} Myr",
                    xy=(0.95, 0.95),
                    xycoords="axes fraction",
                    fontsize=12,
                    ha="right",
                    va="top",
                    color="white",
                )

                # Set labels
                ax2.set_xlabel("$x$ (AU)")
                ax2.set_ylabel("$y$ (AU)")

                # Set axes
                ax2.set_xlim(ax2_xlim_min, ax2_xlim_max)
                ax2.set_ylim(ax2_ylim_min, ax2_ylim_max)

                # Set equal aspect ratio to prevent distortion
                ax2.set_aspect("equal")

                ax2.set_title("2D scatter plot with correction\n(eccentricity = inclination = 0)")

                fig2.tight_layout()

                # Capture the frame
                plt.savefig(
                    file_path
                    / f"visualization_frames_{save_count_visualization:04d}.png",
                    dpi=DPI,
                )
                save_count_visualization += 1

                # Clear the plot to prepare for the next frame
                ax2.clear()

                if data_size is not None:
                    progress_bar.update(task, completed=save_count_visualization)

        plt.close("all")

    progress_bar.stop_task(task)
    print()
    print("Combining frames to gif...")
    semi_major_axes_frames = []
    for i in range(save_count_semi_major_axes):
        semi_major_axes_frames.append(
            PIL.Image.open(file_path / f"semi_major_axes_frames_{i:04d}.png")
        )

    semi_major_axes_frames[0].save(
        file_path / "Kirkwood_gap_semi_major_axes.gif",
        save_all=True,
        append_images=semi_major_axes_frames[1:],
        loop=0,
        duration=(1000 // FPS),
    )

    visualization_frames = []
    for i in range(save_count_visualization):
        visualization_frames.append(
            PIL.Image.open(file_path / f"visualization_frames_{i:04d}.png")
        )

    visualization_frames[0].save(
        file_path / "Kirkwood_gap_visualization.gif",
        save_all=True,
        append_images=visualization_frames[1:],
        loop=0,
        duration=(1000 // FPS),
    )

    print(
        f"Output completed! Please check {file_path / 'Kirkwood_gap_semi_major_axes.gif'} \nand {file_path / 'Kirkwood_gap_visualization.gif'}"
    )
    print()

    print("Deleting files...")
    for i in range(save_count_semi_major_axes):
        (file_path / f"semi_major_axes_frames_{i:04d}.png").unlink()

    for i in range(save_count_visualization):
        (file_path / f"visualization_frames_{i:04d}.png").unlink()

    if get_bool(f"Delete data file? Path: {data_path}"):
        data_path.unlink()

    print("Done! Exiting the program...")


def calculate_eccentricity(x, v, m, G, M):
    ecc_vec = (
        np.cross(v, np.cross(x, v)) / (G * (m + M))
        - x / np.linalg.norm(x, axis=1)[:, np.newaxis]
    )
    ecc = np.linalg.norm(ecc_vec, axis=1)

    return ecc


def calculate_semi_major_axis(x, v, m, G, M):
    E_sp = 0.5 * np.linalg.norm(v, axis=1) ** 2 - G * (m + M) / np.linalg.norm(
        x, axis=1
    )
    a = -0.5 * G * (m + M) / E_sp

    return a

def cartesian_to_orbital_elements(
    mp, ms, position, velocity, G,
):
    """

    Function that computes orbital elements from cartesian coordinates.
    Return values are: mass1, mass2, semimajor axis, eccentricity,
    true anomaly, inclination, longitude of the ascending nodes and the
    argument of pericenter. All angles are returned in radians.
    In case of a perfectly circular orbit the true anomaly
    and argument of pericenter cannot be determined; in this case, the
    return values are 0.0 for both angles.

    Reference: copied from ABIE with slight modification
    https://github.com/MovingPlanetsAround/ABIE/blob/master/abie/tools.py
    """

    total_mass = mp + ms

    ### specific energy ###
    v_sq = np.dot(velocity, velocity)
    r_sq = np.dot(position, position)
    r = np.sqrt(r_sq)

    specific_energy = (1.0 / 2.0) * v_sq - G * total_mass / r
    # if specific_energy >= 0.0:
    #     print 'Not a bound orbit!'

    semimajor_axis = -G * total_mass / (2.0 * specific_energy)

    ### specific angular momentum ###
    specific_angular_momentum = np.cross(position, velocity)
    specific_angular_momentum_norm = np.sqrt(
        np.dot(specific_angular_momentum, specific_angular_momentum)
    )
    specific_angular_momentum_unit = (
        specific_angular_momentum / specific_angular_momentum_norm
    )

    maximum_specific_angular_momentum_norm = (
        G * total_mass / (np.sqrt(-2.0 * specific_energy))
    )
    ell = (
        specific_angular_momentum_norm / maximum_specific_angular_momentum_norm
    )  ### specific AM in units of maximum AM

    ### for e = 0 or e nearly 0, ell can be slightly larger than unity due to numerical reasons ###
    ell_epsilon = 1e-15

    completely_or_nearly_circular = False

    if ell > 1.0:
        if (
            1.0 < ell <= ell + ell_epsilon
        ):  ### still unity within numerical precision
            print(
                "orbit is completely or nearly circular; in this case the LRL vector cannot be used to reliably obtain the argument of pericenter and true anomaly; the output values of the latter will be set to zero; output e will be e = 0"
            )
            ell = 1.0
            completely_or_nearly_circular = True
        else:  ### larger than unity within numerical precision
            raise Exception(
                "angular momentum larger than maximum angular momentum for bound orbit"
            )

    eccentricity = np.sqrt(1.0 - ell ** 2)

    ### Orbital inclination ###
    z_vector = np.array([0.0, 0.0, 1.0])
    inclination = np.arccos(np.dot(z_vector, specific_angular_momentum_unit))

    ### Longitude of ascending nodes, with reference direction along x-axis ###
    ascending_node_vector = np.cross(z_vector, specific_angular_momentum)
    ascending_node_vector_norm = np.sqrt(
        np.dot(ascending_node_vector, ascending_node_vector)
    )
    if ascending_node_vector_norm == 0:
        ascending_node_vector_unit = np.array([1.0, 0.0, 0.0])
    else:
        ascending_node_vector_unit = (
            ascending_node_vector / ascending_node_vector_norm
        )

    long_asc_nodes = np.arctan2(
        ascending_node_vector_unit[1], ascending_node_vector_unit[0]
    )

    ### Argument of periapsis and true anomaly, using eccentricity a.k.a. Laplace-Runge-Lenz (LRL) vector ###
    mu = G * total_mass
    position_unit = position / r
    e_vector = (1.0 / mu) * np.cross(
        velocity, specific_angular_momentum
    ) - position_unit  ### Laplace-Runge-Lenz vector

    if (
        completely_or_nearly_circular == True
    ):  ### orbit is completely or nearly circular; in this case the LRL vector cannot be used to reliably obtain the argument of pericenter and true anomaly; the output values of the latter will be set to zero; output e will be e = 0
        arg_per = 0.0
        true_anomaly = 0.0
    else:
        e_vector_norm = np.sqrt(np.dot(e_vector, e_vector))
        e_vector_unit = e_vector / e_vector_norm

    e_vector_unit_cross_AM_unit = np.cross(
        e_vector_unit, specific_angular_momentum_unit
    )
    sin_arg_per = np.dot(ascending_node_vector_unit, e_vector_unit_cross_AM_unit)
    cos_arg_per = np.dot(e_vector_unit, ascending_node_vector_unit)
    arg_per = np.arctan2(sin_arg_per, cos_arg_per)

    sin_true_anomaly = np.dot(position_unit, -1.0 * e_vector_unit_cross_AM_unit)
    cos_true_anomaly = np.dot(position_unit, e_vector_unit)
    true_anomaly = np.arctan2(sin_true_anomaly, cos_true_anomaly)

    return (
        semimajor_axis,
        eccentricity,
        true_anomaly,
        inclination,
        arg_per,
        long_asc_nodes,
    )

def keplerian_to_cartesian(
    semi_major_axis: float,
    eccentricity: float,
    inclination: float,
    argument_of_periapsis: float,
    longitude_of_ascending_node: float,
    true_anomaly: float,
):
    """
    Convert keplerian elements to cartesian coordinates,
    modified to calculate position only.

    Reference
    ---------
    Moving Planets Around: An Introduction to N-Body
    Simulations Applied to Exoplanetary Systems, Chapter 2
    """

    cos_inc = np.cos(inclination)
    sin_inc = np.sin(inclination)

    cos_arg_periapsis = np.cos(argument_of_periapsis)
    sin_arg_periapsis = np.sin(argument_of_periapsis)

    cos_long_asc_node = np.cos(longitude_of_ascending_node)
    sin_long_asc_node = np.sin(longitude_of_ascending_node)

    cos_true_anomaly = np.cos(true_anomaly)
    sin_true_anomaly = np.sin(true_anomaly)

    # ecc_unit_vec is the unit vector pointing towards periapsis
    ecc_unit_vec = np.zeros(3)
    ecc_unit_vec[0] = (
        cos_long_asc_node * cos_arg_periapsis
        - sin_long_asc_node * sin_arg_periapsis * cos_inc
    )
    ecc_unit_vec[1] = (
        sin_long_asc_node * cos_arg_periapsis
        + cos_long_asc_node * sin_arg_periapsis * cos_inc
    )
    ecc_unit_vec[2] = sin_arg_periapsis * sin_inc

    # q_unit_vec is the unit vector that is perpendicular to ecc_unit_vec and orbital angular momentum vector
    q_unit_vec = np.zeros(3)
    q_unit_vec[0] = (
        -cos_long_asc_node * sin_arg_periapsis
        - sin_long_asc_node * cos_arg_periapsis * cos_inc
    )
    q_unit_vec[1] = (
        -sin_long_asc_node * sin_arg_periapsis
        + cos_long_asc_node * cos_arg_periapsis * cos_inc
    )
    q_unit_vec[2] = cos_arg_periapsis * sin_inc

    # Calculate the position vector
    x = (
        semi_major_axis
        * (1.0 - eccentricity**2)
        / (1.0 + eccentricity * cos_true_anomaly)
        * (cos_true_anomaly * ecc_unit_vec + sin_true_anomaly * q_unit_vec)
    )

    if np.isnan(x).any():
        return np.array([np.nan, np.nan, np.nan])

    return x

if __name__ == "__main__":
    main()
