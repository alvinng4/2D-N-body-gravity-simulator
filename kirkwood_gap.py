"""
Demonstration on using the gravity simulator to simulate the Kirkwood gap.

Note: N = 1000 is enough to observe the gaps, but it may not be very clear.
Since the simulation is O(N), seting N = 50000 would take a few days to 
finish, and reducing N to 25000 will reduce the runtime by half. 

Warning: This script will take a lot of storage space on your computer (probably a few GBs).
         It will attempt to erase the data after the video is generated.
"""

import csv
from pathlib import Path
import sys

import numpy as np
import PIL
import matplotlib.pyplot as plt

from gravity_sim import GravitySimulator
from gravity_sim.common import get_bool
from gravity_sim.common import Progress_bar

N = 50000
FPS = 30
DPI = 300

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
    system.name = "kirkwood_gap_N50000"
    system.save()

    # ---------- Simulation ---------- #
    file_path = Path(__file__).parent / "gravity_sim" / "results"
    file_path.mkdir(parents=True, exist_ok=True)
    data_path = file_path / "kirkwood_gap_sim.csv"
    
    print("Simulating kirkwood gap...")
    # Store about 2000 points in total
    tf = grav_sim.years_to_days(5000000)
    dt = 180.0
    store_every_n = int((grav_sim.years_to_days(5000000) // dt) // 2000)
    grav_sim.launch_simulation(
        "whfast",
        tf=tf,
        dt=dt,
        store_every_n=store_every_n,
        acceleration="massless",
        flush=True,
        flush_results_path=str(data_path),
        no_print=True,
        kepler_auto_remove=True,
        debug=True,
    )

    # ---------- Data Analysis and drawing frames ---------- #

    # In the library, we use PillowWriter to generate animations.
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
                    ).reshape(asteroids_count, 3) - sun_x
                )
                asteroids_v = (
                    np.array(
                        row[inner_objects_count * 6 + (asteroids_count + outer_objects_count) * 3
                            : (inner_objects_count + asteroids_count) * 6 + outer_objects_count * 3
                        ]
                    ).reshape(asteroids_count, 3) - sun_v
                )

                eccentricity = calculate_eccentricity(asteroids_x, asteroids_v, 0.0, G, M)
                semi_major_axes = calculate_semi_major_axis(asteroids_x, asteroids_v, 0.0, G, M)
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

                # Due to removal of asteroids by Kepler auto clear,
                # we need to count the number of asteroids every row
                asteroids_count = int(len(row) // 6) - massive_objects_count
                sun_x = np.array(row[:3])
                asteroids_x = (
                    np.array(
                        row[
                            (inner_objects_count * 3) : (
                                inner_objects_count + asteroids_count
                            )
                            * 3
                        ]
                    ).reshape(asteroids_count, 3)
                    - sun_x
                )

                # Plotting the sun
                ax2.plot(
                    0.0,
                    0.0,
                    "o",
                    label=labels[0],
                    color=colors[0],
                    markersize=marker_sizes[0],
                )

                # Plotting the asteroids
                ax2.scatter(
                    asteroids_x[:, 0],
                    asteroids_x[:, 1],
                    color="white",
                    marker=".",
                    s=0.1,
                    alpha=0.2,
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
    visualization_frames = []
    for i in range(save_count_semi_major_axes):
        semi_major_axes_frames.append(
            PIL.Image.open(file_path / f"semi_major_axes_frames_{i:04d}.png")
        )

    for i in range(save_count_visualization):
        visualization_frames.append(
            PIL.Image.open(file_path / f"visualization_frames_{i:04d}.png")
        )

    semi_major_axes_frames[0].save(
        file_path / "Kirkwood_gap_semi_major_axes.gif",
        save_all=True,
        append_images=semi_major_axes_frames[1:],
        loop=0,
        duration=(1000 // FPS),
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


if __name__ == "__main__":
    main()
