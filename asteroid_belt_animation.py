"""
Demonstration on using the gravity simulator to simulate the asteroid belt

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

N = 50000
FPS = 30
DPI = 300


def main():
    grav_sim = GravitySimulator()
    system = grav_sim.create_system()
    grav_sim.set_current_system(system)

    system.load("solar_system")
    system.remove(name="Mercury")
    system.remove(name="Venus")
    system.remove(name="Neptune")
    system.remove(name="Uranus")
    colors = [
        "orange",
        "skyblue",
        "red",
        "darkgoldenrod",
        "gold",
    ]
    marker_sizes = [6.0, 2.0, 1.5, 4.0, 3.5]

    #################################################
    # Adding a star to the system

    # Star 1
    # system.add_keplerian(
    #     semi_major_axis=5.5,
    #     eccentricity=0.7,
    #     inclination=0.05,
    #     argument_of_periapsis=0.07,
    #     longitude_of_ascending_node=0.07,
    #     true_anomaly=0.35,
    #     m=1.0,
    #     primary_object_name="Sun",
    #     object_name="Added Star",
    # )
    # colors.append("orange")
    # marker_sizes.append(6.0)

    # Star 2
    # system.add_keplerian(
    #     semi_major_axis=5.0,
    #     eccentricity=0.7,
    #     inclination=0.4,
    #     argument_of_periapsis=4.0,
    #     longitude_of_ascending_node=4.0,
    #     true_anomaly=4.0,
    #     m=1.0,
    #     primary_object_name="Sun",
    #     object_name="Added Star",
    # )
    # colors.append("orange")
    # marker_sizes.append(6.0)

    #################################################

    massive_objects_count = system.objects_count

    rng = np.random.default_rng()
    a = rng.uniform(2.1, 3.2, size=N)
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

    system.center_of_mass_correction()

    file_path = Path(__file__).parent / "gravity_sim" / "results"
    file_path.mkdir(parents=True, exist_ok=True)
    data_path = file_path / "asteroid_belt_sim.csv"

    print("Simulating asteroid belt...")
    grav_sim.launch_simulation(
        "rk4",
        grav_sim.years_to_days(5.0),
        dt=grav_sim.years_to_days(0.001),
        store_every_n=10,
        acceleration="massless",
        flush=True,
        flush_results_path=str(data_path),
        no_print=True,
    )

    # Draw frames
    print()
    print("Drawing frames...")
    fig = plt.figure()
    plt.style.use("dark_background")
    ax = fig.add_subplot(111, projection="3d")

    xlim_min = -3
    xlim_max = 3
    ylim_min = -3
    ylim_max = 3
    zlim_min = -3
    zlim_max = 3

    # In the library, we use PillowWriter to generate animations.
    # However, for some reason, the PillowWriter run out of memory
    # in this case. Therefore, we save each frames as images and
    # combine them as gif instead.
    save_count = 0
    new_field_lim = sys.maxsize
    while True:
        try:
            csv.field_size_limit(new_field_lim)
            break
        except OverflowError:
            new_field_lim = new_field_lim // 10

    with open(data_path, "r") as file:
        reader = csv.reader(file)
        for row in reader:
            if len(row) == 0 or row[0].startswith("#"):
                continue

            ax.grid(False)
            ax.xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
            ax.yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
            ax.zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            ax.zaxis.set_visible(False)
            ax.xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
            ax.yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
            ax.zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

            row = row[3:]
            row = list(map(float, row))

            for i in range(0, massive_objects_count):
                ax.plot(
                    np.array(row[i * 3]),
                    np.array(row[1 + i * 3]),
                    "o",
                    label=system.objects_names[i],
                    color=colors[i],
                    markersize=marker_sizes[i],
                )

            x = row[(massive_objects_count * 3) : (system.objects_count * 3) : 3]
            y = row[(massive_objects_count * 3 + 1) : (system.objects_count * 3) : 3]
            z = row[(massive_objects_count * 3 + 2) : (system.objects_count * 3) : 3]
            ax.scatter(
                x,
                y,
                z,
                color="white",
                marker=".",
                s=0.1,
                alpha=0.2,
            )

            # Add legend
            legend = ax.legend(loc="center right", bbox_to_anchor=(1.325, 0.5))
            legend.facecolor = "transparent"

            # Adjust figure for the legend
            if save_count == 0:
                fig.subplots_adjust(right=0.7)
                fig.tight_layout()

            # Set axis labels and setting 3d axes scale before capturing the frame
            # ax.set_xlabel("$x$ (AU)")
            # ax.set_ylabel("$y$ (AU)")
            # ax.set_zlabel("$z$ (AU)")

            ax.set_xlim3d([xlim_min, xlim_max])
            ax.set_ylim3d([ylim_min, ylim_max])
            ax.set_zlim3d([zlim_min, zlim_max])

            # Set equal aspect ratio to prevent distortion
            ax.set_aspect("equal")

            # Capture the frame
            plt.savefig(file_path / f"frames_{save_count:04d}.png", dpi=DPI)
            save_count += 1

            # Clear the plot to prepare for the next frame
            ax.clear()

        plt.close("all")

    print()
    print("Combining frames to gif...")
    frames = []
    for i in range(save_count):
        frames.append(PIL.Image.open(file_path / f"frames_{i:04d}.png"))

    frames[0].save(
        file_path / "asteroid_belt.gif",
        save_all=True,
        append_images=frames[1:],
        loop=0,
        duration=(1000 // FPS),
    )

    for i in range(save_count):
        (file_path / f"frames_{i:04d}.png").unlink()
    data_path.unlink()

    print(f"Output completed! Please check {file_path / 'asteroid_belt.gif'}")
    print()


if __name__ == "__main__":
    main()
