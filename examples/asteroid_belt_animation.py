"""
Demonstration on using the gravity simulator to simulate the asteroid belt

Note: Technically you can also create nice looking solar system animations by setting N = 0 and 
      expanding the axes limits.

Author: Ching Yin Ng
"""

from pathlib import Path
import sys
sys.path.append(str(Path(__file__).parent.parent))

import matplotlib.pyplot as plt
import numpy as np
import PIL
import rich.progress

from gravity_sim import GravitySimulator

N = 50000
FPS = 30
DPI = 200
N_FRAMES = 500

class Progress_bar(rich.progress.Progress):
    def __init__(self):
        super().__init__(
            rich.progress.BarColumn(),
            rich.progress.TextColumn("[green]{task.percentage:>3.0f}%"),
            rich.progress.TextColumn("•"),
            rich.progress.TimeElapsedColumn(),
            rich.progress.TextColumn("•"),
            rich.progress.TimeRemainingColumn(),
            "•[magenta] {task.completed}/{task.total}",
        )

def main():
    # ---------- Initialization ---------- #
    print("Initializing the system...", end="")
    grav_sim = GravitySimulator()
    system = grav_sim.create_system()

    system.load("solar_system")
    system.remove("Mercury")
    system.remove(name="Neptune")
    system.remove(name="Uranus")
    colors = [grav_sim.SOLAR_SYSTEM_COLORS[name] for name in system.objects_names]
    marker_sizes = [6.0, 1.5, 2.0, 1.5, 4.0, 3.5]

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

    # m = 0.0 if we assume asteroids are massless
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
    print("Done!")

    # ---------- Simulation and draw frames ---------- #
    file_path = Path(__file__).parent.parent / "gravity_sim" / "results"
    file_path.mkdir(parents=True, exist_ok=True)

    print()
    print("Simulating asteroid belt and drawing frames...")
    progress_bar = Progress_bar()
    with progress_bar:
        task = progress_bar.add_task("", total=N_FRAMES)
        for i in range(N_FRAMES):
            if i == 0:
                grav_sim.launch_simulation(
                    "rkf78",
                    grav_sim.years_to_days(5.0 / N_FRAMES),
                    tolerance=1e-6,
                    acceleration_method="massless",
                    storing_method="no_store",
                    no_print=True,
                    no_progress_bar=True,
                )
            else:
                grav_sim.resume_simulation(grav_sim.years_to_days(5.0 / N_FRAMES))

            # Drawing frame
            fig = plt.figure()
            plt.style.use("dark_background")
            ax = fig.add_subplot(111, projection="3d")

            xlim_min = -3
            xlim_max = 3
            ylim_min = -3
            ylim_max = 3
            zlim_min = -3
            zlim_max = 3

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

            # Plotting massive objects
            for j in range(massive_objects_count):
                ax.scatter(
                    grav_sim.simulator.x[j, 0],
                    grav_sim.simulator.x[j, 1],
                    grav_sim.simulator.x[j, 2],
                    marker="o",
                    label=system.objects_names[j],
                    color=colors[j],
                    s=marker_sizes[j],
                )

            # Plotting massless objects
            ax.scatter(
                grav_sim.simulator.x[massive_objects_count:, 0],
                grav_sim.simulator.x[massive_objects_count:, 1],
                grav_sim.simulator.x[massive_objects_count:, 2],
                color="white",
                marker=".",
                s=0.1,
                alpha=0.2,
            )

            # Add legend
            legend = ax.legend(loc="center right", bbox_to_anchor=(1.325, 0.5))
            legend.facecolor = "transparent"

            # Adjust figure for the legend
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
            plt.savefig(file_path / f"frames_{i:04d}.png", dpi=DPI)
            plt.close("all")

            progress_bar.update(task, advance=1)

    progress_bar.update(task, completed=N_FRAMES)
    plt.close("all")

    print()
    print("Combining frames to gif...")
    def frames_generator():
        for i in range(N_FRAMES):
            yield PIL.Image.open(file_path / f"frames_{i:04d}.png")

    frames = frames_generator()
    next(frames).save(
        file_path / "asteroid_belt.gif",
        save_all=True,
        append_images=frames,
        loop=0,
        duration=(1000 // FPS)
    )

    for i in range(N_FRAMES):
        (file_path / f"frames_{i:04d}.png").unlink()

    print(f"Output completed! Please check {file_path / 'asteroid_belt.gif'}")
    print()

if __name__ == "__main__":
    main()
