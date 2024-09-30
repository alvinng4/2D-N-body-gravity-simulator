"""
Demonstration on using the gravity simulator to simulate the asteroid belt
You will need to install the `Pillow` library for this script.

Warning: Do not run multiple instances of this program at the same time, unless you made copies
         of the whole directory. Otherwise, the final data may overwrite each other.
"""

from pathlib import Path
import sys
sys.path.append(str(Path(__file__).parent.parent))

import numpy as np
import PIL
import matplotlib.pyplot as plt
import rich.progress
from matplotlib.colors import LinearSegmentedColormap

from gravity_sim import GravitySimulator

N = 50000
FPS = 30
DPI = 200
N_FRAMES = 2000 + 1

# kpc^3 / (Msun * kyr^2)
GM_SUN = 132712440041.279419 # km^3 s^-2 M_sun^-1
G = GM_SUN * (365 * 24 * 3600 * 1000) ** 2 * (3.24077929e-17) ** 3

orange_white_cmap = LinearSegmentedColormap.from_list("orange_white", ["black", "orange", "white"])

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
    print("Initializing the system...", end="", flush=True)
    grav_sim = GravitySimulator()
    system = grav_sim.create_system()
    grav_sim.set_current_system(system)
    system.G = G
    
    # Galaxy 1
    # Central black hole
    system.add(
        x=np.array([-60.0, 0.0, 0.0]),
        v=np.array([0.0, 0.0, 0.0]),
        m=5e6,
        object_name="CentralBH_1",
    )
    # Disk with random radius
    rng = np.random.RandomState(seed=0)
    radius = rng.uniform(1, 30, size=(N - 1))
    argument_of_periapsis = rng.uniform(0, 2 * np.pi, size=(N - 1))
    long_asc_node = rng.uniform(0, 2 * np.pi, size=(N - 1))
    true_anomaly = rng.uniform(0, 2 * np.pi, size=(N - 1))
    masses = rng.uniform(1e-1, 1e1, size=(N - 1))
    for i in range(N - 1):
        system.add_keplerian(
            semi_major_axis=radius[i],
            eccentricity=0.0,
            inclination=0.0,
            longitude_of_ascending_node=long_asc_node[i],
            argument_of_periapsis=argument_of_periapsis[i],
            true_anomaly=true_anomaly[i],
            m=masses[i],
            primary_object_name="CentralBH_1",
        )

    # Galaxy 2
    # Central black hole
    system.add(
        x=np.array([60.0, 0.0, 0.0]),
        v=np.array([0.0, 0.0, 0.0]),
        m=5e6,
        object_name="CentralBH_2",
    )
    # Disk with random radius
    radius = rng.uniform(1, 30, size=(N - 1))
    argument_of_periapsis = rng.uniform(0, 2 * np.pi, size=(N - 1))
    long_asc_node = rng.uniform(0, 2 * np.pi, size=(N - 1))
    true_anomaly = rng.uniform(0, 2 * np.pi, size=(N - 1))
    masses = rng.uniform(1e-1, 1e1, size=(N - 1))
    for i in range(N - 1):
        system.add_keplerian(
            semi_major_axis=radius[i],
            eccentricity=0.0,
            inclination=0.0,
            longitude_of_ascending_node=long_asc_node[i],
            argument_of_periapsis=argument_of_periapsis[i],
            true_anomaly=true_anomaly[i],
            m=masses[i],
            primary_object_name="CentralBH_2",
        )

    system.center_of_mass_correction()
    # system.save("galaxy_collision")

    print("Done!")

    # ---------- Simulation and draw frames ---------- #
    file_path = Path(__file__).parent.parent / "gravity_sim" / "results"
    file_path.mkdir(parents=True, exist_ok=True)

    print()
    print("Simulating galaxy collision and drawing frames...")
    progress_bar = Progress_bar()
    with progress_bar:
        task = progress_bar.add_task("[cyan]Simulation progress", total=N_FRAMES)
        for i in range(N_FRAMES):
            if i == 0:
                grav_sim.launch_simulation(
                    integrator="rk4",
                    tf=0.0,
                    dt=1e4,
                    acceleration_method="barnes-hut",
                    storing_method="no_store",
                    no_print=True,
                    no_progress_bar=True,
                    softening_length=1.0,
                    barnes_hut_theta=1.0,
                )
            else:
                grav_sim.resume_simulation(
                    1e6,
                )

            # Drawing frame
            fig = plt.figure()
            plt.style.use("dark_background")
            ax = fig.add_subplot(111)

            xlim_min = -200
            xlim_max = 200
            ylim_min = -100
            ylim_max = 100

            ax.grid(False)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            ax.spines['top'].set_color('none')
            ax.spines['bottom'].set_color('none')
            ax.spines['left'].set_color('none')
            ax.spines['right'].set_color('none')

            # Galaxy 1
            hist1, _, _ = np.histogram2d(
                grav_sim.simulator.x[:N, 0],
                grav_sim.simulator.x[:N, 1],
                bins=500,
                range=[[xlim_min, xlim_max], [ylim_min, ylim_max]],
            )

            # Galaxy 2
            hist2, _, _ = np.histogram2d(
                grav_sim.simulator.x[N:, 0],
                grav_sim.simulator.x[N:, 1],
                bins=500,
                range=[[xlim_min, xlim_max], [ylim_min, ylim_max]],
            )

            if i == 0:
                print(hist1.max(), hist2.max())
            hist = np.clip(hist1 + hist2, 0.0, 25.0)
            
            ax.imshow(
                hist.T,
                origin="lower",
                cmap=orange_white_cmap,
                extent=[xlim_min, xlim_max, ylim_min, ylim_max],
            )   

            fig.tight_layout()

            ax.set_xlim(xlim_min, xlim_max)
            ax.set_ylim(ylim_min, ylim_max)

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
        file_path / "galaxy_collision.gif",
        save_all=True,
        append_images=frames,
        loop=0,
        duration=(1000 // FPS)
    )

    for i in range(N_FRAMES):
        (file_path / f"frames_{i:04d}.png").unlink()

    print(f"Output completed! Please check {file_path / 'galaxy_collision.gif'}")
    print()


if __name__ == "__main__":
    main()
