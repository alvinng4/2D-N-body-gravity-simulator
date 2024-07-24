"""
Demonstration on using the gravity simulator to simulate the asteroid belt
"""

import os
from pathlib import Path

import numpy as np
import PIL
import matplotlib.pyplot as plt

import gravity_sim


def main():
    grav_sim = gravity_sim.GravitySimulator()
    system = grav_sim.create_system()
    system.load("solar_system")
    system.remove(name="Mercury")
    system.remove(name="Venus")

    massive_objects_count = system.objects_count

    N = 50000

    rng = np.random.default_rng()
    a = rng.uniform(2.1, 3.2, size=N)  # Semi-major axis in AU
    ecc = rng.uniform(0.0, 0.2, size=N)  # Eccentricity
    inc = rng.uniform(-0.5, 0.5, size=N)  # Inclination in radians
    raan = rng.uniform(
        0, 360, size=N
    )  # Right ascension of the ascending node in degrees
    argp = rng.uniform(0, 360, size=N)  # Argument of perigee in degrees
    nu = rng.uniform(0, 360, size=N)  # True anomaly in degrees

    for i in range(N):
        x, v = from_orbital_elements_to_cartesian(
            mp=1.0,
            ms=0.0,
            semimajor_axis=a[i],
            eccentricity=ecc[i],
            true_anomaly=np.radians(nu[i]),
            inclination=np.radians(inc[i]),
            argument_of_periapsis=np.radians(argp[i]),
            longitude_of_ascending_node=np.radians(raan[i]),
            G=system.G,
        )

        system.add(x, v, 0.0)

    system.center_of_mass_correction()

    fig = plt.figure()
    plt.style.use("dark_background")
    ax = fig.add_subplot(111, projection="3d")

    file_path = Path(__file__).parent / "gravity_sim" / "results"
    file_path.mkdir(parents=True, exist_ok=True)

    xlim_min = -3
    xlim_max = 3
    ylim_min = -3
    ylim_max = 3
    zlim_min = -3
    zlim_max = 3

    colors = ["orange", "skyblue", "red", "darkgoldenrod", "gold", "paleturquoise", "blue"]

    # In the package, we use PillowWriter to generate animation
    # However, for some reason, the PillowWriter run out of memory
    # in this case. Therefore, we save each frames as images and
    # combine them as gif instead.
    save_count = 0
    num_loops = 200
    for k in range(num_loops):
        print(f"Loop {k+1}/{num_loops}")
        grav_sim.simulator.launch_simulation(
            system,
            "rk4",
            grav_sim.years_to_days(0.02),
            dt=grav_sim.years_to_days(0.0005),
            store_every_n=20,
            acceleration="massless",
        )

        # Plot once every nth point
        for j in range(len(grav_sim.simulator.sol_state)):
            if j == 0:
                continue

            # Plot the trajectory from the beginning to current position
            for i in range(0, massive_objects_count):
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

                ax.plot(
                    grav_sim.simulator.sol_state[j, i * 3],
                    grav_sim.simulator.sol_state[j, 1 + i * 3],
                    "o",
                    label=system.objects_names[i],
                    color=colors[i],
                    markersize=6,
                )

            x = grav_sim.simulator.sol_state[
                j, (massive_objects_count * 3) : (grav_sim.simulator.objects_count * 3) : 3
            ]
            y = grav_sim.simulator.sol_state[
                j, (massive_objects_count * 3 + 1) : (grav_sim.simulator.objects_count * 3) : 3
            ]
            z = grav_sim.simulator.sol_state[
                j, (massive_objects_count * 3 + 2) : (grav_sim.simulator.objects_count * 3) : 3
            ]
            ax.scatter(
                x,
                y,
                z,
                color="white",
                marker=".",
                s=0.1,
                alpha=0.1,
            )

            # Add legend
            ax.legend(loc="center right", bbox_to_anchor=(1.325, 0.5))
            # Adjust figure for the legend
            if j == 0:
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
            plt.savefig(file_path / f"frames_{save_count:04d}.png", dpi=300)
            save_count += 1

            # Clear the plot to prepare for the next frame
            ax.clear()

        system = grav_sim.sol_state_to_system(objects_names=system.objects_names)

    plt.close("all")

    print("Combining frames to gif...")
    frames = []
    for i in range(save_count):
        frames.append(PIL.Image.open(file_path / f"frames_{i:04d}.png"))

    frames[0].save(
        file_path / "astroid_belt.gif",
        save_all=True,
        append_images=frames[1:],
        loop=0,
        duration=34,
    )

    for i in range(save_count):
        os.remove(file_path / f"frames_{i:04d}.png")
    print(f"Output completed! Please check {file_path}")
    print()


def from_orbital_elements_to_cartesian(
    mp,
    ms,
    semimajor_axis,
    eccentricity,
    true_anomaly,
    inclination,
    argument_of_periapsis,
    longitude_of_ascending_node,
    G,
):
    """

    Function that returns position and velocities computed from the input orbital
    elements. Angles in radians, inclination between 0 and 180

    Reference
    ---------
    https://github.com/MovingPlanetsAround/ABIE
    """

    cos_true_anomaly = np.cos(true_anomaly)
    sin_true_anomaly = np.sin(true_anomaly)

    cos_inclination = np.cos(inclination)
    sin_inclination = np.sin(inclination)

    cos_arg_per = np.cos(argument_of_periapsis)
    sin_arg_per = np.sin(argument_of_periapsis)

    cos_long_asc_nodes = np.cos(longitude_of_ascending_node)
    sin_long_asc_nodes = np.sin(longitude_of_ascending_node)

    ### e_vec is a unit vector directed towards periapsis ###
    e_vec_x = (
        cos_long_asc_nodes * cos_arg_per
        - sin_long_asc_nodes * sin_arg_per * cos_inclination
    )
    e_vec_y = (
        sin_long_asc_nodes * cos_arg_per
        + cos_long_asc_nodes * sin_arg_per * cos_inclination
    )
    e_vec_z = sin_arg_per * sin_inclination
    e_vec = np.array([e_vec_x, e_vec_y, e_vec_z])

    ### q is a unit vector perpendicular to e_vec and the orbital angular momentum vector ###
    q_vec_x = (
        -cos_long_asc_nodes * sin_arg_per
        - sin_long_asc_nodes * cos_arg_per * cos_inclination
    )
    q_vec_y = (
        -sin_long_asc_nodes * sin_arg_per
        + cos_long_asc_nodes * cos_arg_per * cos_inclination
    )
    q_vec_z = cos_arg_per * sin_inclination
    q_vec = np.array([q_vec_x, q_vec_y, q_vec_z])

    #    print 'alpha',alphax**2+alphay**2+alphaz**2 # For debugging; should be 1
    #    print 'beta',betax**2+betay**2+betaz**2 # For debugging; should be 1

    ### Relative position and velocity ###
    separation = (
        semimajor_axis
        * (1.0 - eccentricity**2)
        / (1.0 + eccentricity * cos_true_anomaly)
    )  # Compute the relative separation
    position_vector = (
        separation * cos_true_anomaly * e_vec + separation * sin_true_anomaly * q_vec
    )
    velocity_tilde = np.sqrt(
        G * (mp + ms) / (semimajor_axis * (1.0 - eccentricity**2))
    )  # Common factor
    velocity_vector = (
        -1.0 * velocity_tilde * sin_true_anomaly * e_vec
        + velocity_tilde * (eccentricity + cos_true_anomaly) * q_vec
    )

    return position_vector, velocity_vector


def set_3d_axes_equal(ax: plt.Axes) -> None:
    """
    Make axes of 3D plot have equal scale

    Parameters
    ----------
    ax : matplotlib axis
        The axis to set equal scale

    Reference
    ---------
    karlo, https://stackoverflow.com/questions/13685386/how-to-set-the-equal-aspect-ratio-for-all-axes-x-y-z
    """

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


if __name__ == "__main__":
    main()
