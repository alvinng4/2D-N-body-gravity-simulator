import csv
import math
from pathlib import Path
import sys
import warnings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter

import common


def plot_2d_trajectory(
    objects_count,
    sol_state,
    colors=None,
    labels=None,
    legend=False,
    xlabel="$x$ (AU)",
    ylabel="$y$ (AU)",
    marker="o",
    markersize=6,
):
    if colors is None:
        colors = [None for _ in range(objects_count)]
    if labels is None:
        labels = [None for _ in range(objects_count)]

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect="equal")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    for i in range(objects_count):
        traj = ax.plot(
            sol_state[:, i * 3],
            sol_state[:, 1 + i * 3],
            color=colors[i],
        )
        # Plot the last position with marker
        ax.plot(
            sol_state[-1, i * 3],
            sol_state[-1, 1 + i * 3],
            color=traj[0].get_color(),
            label=labels[i],
            marker=marker,
            markersize=markersize,
        )

    if legend:
        fig.legend(loc="center right", borderaxespad=0.2)
        fig.tight_layout()

    plt.show()


def plot_3d_trajectory(
    objects_count,
    sol_state,
    colors=None,
    labels=None,
    legend=False,
    xlabel="$x$ (AU)",
    ylabel="$y$ (AU)",
    zlabel="$z$ (AU)",
    marker="o",
    markersize=6,
):
    if colors is None:
        colors = [None for _ in range(objects_count)]
    if labels is None:
        labels = [None for _ in range(objects_count)]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.set_box_aspect([1.0, 1.0, 1.0])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)

    # Plot with solar system colors
    for i in range(objects_count):
        traj = ax.plot(
            sol_state[:, i * 3],
            sol_state[:, i * 3 + 1],
            sol_state[:, i * 3 + 2],
            color=colors[i],
        )
        # Plot the last position with marker
        ax.plot(
            sol_state[-1, i * 3],
            sol_state[-1, i * 3 + 1],
            sol_state[-1, i * 3 + 2],
            color=traj[0].get_color(),
            label=labels[i],
            marker=marker,
            markersize=markersize,
        )

    set_3d_axes_equal(ax)

    if legend:
        ax.legend(loc="center right", bbox_to_anchor=(1.325, 0.5))
        fig.subplots_adjust(
            right=0.7
        )  # Adjust the right boundary of the plot to make room for the legend

    fig.tight_layout()
    plt.show()


def animate_2d_traj_gif(
    objects_count,
    sol_state,
    fps,
    plot_every_nth_point,
    dpi,
    is_dynamic_axes=False,
    axes_lim=None,
    is_maintain_fixed_dt=True,
    traj_len=-1,
    colors=None,
    labels=None,
    legend=False,
    xlabel="$x$ (AU)",
    ylabel="$y$ (AU)",
    marker="o",
    markersize=6,
    file_name=None,
    file_path=None,
    sol_time=None,
):
    if colors is None:
        colors = [None for _ in range(objects_count)]
    if labels is None:
        labels = [None for _ in range(objects_count)]

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect="equal")

    if file_name is None and file_path is None:
        warnings.warn(
            "Either file_name or file_path should be provided. "
            + 'Setting file_name to "2d_traj".'
        )
        file_name = "2d_traj"

    if file_path is None:
        file_path = Path(__file__).parent / "results"
        file_path.mkdir(parents=True, exist_ok=True)
        if not file_name.endswith(".gif"):
            file_name += ".gif"
        file_path /= file_name

    if not is_dynamic_axes:
        if axes_lim is not None:
            xlim_min = axes_lim[0]
            xlim_max = axes_lim[1]
            ylim_min = axes_lim[2]
            ylim_max = axes_lim[3]
        else:
            xlim_min, xlim_max = get_axes_lim(2, objects_count, sol_state)
            ylim_min = xlim_min
            ylim_max = xlim_max

    data_size = len(sol_state)
    writer = PillowWriter(fps=fps)
    progress_bar = common.Progress_bar()
    if not is_maintain_fixed_dt:
        with writer.saving(fig, file_path, dpi):
            with progress_bar:
                # Plot once every nth point
                for i in progress_bar.track(range(data_size)):
                    if i % plot_every_nth_point != 0 and i != (data_size - 1):
                        continue

                    # Plot the trajectory from the beginning to current position
                    if traj_len == -1:
                        start_index = 0
                    else:
                        start_index = np.clip(i + 1 - traj_len, 0, None)

                    for j in range(objects_count):
                        traj = ax.plot(
                            sol_state[start_index : (i + 1), j * 3],
                            sol_state[start_index : (i + 1), j * 3 + 1],
                            color=colors[j],
                        )
                        # Plot the current position as a filled circle
                        ax.plot(
                            sol_state[i, j * 3],
                            sol_state[i, j * 3 + 1],
                            color=traj[0].get_color(),
                            label=labels[j],
                            marker=marker,
                            markersize=markersize,
                        )
                    # Add legend at the beginning
                    if legend:
                        fig.legend(loc="center right", borderaxespad=0.2)

                    # Set axis labels before capturing the frame
                    ax.set_xlabel(xlabel)
                    ax.set_ylabel(ylabel)

                    if not is_dynamic_axes:
                        ax.set_xlim([xlim_min, xlim_max])
                        ax.set_ylim([ylim_min, ylim_max])

                    # Capture the frame
                    writer.grab_frame()

                    # Clear the plot to prepare for the next frame
                    ax.clear()

            print("Saving the file...")

    else:
        # Attempt to maintain fixed dt for the animation
        with writer.saving(fig, file_path, dpi):
            with progress_bar:
                frame_size = math.floor(data_size / plot_every_nth_point) + 1
                plot_time = np.linspace(
                    sol_time[0],
                    sol_time[-1],
                    frame_size,
                )
                # Plot once every nth point
                for i in progress_bar.track(range(frame_size)):
                    # Search the index with the closest value of time
                    index = np.searchsorted(sol_time, plot_time[i])

                    if traj_len == -1:
                        start_index = 0
                    else:
                        start_index = np.searchsorted(
                            sol_time, plot_time[np.clip(i - traj_len, 0, None)]
                        )

                    # Plot the trajectory from the beginning to current position
                    for j in range(objects_count):
                        traj = ax.plot(
                            sol_state[(start_index + 1) : (index + 1), j * 3],
                            sol_state[(start_index + 1) : (index + 1), j * 3 + 1],
                            color=colors[j],
                        )
                        # Plot the current position as a filled circle
                        ax.plot(
                            sol_state[index, j * 3],
                            sol_state[index, j * 3 + 1],
                            color=traj[0].get_color(),
                            label=labels[j],
                            marker=marker,
                            markersize=markersize,
                        )
                    # Add legend at the beginning
                    if legend:
                        fig.legend(loc="center right", borderaxespad=0.2)

                    # Set axis labels before capturing the frame
                    ax.set_xlabel(xlabel)
                    ax.set_ylabel(ylabel)

                    if not is_dynamic_axes:
                        ax.set_xlim([xlim_min, xlim_max])
                        ax.set_ylim([ylim_min, ylim_max])

                    # Capture the frame
                    writer.grab_frame()

                    # Clear the plot to prepare for the next frame
                    ax.clear()

            print("Saving the file...")

    plt.close("all")
    print(f"Output completed! Please check {file_path}")


def animate_3d_traj_gif(
    objects_count: int,
    sol_state: np.ndarray,
    fps: int,
    plot_every_nth_point,
    dpi,
    is_dynamic_axes=False,
    axes_lim=None,
    is_maintain_fixed_dt=True,
    traj_len=-1,
    colors=None,
    labels=None,
    legend=False,
    xlabel="$x$ (AU)",
    ylabel="$y$ (AU)",
    zlabel="$z$ (AU)",
    marker="o",
    markersize=6,
    file_name=None,
    file_path=None,
    sol_time=None,
) -> None:
    if colors is None:
        colors = [None for _ in range(objects_count)]
    if labels is None:
        labels = [None for _ in range(objects_count)]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    if file_name is None and file_path is None:
        warnings.warn(
            "Either file_name or file_path should be provided. "
            + 'Setting file_name to "3d_traj".'
        )
        file_name = "3d_traj"

    if file_path is None:
        file_path = Path(__file__).parent / "results"
        file_path.mkdir(parents=True, exist_ok=True)
        if not file_name.endswith(".gif"):
            file_name += ".gif"
        file_path /= file_name

    if not is_dynamic_axes:
        if axes_lim is not None:
            xlim_min = axes_lim[0]
            xlim_max = axes_lim[1]
            ylim_min = axes_lim[2]
            ylim_max = axes_lim[3]
            zlim_min = axes_lim[4]
            zlim_max = axes_lim[5]
        else:
            xlim_min, xlim_max = get_axes_lim(3, objects_count, sol_state)
            ylim_min = xlim_min
            ylim_max = xlim_max
            zlim_min = xlim_min
            zlim_max = xlim_max

    data_size = len(sol_state)
    writer = PillowWriter(fps=fps)
    progress_bar = common.Progress_bar()
    if not is_maintain_fixed_dt:
        with writer.saving(fig, file_path, dpi):
            with progress_bar:
                # Plot once every nth point
                for i in progress_bar.track(range(data_size)):
                    if i % plot_every_nth_point != 0 and i != (data_size - 1):
                        continue

                    if traj_len == -1:
                        start_index = 0
                    else:
                        start_index = np.clip(i + 1 - traj_len, 0, None)

                    # Plot the trajectory from the beginning to current position
                    for j in range(objects_count):
                        traj = ax.plot(
                            sol_state[start_index : (i + 1), j * 3],
                            sol_state[start_index : (i + 1), j * 3 + 1],
                            sol_state[start_index : (i + 1), j * 3 + 2],
                            color=colors[j],
                        )
                        # Plot the current position as a filled circle
                        ax.plot(
                            sol_state[i, j * 3],
                            sol_state[i, j * 3 + 1],
                            sol_state[i, j * 3 + 2],
                            color=traj[0].get_color(),
                            label=labels[j],
                            marker=marker,
                            markersize=markersize,
                        )

                    # Add legend
                    if legend:
                        ax.legend(loc="center right", bbox_to_anchor=(1.325, 0.5))

                        # Adjust figure for the legend
                        if i == 0:
                            fig.subplots_adjust(right=0.7)
                            fig.tight_layout()

                    # Set axis labels and setting 3d axes scale before capturing the frame
                    ax.set_xlabel(xlabel)
                    ax.set_ylabel(ylabel)
                    ax.set_zlabel(zlabel)

                    if not is_dynamic_axes:
                        ax.set_xlim3d([xlim_min, xlim_max])
                        ax.set_ylim3d([ylim_min, ylim_max])
                        ax.set_zlim3d([zlim_min, zlim_max])
                    else:
                        set_3d_axes_equal(ax)

                    # Set equal aspect ratio to prevent distortion
                    ax.set_aspect("equal")

                    # Capture the frame
                    writer.grab_frame()

                    # Clear the plot to prepare for the next frame
                    ax.clear()

            print("Saving the file...")

    else:
        # Attempt to maintain fixed dt for the animation
        frame_size = math.floor(data_size / plot_every_nth_point) + 1
        plot_time = np.linspace(
            sol_time[0],
            sol_time[-1],
            frame_size,
        )
        with writer.saving(fig, file_path, dpi):
            with progress_bar:
                # Plot once every nth point
                for i in progress_bar.track(range(frame_size)):
                    # Search the index with the closest value of time
                    index = np.searchsorted(sol_time, plot_time[i])

                    if traj_len == -1:
                        start_index = 0
                    else:
                        start_index = np.searchsorted(
                            sol_time, plot_time[np.clip(i - traj_len, 0, None)]
                        )

                    # Plot the trajectory from the beginning to current position
                    for j in range(objects_count):
                        traj = ax.plot(
                            sol_state[start_index : (index + 1), j * 3],
                            sol_state[start_index : (index + 1), j * 3 + 1],
                            sol_state[start_index : (index + 1), j * 3 + 2],
                            color=colors[j],
                        )
                        # Plot the current position as a filled circle
                        ax.plot(
                            sol_state[index, j * 3],
                            sol_state[index, j * 3 + 1],
                            sol_state[index, j * 3 + 2],
                            color=traj[0].get_color(),
                            label=labels[j],
                            marker=marker,
                            markersize=markersize,
                        )

                    # Add legend
                    if legend:
                        ax.legend(loc="center right", bbox_to_anchor=(1.325, 0.5))

                        # Adjust figure for the legend
                        if i == 0:
                            fig.subplots_adjust(right=0.7)
                            fig.tight_layout()

                    # Set axis labels and setting 3d axes scale before capturing the frame
                    ax.set_xlabel("$x$ (AU)")
                    ax.set_ylabel("$y$ (AU)")
                    ax.set_zlabel("$z$ (AU)")

                    if not is_dynamic_axes:
                        ax.set_xlim3d([xlim_min, xlim_max])
                        ax.set_ylim3d([ylim_min, ylim_max])
                        ax.set_zlim3d([zlim_min, zlim_max])
                    else:
                        set_3d_axes_equal(ax)

                    # Set equal aspect ratio to prevent distortion
                    ax.set_aspect("equal")

                    # Capture the frame
                    writer.grab_frame()

                    # Clear the plot to prepare for the next frame
                    ax.clear()

            print("Saving the file...")

    plt.close("all")
    print(f"Output completed! Please check {file_path}")
    print()


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


def get_axes_lim(
    dim: int, objects_count: int, sol_state: np.ndarray
) -> tuple[float, float]:
    """
    Calculate axes limits if user chooses static axes limits

    Parameters
    ----------
    dim : int
        Dimension of the system
    objects_count : int
        Number of objects in the system
    sol_state : np.ndarray
        State of the system

    Returns
    -------
    lim_min : float
        Minimum axes limits
    lim_max : float
        Maximum axes limits
    """
    if dim == 2:
        xmax = np.max(sol_state[:, 0 * 3])
        xmin = np.min(sol_state[:, 0 * 3])
        ymax = np.max(sol_state[:, 0 * 3 + 1])
        ymin = np.min(sol_state[:, 0 * 3 + 1])

        for i in range(objects_count):
            xmax = max(xmax, np.max(sol_state[:, i * 3]))
            xmin = min(xmin, np.min(sol_state[:, i * 3]))
            ymax = max(ymax, np.max(sol_state[:, i * 3 + 1]))
            ymin = min(ymin, np.min(sol_state[:, i * 3 + 1]))

        lim_max = max(xmax, ymax)
        lim_min = min(xmin, ymin)

        if lim_max > 0:
            lim_max *= 1.1
        else:
            lim_max *= 0.9

        if lim_min < 0:
            lim_min *= 1.1
        else:
            lim_min *= 0.9

    elif dim == 3:
        xmax = np.max(sol_state[:, 0 * 3])
        xmin = np.min(sol_state[:, 0 * 3])
        ymax = np.max(sol_state[:, 0 * 3 + 1])
        ymin = np.min(sol_state[:, 0 * 3 + 1])
        zmax = np.max(sol_state[:, 0 * 3 + 2])
        zmin = np.min(sol_state[:, 0 * 3 + 2])

        for i in range(1, objects_count):
            xmax = max(xmax, np.max(sol_state[:, i * 3]))
            xmin = min(xmin, np.min(sol_state[:, i * 3]))
            ymax = max(ymax, np.max(sol_state[:, i * 3 + 1]))
            ymin = min(ymin, np.min(sol_state[:, i * 3 + 1]))
            zmax = max(zmax, np.max(sol_state[:, i * 3 + 2]))
            zmin = min(zmin, np.min(sol_state[:, i * 3 + 2]))

        lim_max = max(xmax, ymax, zmax)
        lim_min = min(xmin, ymin, zmin)

        if lim_max > 0:
            lim_max *= 1.1
        else:
            lim_max *= 0.9

        if lim_min < 0:
            lim_min *= 1.1
        else:
            lim_min *= 0.9

    else:
        raise ValueError("Invalid dimension")

    return lim_min, lim_max


def plot_rel_energy(
    energy: np.ndarray,
    sol_time: np.ndarray,
    title="Relative energy error against time",
    xlabel=f"Time",
    ylabel="$|(E(t)-E_0)/E_0|$",
) -> None:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.semilogy(
        sol_time,
        np.abs((energy - energy[0]) / energy[0]),
    )
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    plt.show()


def plot_tot_energy(
    energy: np.ndarray,
    sol_time: np.ndarray,
    title="Total energy against time",
    xlabel="Time",
    ylabel="$E(t)$",
):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.semilogy(sol_time, np.abs(energy))
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    plt.show()


def plot_rel_linear_momentum(
    linear_momentum: np.ndarray,
    sol_time: np.ndarray,
    title: str = "Relative linear momentum error against time",
    xlabel: str = "Time",
    ylabel: str = "$|(p(t)-p_0)/p_0|$",
) -> None:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.semilogy(
        sol_time,
        np.abs((linear_momentum - linear_momentum[0]) / linear_momentum[0]),
    )
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    plt.show()


def plot_rel_angular_momentum(
    angular_momentum: np.ndarray,
    sol_time: np.ndarray,
    title: str = "Relative angular momentum error against time",
    xlabel: str = "Time",
    ylabel: str = "$|(L(t)-L_0)/L_0|$",
) -> None:
    if angular_momentum[0] == 0.0:
        warnings.warn("Initial angular momentum is zero.")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.semilogy(
        sol_time,
        np.abs((angular_momentum - angular_momentum[0]) / angular_momentum[0]),
    )
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    plt.show()


def plot_dt(
    sol_dt: np.ndarray,
    sol_time: np.ndarray,
    title: str = "dt against time",
    xlabel: str = "Time",
    ylabel: str = "$dt$ (days)",
    yscale: str = "log",
    marker_size: float = 0.1,
) -> None:
    """
    Plot dt(days)
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(sol_time, sol_dt, s=marker_size)
    ax.set_yscale(yscale)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    plt.show()


def plot_eccentricity(
    eccentricity: np.ndarray,
    sol_time: np.ndarray,
    colors=None,
    labels=None,
    legend=False,
    title: str = "Eccentricity against time",
    xlabel: str = "Time",
    ylabel: str = "Eccentricity",
) -> None:
    objects_count = len(eccentricity[0])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(objects_count):
        ax.plot(sol_time, eccentricity[:, i], color=colors[i], label=labels[i])
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    if legend:
        fig.legend(loc="center right", borderaxespad=0.2)
        fig.tight_layout()

    plt.show()


def plot_inclination(
    inclination: np.ndarray,
    sol_time: np.ndarray,
    colors=None,
    labels=None,
    legend=False,
    title: str = "Inclination against time",
    xlabel: str = "Time",
    ylabel: str = "Inclination (radians)",
) -> None:
    objects_count = len(inclination[0])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(objects_count):
        ax.plot(sol_time, inclination[:, i], color=colors[i], label=labels[i])
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    if legend:
        fig.legend(loc="center right", borderaxespad=0.2)
        fig.tight_layout()

    plt.show()


def plot_compare_rel_energy(grav_plot):
    filepaths_to_labels = {}
    filepaths_to_labels["current"] = input(
        "Enter plot label for current set of data in memory: "
    )

    read_folder_path = Path(__file__).parent / "results"

    while True:
        num_of_data = len(filepaths_to_labels)
        if num_of_data >= 2:
            if not common.get_bool(
                f"There are {num_of_data} sets of data. Continue adding new files?"
            ):
                break
        while True:
            print()
            read_file_path = input(
                "Enter absolute path of new file to compare, or the complete file name if it is inside gravity_plot/results: "
            )
            read_file_path = read_folder_path / read_file_path

            if read_file_path.is_file():
                label = input("Enter plot label for this file: ")
                filepaths_to_labels[read_file_path] = label
                break
            else:
                print("File not found! Please try again.")

    print()
    print("Plotting relative energy error...(Please check the window)")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(f"Time ({grav_plot.tf_unit})")
    ax.set_ylabel("$|(E(t)-E_0)/E_0|$")

    for filepath in filepaths_to_labels:
        if filepath == "current":
            ax.semilogy(
                grav_plot.sol_time_in_tf_unit,
                np.abs(
                    (grav_plot.simulator.energy - grav_plot.simulator.energy[0])
                    / grav_plot.simulator.energy[0]
                ),
                label=filepaths_to_labels[filepath],
            )
        else:
            try:
                with open(filepath, "r") as file:
                    reader = csv.reader(file)

                    # Allocate memory
                    sol_time = np.zeros(50000)
                    energy = np.zeros(50000)

                    i = 0
                    for row in reader:
                        sol_time[i] = row[0]
                        energy[i] = row[2]
                        i += 1

                        # Extending memory buffer
                        if i % 50000 == 0:
                            sol_time = np.concatenate((sol_time, np.zeros(50000)))
                            energy = np.concatenate((energy, np.zeros(50000)))

                    sol_time = sol_time[:i]
                    energy = energy[:i]

                    if grav_plot.tf_unit == "years":
                        sol_time /= grav_plot.DAYS_PER_YEAR

                    ax.semilogy(
                        sol_time,
                        np.abs((energy - energy[0]) / energy[0]),
                        label=filepaths_to_labels[filepath],
                    )

            except FileNotFoundError:
                sys.exit("Error: file is not found. Exiting the program")

    fig.legend(loc=7)
    fig.tight_layout()
    fig.subplots_adjust(right=0.8)
    plt.show()

    print()
