from pathlib import Path
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
import PIL

from . import utils

SOLAR_SYSTEM_COLORS = {
    "Sun": "orange",
    "Mercury": "slategrey",
    "Venus": "wheat",
    "Earth": "skyblue",
    "Mars": "red",
    "Jupiter": "darkgoldenrod",
    "Saturn": "gold",
    "Uranus": "paleturquoise",
    "Neptune": "blue",
    "Moon": "grey",
    "Pluto": None,
    "Ceres": None,
    "Vesta": None,
}


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

    x_limits = ax.get_xlim3d()  # type: ignore
    y_limits = ax.get_ylim3d()  # type: ignore
    z_limits = ax.get_zlim3d()  # type: ignore

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])  # type: ignore
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])  # type: ignore
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])  # type: ignore


def plot_2d_trajectory(
    sol_state: np.ndarray,
    colors: Optional[list[str]] = None,
    labels: Optional[list[str]] = None,
    legend: bool = False,
    xlabel: str = "$x$ (AU)",
    ylabel: str = "$y$ (AU)",
    title: Optional[str] = None,
    marker: str = "o",
    markersize: int = 6,
    save_fig: bool = False,
    save_fig_path: Optional[str | Path] = None,
) -> None:
    """Plot 2D trajectory of objects

    Parameters
    ----------
    sol_state : np.ndarray
        Solution state of the system
    colors : Optional[list[str]], optional
        Colors of the trajectories, by default None
    labels : Optional[list[str]], optional
        Labels of the objects, used for legend, by default None
    legend : bool, optional
        Flag to check whether to show legend, by default False
    xlabel : str, optional
        Label of x-axis, by default "$x$ (AU)"
    ylabel : str, optional
        Label of y-axis, by default "$y$ (AU)"
    marker : str, optional
        Marker for the last position, by default "o"
    markersize : int, optional
        Marker size for the last position, by default 6
    save_fig : bool, optional
        Flag to check whether to save the figure, by default False
    save_fig_path : Optional[str], optional
        Path to save the figure, by default None
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect="equal")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    for i in range(sol_state.shape[1] // 6):
        if colors is not None:
            colors_i = colors[i]
        else:
            colors_i = None

        if labels is not None:
            labels_i = labels[i]
        else:
            labels_i = None

        traj = ax.plot(
            sol_state[:, i * 3],
            sol_state[:, 1 + i * 3],
            color=colors_i,
        )
        # Plot the last position with marker
        ax.plot(
            sol_state[-1, i * 3],
            sol_state[-1, 1 + i * 3],
            color=traj[0].get_color(),
            label=labels_i,
            marker=marker,
            markersize=markersize,
        )

    if title is not None:
        ax.set_title(title)

    if legend:
        fig.legend(loc="center right", borderaxespad=0.2)
        fig.tight_layout()

    if save_fig:
        plt.savefig(save_fig_path, dpi=300)
    else:
        plt.show()
    plt.close("all")


def plot_3d_trajectory(
    sol_state: np.ndarray,
    colors: Optional[list[str]] = None,
    labels: Optional[list[str]] = None,
    legend: bool = False,
    xlabel: str = "$x$ (AU)",
    ylabel: str = "$y$ (AU)",
    zlabel: str = "$z$ (AU)",
    title: Optional[str] = None,
    marker: str = "o",
    markersize: int = 6,
    save_fig: bool = False,
    save_fig_path: Optional[str | Path] = None,
) -> None:
    """Plot 3D trajectory of objects

    Parameters
    ----------
    sol_state : np.ndarray
        Solution state of the system
    colors : Optional[list[str]], optional
        Colors of the trajectories, by default None
    labels : Optional[list[str]], optional
        Labels of the objects, used for legend, by default None
    legend : bool, optional
        Flag to check whether to show legend, by default False
    xlabel : str, optional
        Label of x-axis, by default "$x$ (AU)"
    ylabel : str, optional
        Label of y-axis, by default "$y$ (AU)"
    zlabel : str, optional
        Label of z-axis, by default "$z$ (AU)"
    marker : str, optional
        Marker for the last position, by default "o"
    markersize : int, optional
        Marker size for the last position, by default 6
    save_fig : bool, optional
        Flag to check whether to save the figure, by default False
    save_fig_path : Optional[str], optional
        Path to save the figure, by default None
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.set_box_aspect([1.0, 1.0, 1.0])  # type: ignore
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)  # type: ignore

    for i in range(sol_state.shape[1] // 6):
        if colors is not None:
            colors_i = colors[i]
        else:
            colors_i = None

        if labels is not None:
            labels_i = labels[i]
        else:
            labels_i = None

        traj = ax.plot(
            sol_state[:, i * 3],
            sol_state[:, i * 3 + 1],
            sol_state[:, i * 3 + 2],
            color=colors_i,
        )
        # Plot the last position with marker
        ax.plot(
            sol_state[-1, i * 3],
            sol_state[-1, i * 3 + 1],
            sol_state[-1, i * 3 + 2],
            color=traj[0].get_color(),
            label=labels_i,
            marker=marker,
            markersize=markersize,
        )

    set_3d_axes_equal(ax)

    if title is not None:
        ax.set_title(title)

    if legend:
        ax.legend(loc="center right", bbox_to_anchor=(1.325, 0.5))
        fig.subplots_adjust(
            right=0.7
        )  # Adjust the right boundary of the plot to make room for the legend

    fig.tight_layout()

    if save_fig:
        plt.savefig(save_fig_path, dpi=300)
    else:
        plt.show()
    plt.close("all")


def _animation_get_axes_lim(
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


def animate_2d_traj_gif(
    file_path: str | Path,
    sol_state: np.ndarray,
    fps: int | float,
    plotting_freq: int,
    dpi: int | float,
    is_dynamic_axes: bool = False,
    axes_limits: Optional[list[float]] = None,
    is_maintain_fixed_dt: bool = False,
    sol_time: Optional[np.ndarray] = None,
    traj_len: int = -1,
    colors: Optional[list[str]] = None,
    labels: Optional[list[str]] = None,
    legend: bool = False,
    xlabel: str = "$x$ (AU)",
    ylabel: str = "$y$ (AU)",
    marker: str = "o",
    markersize: int = 6,
) -> None:
    """Animate 2D trajectory of objects and save as GIF

    Parameters
    ----------
    file_path : str | Path
        Path to save the GIF
    sol_state : np.ndarray
        Solution state of the system
    fps : int
        Frames per second
    plotting_freq : int
        Frequency of plotting the data
    dpi : int
        Dots per inch
    is_dynamic_axes : bool, optional
        Flag to check whether to use dynamic axes, by default False
    axes_limits : Optional[list[float]], optional
        Limits of the axes, by default None
    is_maintain_fixed_dt : bool, optional
        Flag to check whether to maintain fixed time interval, by default False
    sol_time : Optional[np.ndarray], optional
        Solution time, by default None
    traj_len : int, optional
        Length of the trajectory to be plotted, by default -1
    colors : Optional[list[str]], optional
        Colors of the trajectories, by default None
    labels : Optional[list[str]], optional
        Labels of the objects, used for legend, by default None
    legend : bool, optional
        Flag to check whether to show legend, by default False
    xlabel : str, optional
        Label of x-axis, by default "$x$ (AU)"
    ylabel : str, optional
        Label of y-axis, by default "$y$ (AU)"
    marker : str, optional
        Marker for the last position, by default "o"
    markersize : int, optional
        Marker size for the last position, by default 6

    Raises
    ------
    TypeError
        Invalid file path type
    ValueError
        Invalid axes limits provided
    """
    if isinstance(file_path, str):
        Path(file_path).parent.mkdir(parents=True, exist_ok=True)
    elif isinstance(file_path, Path):
        file_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        raise TypeError("Invalid file path type")

    objects_count = sol_state.shape[1] // 6

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect="equal")

    if not is_dynamic_axes:
        if axes_limits is not None:
            if len(axes_limits) != 4:
                raise ValueError(
                    "Invalid axes limits, should be [xmin, xmax, ymin, ymax]"
                )
            xlim_min = axes_limits[0]
            xlim_max = axes_limits[1]
            ylim_min = axes_limits[2]
            ylim_max = axes_limits[3]
        else:
            xlim_min, xlim_max = _animation_get_axes_lim(2, objects_count, sol_state)
            ylim_min, ylim_max = xlim_min, xlim_max

        data_size = len(sol_state)
        progress_bar = utils.Progress_bar()
        num_frames_count = 0

        if not is_maintain_fixed_dt:
            with progress_bar:
                for i in progress_bar.track(range(data_size)):
                    if i % plotting_freq != 0:
                        continue

                    # Plot the trajectory from the beginning to current position
                    if traj_len == -1:
                        start_index = 0
                    else:
                        start_index = np.clip(i + 1 - traj_len, 0, None)

                    for j in range(objects_count):
                        colors_j = None
                        labels_j = None
                        if colors is not None:
                            colors_j = colors[j]
                        if labels is not None:
                            labels_j = labels[j]
                        traj = ax.plot(
                            sol_state[start_index : (i + 1), j * 3],
                            sol_state[start_index : (i + 1), j * 3 + 1],
                            color=colors_j,
                        )
                        # Plot the current position as a filled circle
                        ax.plot(
                            sol_state[i, j * 3],
                            sol_state[i, j * 3 + 1],
                            color=traj[0].get_color(),
                            label=labels_j,
                            marker=marker,
                            markersize=markersize,
                        )
                    # Add legend at the beginning
                    # Warning: if we keep adding legend, the memory consumption will keep increasing
                    if legend and i == 0:
                        fig.legend(loc="center right", borderaxespad=0.2)

                    # Set axis labels before capturing the frame
                    ax.set_xlabel(xlabel)
                    ax.set_ylabel(ylabel)

                    if not is_dynamic_axes:
                        ax.set_xlim((xlim_min, xlim_max))
                        ax.set_ylim((ylim_min, ylim_max))

                    # Capture the frame
                    plt.savefig(
                        Path(file_path).parent / f"frames_{num_frames_count:06d}.png",
                        dpi=dpi,
                    )
                    num_frames_count += 1

                    # Clear the plot to prepare for the next frame
                    ax.clear()
        else:
            if sol_time is None:
                raise ValueError("Solution time is required to maintain fixed dt")

            # Attempt to maintain fixed dt for the animation
            frame_size = int(data_size / plotting_freq) + 1
            plot_time = np.linspace(
                sol_time[0],
                sol_time[-1],
                frame_size,
            )
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
                        colors_j = None
                        labels_j = None
                        if colors is not None:
                            colors_j = colors[j]
                        if labels is not None:
                            labels_j = labels[j]
                        traj = ax.plot(
                            sol_state[(start_index + 1) : (index + 1), j * 3],
                            sol_state[(start_index + 1) : (index + 1), j * 3 + 1],
                            color=colors_j,
                        )
                        # Plot the current position as a filled circle
                        ax.plot(
                            sol_state[index, j * 3],
                            sol_state[index, j * 3 + 1],
                            color=traj[0].get_color(),
                            label=labels_j,
                            marker=marker,
                            markersize=markersize,
                        )
                    # Add legend at the beginning
                    # Warning: if we keep adding legend, the memory consumption will keep increasing
                    if legend and i == 0:
                        fig.legend(loc="center right", borderaxespad=0.2)

                    # Set axis labels before capturing the frame
                    ax.set_xlabel(xlabel)
                    ax.set_ylabel(ylabel)

                    if not is_dynamic_axes:
                        ax.set_xlim((xlim_min, xlim_max))
                        ax.set_ylim((ylim_min, ylim_max))

                    # Capture the frame
                    plt.savefig(
                        Path(file_path).parent / f"frames_{num_frames_count:06d}.png",
                        dpi=dpi,
                        format="png",
                    )
                    num_frames_count += 1

                    # Clear the plot to prepare for the next frame
                    ax.clear()

    print("Combining frames to gif...")

    def frames_generator(num_frames_count):
        for i in range(num_frames_count):
            yield PIL.Image.open(file_path.parent / f"frames_{i:06d}.png")

    frames = frames_generator(num_frames_count)
    next(frames).save(
        file_path,
        format="GIF",
        append_images=frames,
        save_all=True,
        duration=(1000 / fps),
        loop=0,
    )

    plt.close("all")

    for i in range(num_frames_count):
        (Path(file_path).parent / f"frames_{i:06d}.png").unlink()


def animate_3d_traj_gif(
    file_path: str | Path,
    sol_state: np.ndarray,
    fps: int | float,
    plotting_freq: int,
    dpi: int | float,
    is_dynamic_axes: bool = False,
    axes_limits: Optional[list[float]] = None,
    is_maintain_fixed_dt: bool = False,
    sol_time: Optional[np.ndarray] = None,
    traj_len: int = -1,
    colors: Optional[list[str]] = None,
    labels: Optional[list[str]] = None,
    legend: bool = False,
    xlabel: str = "$x$ (AU)",
    ylabel: str = "$y$ (AU)",
    zlabel: str = "$z$ (AU)",
    marker: str = "o",
    markersize: int = 6,
) -> None:
    """Animate 3D trajectory of objects and save as GIF

    Parameters
    ----------
    file_path : str | Path
        Path to save the GIF
    sol_state : np.ndarray
        Solution state of the system
    fps : int
        Frames per second
    plotting_freq : int
        Frequency of plotting the data
    dpi : int
        Dots per inch
    is_dynamic_axes : bool, optional
        Flag to check whether to use dynamic axes, by default False
    axes_limits : Optional[list[float]], optional
        Limits of the axes, by default None
    is_maintain_fixed_dt : bool, optional
        Flag to check whether to maintain fixed time interval, by default False
    sol_time : Optional[np.ndarray], optional
        Solution time, by default None
    traj_len : int, optional
        Length of the trajectory to be plotted, by default -1
    colors : Optional[list[str]], optional
        Colors of the trajectories, by default None
    labels : Optional[list[str]], optional
        Labels of the objects, used for legend, by default None
    legend : bool, optional
        Flag to check whether to show legend, by default False
    xlabel : str, optional
        Label of x-axis, by default "$x$ (AU)"
    ylabel : str, optional
        Label of y-axis, by default "$y$ (AU)"
    zlabel : str, optional
        Label of z-axis, by default "$z$ (AU)"
    marker : str, optional
        Marker for the last position, by default "o"
    markersize : int, optional
        Marker size for the last position, by default 6

    Raises
    ------
    TypeError
        Invalid file path type
    ValueError
        Invalid axes limits provided
    """
    if isinstance(file_path, str):
        Path(file_path).parent.mkdir(parents=True, exist_ok=True)
    elif isinstance(file_path, Path):
        file_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        raise TypeError("Invalid file path type")

    objects_count = sol_state.shape[1] // 6

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    if not is_dynamic_axes:
        if axes_limits is not None:
            if len(axes_limits) != 6:
                raise ValueError(
                    "Invalid axes limits, should be [xmin, xmax, ymin, ymax, zmin, zmax]"
                )
            xlim_min = axes_limits[0]
            xlim_max = axes_limits[1]
            ylim_min = axes_limits[2]
            ylim_max = axes_limits[3]
            zlim_min = axes_limits[4]
            zlim_max = axes_limits[5]
        else:
            xlim_min, xlim_max = _animation_get_axes_lim(3, objects_count, sol_state)
            ylim_min, ylim_max = xlim_min, xlim_max
            zlim_min, zlim_max = xlim_min, xlim_max

    data_size = len(sol_state)
    progress_bar = utils.Progress_bar()
    num_frames_count = 0
    if not is_maintain_fixed_dt:
        with progress_bar:
            for i in progress_bar.track(range(data_size)):
                if i % plotting_freq != 0:
                    continue

                # Plot the trajectory from the beginning to current position
                if traj_len == -1:
                    start_index = 0
                else:
                    start_index = np.clip(i + 1 - traj_len, 0, None)

                for j in range(objects_count):
                    colors_j = None
                    labels_j = None
                    if colors is not None:
                        colors_j = colors[j]
                    if labels is not None:
                        labels_j = labels[j]
                    traj = ax.plot(
                        sol_state[start_index : (i + 1), j * 3],
                        sol_state[start_index : (i + 1), j * 3 + 1],
                        sol_state[start_index : (i + 1), j * 3 + 2],
                        color=colors_j,
                    )
                    # Plot the current position as a filled circle
                    ax.plot(
                        sol_state[i, j * 3],
                        sol_state[i, j * 3 + 1],
                        sol_state[i, j * 3 + 2],
                        color=traj[0].get_color(),
                        label=labels_j,
                        marker=marker,
                        markersize=markersize,
                    )
                # Add legend at the beginning
                # Warning: if we keep adding legend, the memory consumption will keep increasing
                if legend and i == 0:
                    fig.legend(loc="center right", borderaxespad=0.2)

                    # Adjust figure for the legend
                    if i == 0:
                        fig.subplots_adjust(right=0.7)
                        fig.tight_layout()

                # Set axis labels before capturing the frame
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                ax.set_zlabel(zlabel)  # type: ignore

                if not is_dynamic_axes:
                    ax.set_xlim3d((xlim_min, xlim_max))  # type: ignore
                    ax.set_ylim3d((ylim_min, ylim_max))  # type: ignore
                    ax.set_zlim3d((zlim_min, zlim_max))  # type: ignore
                else:
                    set_3d_axes_equal(ax)

                # Set equal aspect ratio to prevent distortion
                ax.set_aspect("equal")

                # Capture the frame
                plt.savefig(
                    Path(file_path).parent / f"frames_{num_frames_count:06d}.png",
                    dpi=dpi,
                )
                num_frames_count += 1

                # Clear the plot to prepare for the next frame
                ax.clear()

    else:
        if sol_time is None:
            raise ValueError("Solution time is required to maintain fixed dt")

        # Attempt to maintain fixed dt for the animation
        frame_size = int(data_size / plotting_freq) + 1
        plot_time = np.linspace(
            sol_time[0],
            sol_time[-1],
            frame_size,
        )
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
                    colors_j = None
                    labels_j = None
                    if colors is not None:
                        colors_j = colors[j]
                    if labels is not None:
                        labels_j = labels[j]
                    traj = ax.plot(
                        sol_state[(start_index + 1) : (index + 1), j * 3],
                        sol_state[(start_index + 1) : (index + 1), j * 3 + 1],
                        sol_state[(start_index + 1) : (index + 1), j * 3 + 2],
                        color=colors_j,
                    )
                    # Plot the current position as a filled circle
                    ax.plot(
                        sol_state[index, j * 3],
                        sol_state[index, j * 3 + 1],
                        sol_state[index, j * 3 + 2],
                        color=traj[0].get_color(),
                        label=labels_j,
                        marker=marker,
                        markersize=markersize,
                    )
                # Add legend at the beginning
                # Warning: if we keep adding legend, the memory consumption will keep increasing
                if legend and i == 0:
                    fig.legend(loc="center right", bbox_to_anchor=(1.325, 0.5))

                # Set axis labels before capturing the frame
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                ax.set_zlabel(zlabel)  # type: ignore

                if not is_dynamic_axes:
                    ax.set_xlim3d((xlim_min, xlim_max))  # type: ignore
                    ax.set_ylim3d((ylim_min, ylim_max))  # type: ignore
                    ax.set_zlim3d((zlim_min, zlim_max))  # type: ignore
                else:
                    set_3d_axes_equal(ax)

                # Set equal aspect ratio to prevent distortion
                ax.set_aspect("equal")

                # Capture the frame
                plt.savefig(
                    Path(file_path).parent / f"frames_{num_frames_count:06d}.png",
                    dpi=dpi,
                )
                num_frames_count += 1

                # Clear the plot to prepare for the next frame
                ax.clear()

    print("Combining frames to gif...")

    def frames_generator(num_frames_count):
        for i in range(num_frames_count):
            yield PIL.Image.open(file_path.parent / f"frames_{i:06d}.png")

    frames = frames_generator(num_frames_count)
    next(frames).save(
        file_path,
        format="GIF",
        append_images=frames,
        save_all=True,
        duration=(1000 / fps),
        loop=0,
    )

    plt.close("all")

    for i in range(num_frames_count):
        (Path(file_path).parent / f"frames_{i:06d}.png").unlink()


def plot_quantity_against_time(
    quantity: np.ndarray,
    sol_time: np.ndarray,
    is_log_y: bool = False,
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    save_fig: bool = False,
    save_fig_path: Optional[str | Path] = None,
) -> None:
    """Plot a quantity against time

    Parameters
    ----------
    quantity : np.ndarray
        Quantity to be plotted
    sol_time : np.ndarray
        Solution time
    is_log_y : bool, optional
        Flag to check whether y-axis is in log scale, by default False
    title : Optional[str], optional
        Title of the plot, by default None
    xlabel : Optional[str], optional
        Label of x-axis, by default None
    ylabel : Optional[str], optional
        Label of y-axis, by default None
    save_fig : bool, optional
        Flag to check whether to save the figure, by default False
    save_fig_path : Optional[str], optional
        Path to save the figure, by default None
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if is_log_y:
        ax.semilogy(sol_time, quantity)
    else:
        ax.plot(sol_time, quantity)

    if title is not None:
        ax.set_title(title)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)

    if save_fig:
        plt.savefig(save_fig_path, dpi=300)
    else:
        plt.show()
    plt.close("all")


def plot_eccentricity_or_inclination(
    eccentricity_or_inclination: np.ndarray,
    sol_time: np.ndarray,
    colors: Optional[list[str]] = None,
    labels: Optional[list[str]] = None,
    legend: bool = False,
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    save_fig: bool = False,
    save_fig_path: Optional[str | Path] = None,
) -> None:
    """Plot eccentricity or inclination against time

    Parameters
    ----------
    eccentricity_or_inclination : np.ndarray
        Eccentricity or inclination data
    sol_time : np.ndarray
        Solution time
    colors : Optional[list[str]], optional
        Colors of different data, by default None
    labels : Optional[list[str]], optional
        Labels of different data, by default None
    legend : bool, optional
        Flag to check whether to show legend, by default False
    title : Optional[str], optional
        Title of the plot, by default None
    xlabel : Optional[str], optional
        Label of x-axis, by default None
    ylabel : Optional[str], optional
        Label of y-axis, by default None
    save_fig : bool, optional
        Flag to check whether to save the figure, by default False
    save_fig_path : Optional[str], optional
        Path to save the figure, by default None
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(eccentricity_or_inclination.shape[1]):
        if colors is not None:
            colors_i = colors[i]
        else:
            colors_i = None

        if labels is not None:
            labels_i = labels[i]
        else:
            labels_i = None

        ax.plot(
            sol_time, eccentricity_or_inclination[:, i], color=colors_i, label=labels_i
        )

    if title is not None:
        ax.set_title(title)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)

    if legend:
        fig.legend(loc="center right", borderaxespad=0.2)
        fig.tight_layout()

    if save_fig:
        plt.savefig(save_fig_path, dpi=300)
    else:
        plt.show()
    plt.close("all")
