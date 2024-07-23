"""
API for N-body gravity simulator
"""

import ctypes
import csv
import datetime
import math
from pathlib import Path
import platform
import sys
import timeit
import typing
import warnings

import numpy as np

sys.path.append(str(Path(__file__).parent))

from common import get_bool
from common import get_int
from common import get_float
from gravitational_system import GravitationalSystem
from simulator import Simulator
import plotting
from progress_bar import Progress_bar


class GravitySimulator:
    DEFAULT_SYSTEMS = GravitationalSystem.DEFAULT_SYSTEMS
    AVAILABLE_INTEGRATORS = Simulator.AVAILABLE_INTEGRATORS
    DAYS_PER_YEAR = 365.242189

    def __init__(self):
        self.c_lib = None
        try:
            if platform.system() == "Windows":
                self.c_lib = ctypes.cdll.LoadLibrary(
                    str(Path(__file__).parent / "c_lib.dll")
                )
            elif platform.system() == "Darwin":
                self.c_lib = ctypes.cdll.LoadLibrary(
                    str(Path(__file__).parent / "c_lib.dylib")
                )
            elif platform.system() == "Linux":
                self.c_lib = ctypes.cdll.LoadLibrary(
                    str(Path(__file__).parent / "c_lib.so")
                )
        except Exception as e:
            warnings.warn(
                f"Failed to load C library (Exception: {e}).\n"
                + "NumPy is still available.\n"
                + "To use C library, try to compile the C files provided in the src folders."
            )

        self.simulator = Simulator(self.c_lib)

    def create_system(self, system: str = None):
        """
        Create a gravitational system

        Parameters
        ----------
        system : str (optional)
            Name of the system

        Returns
        -------
        GravitationalSystem object
        """
        return GravitationalSystem()

    def days_to_years(
        self, days: typing.Union[float, np.ndarray]
    ) -> typing.Union[float, np.ndarray]:
        return days / self.DAYS_PER_YEAR

    def years_to_days(
        self, years: typing.Union[float, np.ndarray]
    ) -> typing.Union[float, np.ndarray]:
        return years * self.DAYS_PER_YEAR

    def plot_2d_trajectory(
        self,
        colors=None,
        labels=None,
        legend=False,
        xlabel="$x$ (AU)",
        ylabel="$y$ (AU)",
        marker="o",
        markersize=6,
    ):
        plotting.plot_2d_trajectory(
            self.simulator.objects_count,
            self.simulator.sol_state,
            colors=colors,
            labels=labels,
            legend=legend,
            xlabel=xlabel,
            ylabel=ylabel,
            marker=marker,
            markersize=markersize,
        )

    def plot_3d_trajectory(
        self,
        colors=None,
        labels=None,
        legend=False,
        xlabel="$x$ (AU)",
        ylabel="$y$ (AU)",
        zlabel="$z$ (AU)",
        marker="o",
        markersize=6,
    ):
        plotting.plot_3d_trajectory(
            self.simulator.objects_count,
            self.simulator.sol_state,
            colors=colors,
            labels=labels,
            legend=legend,
            xlabel=xlabel,
            ylabel=ylabel,
            zlabel=zlabel,
            marker=marker,
            markersize=markersize,
        )

    def animate_2d_traj_gif(
        self,
        fps: int = 30,
        animation_length: float = None,
        plot_every_nth_point: int = None,
        dpi=200,
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
        if animation_length is not None and animation_length <= 0.0:
            raise ValueError("gif_time must be greater than 0.")

        if animation_length is None and plot_every_nth_point is None:
            animation_length = 10.0
            plot_every_nth_point = math.floor(
                len(self.simulator.sol_state) / (animation_length * fps)
            )
            if plot_every_nth_point <= 0:
                plot_every_nth_point = 1

        elif animation_length is None and plot_every_nth_point is not None:
            if plot_every_nth_point <= 0:
                warnings.warn(
                    "plot_every_nth_point must be greater than 0. Using default value of 1."
                )
                plot_every_nth_point = 1
            elif plot_every_nth_point > len(self.simulator.sol_state):
                warnings.warn(
                    f"plot_every_nth_point is greater than data size ({len(self.simulator.sol_state)}). Using default value of 1."
                )
                plot_every_nth_point = 1

        elif animation_length is not None and plot_every_nth_point is None:
            plot_every_nth_point = math.floor(
                len(self.simulator.sol_state) / (animation_length * fps)
            )
            if plot_every_nth_point <= 0:
                plot_every_nth_point = 1

        elif animation_length is not None and plot_every_nth_point is not None:
            warnings.warn(
                "Both animation_length and plot_every_nth_point are provided. plot_every_nth_point will be ignored."
            )
            plot_every_nth_point = math.floor(
                len(self.simulator.sol_state) / (animation_length * fps)
            )
            if plot_every_nth_point <= 0:
                plot_every_nth_point = 1

        plotting.animate_2d_traj_gif(
            self.simulator.objects_count,
            self.simulator.sol_state,
            fps=fps,
            plot_every_nth_point=plot_every_nth_point,
            dpi=dpi,
            is_dynamic_axes=is_dynamic_axes,
            axes_lim=axes_lim,
            is_maintain_fixed_dt=is_maintain_fixed_dt,
            traj_len=traj_len,
            colors=colors,
            labels=labels,
            legend=legend,
            xlabel=xlabel,
            ylabel=ylabel,
            marker=marker,
            markersize=markersize,
            file_name=file_name,
            file_path=file_path,
            sol_time=self.simulator.sol_time,
        )

    def animate_3d_traj_gif(
        self,
        fps: int = 30,
        animation_length: float = None,
        plot_every_nth_point: int = None,
        dpi=200,
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
        if animation_length is None and plot_every_nth_point is None:
            animation_length = 10.0
            plot_every_nth_point = math.floor(
                len(self.simulator.sol_state) / (animation_length * fps)
            )
            if plot_every_nth_point <= 0:
                plot_every_nth_point = 1

        elif animation_length is None and plot_every_nth_point is not None:
            if plot_every_nth_point <= 0:
                warnings.warn(
                    "plot_every_nth_point must be greater than 0. Using default value of 1."
                )
                plot_every_nth_point = 1
            elif plot_every_nth_point > len(self.simulator.sol_state):
                warnings.warn(
                    f"plot_every_nth_point is greater than data size ({len(self.simulator.sol_state)}). Using default value of 1."
                )
                plot_every_nth_point = 1

        elif animation_length is not None and plot_every_nth_point is None:
            plot_every_nth_point = math.floor(
                len(self.simulator.sol_state) / (animation_length * fps)
            )
            if plot_every_nth_point <= 0:
                plot_every_nth_point = 1

        elif animation_length is not None and plot_every_nth_point is not None:
            warnings.warn(
                "Both animation_length and plot_every_nth_point are provided. plot_every_nth_point will be ignored."
            )
            plot_every_nth_point = math.floor(
                len(self.simulator.sol_state) / (animation_length * fps)
            )
            if plot_every_nth_point <= 0:
                plot_every_nth_point = 1

        plotting.animate_3d_traj_gif(
            self.simulator.objects_count,
            self.simulator.sol_state,
            fps=fps,
            plot_every_nth_point=plot_every_nth_point,
            dpi=dpi,
            is_dynamic_axes=is_dynamic_axes,
            axes_lim=axes_lim,
            is_maintain_fixed_dt=is_maintain_fixed_dt,
            traj_len=traj_len,
            colors=colors,
            labels=labels,
            legend=legend,
            xlabel=xlabel,
            ylabel=ylabel,
            zlabel=zlabel,
            marker=marker,
            markersize=markersize,
            file_name=file_name,
            file_path=file_path,
            sol_time=self.simulator.sol_time,
        )
