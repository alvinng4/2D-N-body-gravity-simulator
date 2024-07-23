"""
API for N-body gravity simulator
"""

import argparse
import ctypes
import csv
import datetime
import math
from pathlib import Path
import platform
import re
import sys
import timeit
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
            self.has_c_lib = True
        except Exception as e:
            warnings.warn(
                f"Failed to load C library (Exception: {e}).\n"
                + "NumPy is still available.\n"
                + "To use C library, try to compile the C files provided in the src folders."
            )
            self.has_c_lib = False

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
        return GravitationalSystem(system)

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
        fps,
        plot_every_nth_point,
        dpi,
        is_dynamic_axes,
        is_custom_axes,
        axes_lim,
        is_maintain_fixed_dt,
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
        plotting.animate_2d_traj_gif(
            self.simulator.objects_count,
            self.simulator.sol_state,
            fps,
            plot_every_nth_point,
            dpi,
            is_dynamic_axes,
            is_custom_axes,
            axes_lim,
            is_maintain_fixed_dt,
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
        fps: int,
        plot_every_nth_point,
        dpi,
        is_dynamic_axes,
        is_custom_axes,
        axes_lim,
        is_maintain_fixed_dt,
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
        plotting.animate_3d_traj_gif(
            self.simulator.objects_count,
            self.simulator.sol_state,
            fps,
            plot_every_nth_point,
            dpi,
            is_dynamic_axes,
            is_custom_axes,
            axes_lim,
            is_maintain_fixed_dt,
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
