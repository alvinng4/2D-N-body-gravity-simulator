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
    AVAILABLE_INTEGRATORS_TO_PRINTABLE_NAMES = {
        "euler": "Euler",
        "euler_cromer": "Euler_Cromer",
        "rk4": "RK4",
        "leapfrog": "LeapFrog",
        "rkf45": "RKF45",
        "dopri": "DOPRI",
        "dverk": "DVERK",
        "rkf78": "RKF78",
        "ias15": "IAS15",
    }

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

    def sol_state_to_system(
        self, index: int = -1, system_name: str = None, objects_names: list = None
    ) -> GravitationalSystem:
        """
        Convert the latest state of the solution to a new GravitationalSystem object

        Parameters
        ----------
        index : int (optional)
            Index of the solution state. Default is the latest state.
        system_name : str (optional)
            Name of the system.
        objects_names : list (optional)
            List of names of the objects in the system.

        Returns
        -------
        GravitationalSystem object
        """
        system = GravitationalSystem()
        system.name = system_name

        if objects_names is not None:
            if len(objects_names) < self.simulator.objects_count:
                temp = [
                    None
                    for _ in range(self.simulator.objects_count - len(objects_names))
                ]
                objects_names += temp
            elif len(objects_names) > self.simulator.objects_count:
                warnings.warn(
                    "Number of names provided is greater than number of objects. Ignoring extra names."
                )

        for i in range(self.simulator.objects_count):
            x = self.simulator.sol_state[index, (i * 3) : (i * 3) + 3]
            v = self.simulator.sol_state[
                index,
                ((self.simulator.objects_count + i) * 3) : (
                    self.simulator.objects_count + i
                )
                * 3
                + 3,
            ]
            m = self.simulator.m[i]

            if objects_names is None:
                system.add(x, v, m)
            else:
                system.add(x, v, m, objects_name=objects_names[i])

        return system

    def save_results(
        self,
        system_name: str = None,
        path: str = None,
        computed_energy: bool = False,
        store_energy_as_zeros: bool = False,
    ) -> None:
        """
        Save the results in a csv file
        Unit: Solar masses, AU, day
        Format: time, G, dt, total energy, x1, y1, z1, x2, y2, z2, ... vx1, vy1, vz1, vx2, vy2, vz2, ...
        """
        data_size = len(self.simulator.sol_time)

        if not computed_energy or store_energy_as_zeros:
            print("Computing energy...")
            self.simulator.compute_energy()
            energy = self.simulator.energy
        elif store_energy_as_zeros:
            energy = np.zeros(data_size)

        # Storing the result
        print("Storing simulation results...")
        file_path = Path(__file__).parent / "results"
        file_path.mkdir(parents=True, exist_ok=True)
        file_path = (
            Path(__file__).parent
            / "results"
            / (
                str(datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
                + "_result.csv"
            )
        )

        # Storing metadata
        with open(file_path, "w", newline="") as file:
            writer = csv.writer(file, quoting=csv.QUOTE_NONE)
            writer.writerow(
                [
                    f"# Data saved on (YYYY-MM-DD): {str(datetime.datetime.now().strftime('%Y-%m-%d'))}"
                ]
            )
            writer.writerow([f"# System Name: {system_name}"])

            try:
                integrator_name = (
                    GravitySimulator.AVAILABLE_INTEGRATORS_TO_PRINTABLE_NAMES[
                        self.simulator.integrator
                    ]
                )
            except KeyError:
                integrator_name = None

            writer.writerow([f"# Integrator: {integrator_name}"])
            writer.writerow([f"# Number of objects: {self.simulator.objects_count}"])
            writer.writerow([f"# Simulation time (days): {self.simulator.tf}"])
            writer.writerow([f"# dt (days): {self.simulator.dt}"])
            writer.writerow([f"# Tolerance: {self.simulator.tolerance}"])
            writer.writerow([f"# Data size: {data_size}"])
            writer.writerow(
                [f"# Store every nth point: {self.simulator.store_every_n}"]
            )
            writer.writerow([f"# Run time (s): {self.simulator.run_time}"])
            masses_str = " ".join(map(str, self.simulator.m))
            writer.writerow([f"# masses: {masses_str}"])

        progress_bar = Progress_bar()
        with progress_bar:
            with open(file_path, "a", newline="") as file:
                writer = csv.writer(file)
                for count in progress_bar.track(range(data_size)):
                    row = np.insert(
                        self.simulator.sol_state[count],
                        0,
                        energy[count],
                    )
                    row = np.insert(row, 0, self.simulator.sol_dt[count])
                    row = np.insert(row, 0, self.simulator.sol_time[count])
                    writer.writerow(row.tolist())
        print(f"Storing completed. Please check {file_path}")
