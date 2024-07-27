"""
N-body gravity simulator API
"""

import ctypes
import datetime
import math
from pathlib import Path
import platform
import sys
import typing
import warnings

import numpy as np

sys.path.append(str(Path(__file__).parent))

import common
from gravitational_system import GravitationalSystem
from simulator import Simulator
import plotting


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

        self.is_exit = ctypes.c_bool(False)
        self.simulator = Simulator(c_lib=self.c_lib, is_exit_ctypes_bool=self.is_exit)
        self.compute_energy = self.simulator.compute_energy
        self.compute_angular_momentum = self.simulator.compute_angular_momentum

    def create_system(self, system_name: str = None):
        """
        Create a gravitational system

        Parameters
        ----------
        system_name : str (optional)
            Name of the system

        Returns
        -------
        GravitationalSystem object
        """
        return GravitationalSystem(name=system_name)

    def set_current_system(self, system: GravitationalSystem) -> None:
        """
        Set the current gravitational system

        Parameters
        ----------
        system : GravitationalSystem object
            GravitationalSystem object
        """
        self.current_system = system

    def days_to_years(
        self, days: typing.Union[float, np.ndarray]
    ) -> typing.Union[float, np.ndarray]:
        return days / self.DAYS_PER_YEAR

    def years_to_days(
        self, years: typing.Union[float, np.ndarray]
    ) -> typing.Union[float, np.ndarray]:
        return years * self.DAYS_PER_YEAR

    def launch_simulation(
        self,
        integrator: str,
        tf: float,
        dt: float = None,
        tolerance: float = None,
        store_every_n: int = 1,
        acceleration: str = "pairwise",
        flush: bool = False,
        flush_results_path: str = None,
        no_progress_bar: bool = False,
        no_print: bool = False,
    ) -> None:
        try:
            self.simulator.launch_simulation(
                self.current_system,
                integrator,
                tf,
                dt,
                tolerance,
                store_every_n,
                acceleration,
                flush,
                flush_results_path,
                no_progress_bar,
                no_print,
            )
        except KeyboardInterrupt:
            if not no_print:
                print("KeyboardInterrupt detected. Exiting simulation...")
            self.is_exit.value = True  # Exit simulation in c_lib
            sys.exit(0)

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

    def plot_rel_energy(
        self,
        energy: np.ndarray = None,
        sol_time: np.ndarray = None,
        title="Relative energy error against time",
        xlabel=f"Time",
        ylabel="$|(E(t)-E_0)/E_0|$",
    ):
        if energy is None:
            try:
                energy = self.simulator.energy
            except AttributeError:
                self.compute_energy()
                energy = self.simulator.energy
        if sol_time is None:
            sol_time = self.simulator.sol_time

        plotting.plot_rel_energy(energy, sol_time, title, xlabel, ylabel)

    def plot_rel_angular_momentum(
        self,
        angular_momentum: np.ndarray = None,
        sol_time: np.ndarray = None,
        title: str = "Relative angular momentum error against time",
        xlabel: str = "Time",
        ylabel: str = "$|(L(t)-L_0)/L_0|$",
    ):
        if angular_momentum is None:
            try:
                angular_momentum = self.simulator.angular_momentum
            except AttributeError:
                self.compute_angular_momentum()
                angular_momentum = self.simulator.angular_momentum
        if sol_time is None:
            sol_time = self.simulator.sol_time

        plotting.plot_rel_angular_momentum(
            angular_momentum, sol_time, title, xlabel, ylabel
        )

    def plot_dt(
        self,
        sol_dt: np.ndarray = None,
        sol_time: np.ndarray = None,
        title: str = "dt against time",
        xlabel: str = "Time",
        ylabel: str = "$dt$ (days)",
        yscale: str = "log",
        marker_size: float = 0.1,
    ) -> None:
        if sol_dt is None:
            sol_dt = self.simulator.sol_dt
        if sol_time is None:
            sol_time = self.simulator.sol_time

        plotting.plot_dt(sol_dt, sol_time, title, xlabel, ylabel, yscale, marker_size)

    def sol_state_to_system(self, index: int = -1) -> GravitationalSystem:
        """
        Convert the latest state of the solution to a new GravitationalSystem object

        Parameters
        ----------
        index : int (optional)
            Index of the solution state. Default is the latest state.

        Returns
        -------
        GravitationalSystem object
        """
        new_system = GravitationalSystem()
        new_system.name = self.current_system.name

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
            new_system.add(x, v, m, object_name=self.current_system.objects_names[i])

        return new_system

    def simulator_to_system(self) -> GravitationalSystem:
        """
        Convert the current state of the simulator to a new GravitationalSystem object

        Returns
        -------
        GravitationalSystem object
        """
        new_system = GravitationalSystem()
        new_system.name = self.current_system.name

        for i in range(self.simulator.objects_count):
            new_system.add(
                self.simulator.x,
                self.simulator.v,
                self.simulator.m,
                object_name=self.current_system.objects_names[i],
            )

        return new_system

    def save_results(
        self,
        file_path: str = None,
        computed_energy: bool = False,
        store_energy_as_zeros: bool = False,
        no_progress_bar: bool = False,
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

        try:
            integrator_name = GravitySimulator.AVAILABLE_INTEGRATORS_TO_PRINTABLE_NAMES[
                self.simulator.integrator
            ]
        except KeyError:
            integrator_name = None

        common.save_results(
            file_path=file_path,
            system_name=self.current_system.name,
            integrator_name=integrator_name,
            objects_count=self.simulator.objects_count,
            G=self.simulator.G,
            tf=self.simulator.tf,
            dt=self.simulator.dt,
            tolerance=self.simulator.tolerance,
            data_size=data_size,
            store_every_n=self.simulator.store_every_n,
            run_time=self.simulator.run_time,
            masses=self.simulator.m,
            sol_state=self.simulator.sol_state,
            sol_time=self.simulator.sol_time,
            sol_dt=self.simulator.sol_dt,
            energy=energy,
            no_progress_bar=no_progress_bar,
        )

    def read_results(
        self,
        file_path: str,
        start: int = 0,
        end: int = -1,
        step: int = 1,
        memory_buffer_size: int = 50000,
        no_print: bool = False,
        no_progress_bar: bool = False,
    ) -> None:
        """
        Read the results from a csv file,
        and load the data to the simulator

        Parameters
        ----------
        file_path : str
            Path to the csv file
        start : int (optional)
            Start index of the data to read
        end : int (optional)
            End index of the data to read
        step : int (optional)
            Step size to read the data
        memory_buffer_size : int (optional)
            Memory buffer size for storing data
        no_print : bool (optional)
            Disable print statements
        no_progress_bar : bool (optional)
            Disable progress bar
        """

        (
            self.simulator.sol_state,
            self.simulator.sol_time,
            self.simulator.sol_dt,
            self.simulator.energy,
            self.simulator.system_name,
            self.simulator.integrator,
            self.simulator.objects_count,
            self.simulator.G,
            self.simulator.tf,
            self.simulator.dt,
            self.simulator.tolerance,
            self.store_every_n,
            self.simulator.run_time,
            self.simulator.m,
            self.data_size,
        ) = common.read_results(
            file_path=file_path,
            start=start,
            end=end,
            step=step,
            memory_buffer_size=memory_buffer_size,
            no_print=no_print,
            no_progress_bar=no_progress_bar,
        )

    @property
    def integration_mode(self) -> str:
        return self.simulator.integration_mode

    @integration_mode.setter
    def integration_mode(self, value: str) -> None:
        match value:
            case "c_lib":
                if self.c_lib is not None:
                    self.simulator.integration_mode = value
                else:
                    raise ValueError("C library is not available.")

            case "numpy":
                self.simulator.integration_mode = value

            case _:
                raise ValueError('integration_mode must be either "c_lib" or "numpy".')

    @property
    def objects_count(self) -> int:
        return self.simulator.objects_count

    @objects_count.setter
    def objects_count(self, value: int) -> None:
        self.simulator.objects_count = value

    @property
    def sol_time(self) -> np.ndarray:
        return self.simulator.sol_time

    @sol_time.setter
    def sol_time(self, value: np.ndarray) -> None:
        self.simulator.sol_time = value

    @property
    def sol_state(self) -> np.ndarray:
        return self.simulator.sol_state

    @sol_state.setter
    def sol_state(self, value: np.ndarray) -> None:
        self.simulator.sol_state = value

    @property
    def sol_dt(self) -> np.ndarray:
        return self.simulator.sol_dt

    @sol_dt.setter
    def sol_dt(self, value: np.ndarray) -> None:
        self.simulator.sol_dt = value

    @property
    def data_size(self) -> int:
        return self.simulator.data_size

    @data_size.setter
    def data_size(self, value: int) -> None:
        self.simulator.data_size = value

    @property
    def energy(self) -> np.ndarray:
        return self.simulator.energy

    @energy.setter
    def energy(self, value: np.ndarray) -> None:
        self.simulator.energy = value

    @property
    def angular_momentum(self) -> np.ndarray:
        return self.simulator.angular_momentum

    @angular_momentum.setter
    def angular_momentum(self, value: np.ndarray) -> None:
        self.simulator.angular_momentum = value
