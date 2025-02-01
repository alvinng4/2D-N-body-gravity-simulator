"""
Command line interface (CLI) for gravity simulator

Usage:
    python3 -m gravity_sim [options]

Author: Ching Yin Ng
"""

import argparse
import ctypes
import datetime
import re
import sys
from pathlib import Path
from typing import Any, Optional

import numpy as np

from .gravitational_system import GravitationalSystem
from .simulator import Simulator
from . import plotting
from . import utils


class GravitySimulatorCLI:
    C_LIB_BUILT_IN_SYSTEMS = [
        "circular_binary_orbit",
        "eccentric_binary_orbit",
        "3d_helix",
        "sun_earth_moon",
        "figure-8",
        "pyth-3-body",
        "solar_system",
        "solar_system_plus",
    ]
    SOLAR_LIKE_SYSTEM_OBJECTS_PAIRS = {
        "sun_earth_moon": ["Sun", "Earth", "Moon"],
        "solar_system": [
            "Sun",
            "Mercury",
            "Venus",
            "Earth",
            "Mars",
            "Jupiter",
            "Saturn",
            "Uranus",
            "Neptune",
        ],
        "solar_system_plus": [
            "Sun",
            "Mercury",
            "Venus",
            "Earth",
            "Mars",
            "Jupiter",
            "Saturn",
            "Uranus",
            "Neptune",
            "Pluto",
            "Ceres",
            "Vesta",
        ],
    }

    def __init__(self) -> None:
        ### Read command line arguments ###
        args = self._read_command_line_arg()
        self.is_save_plots = args.save_plots
        if args.c_lib_path is not None:
            self.c_lib_path = Path(args.c_lib_path)
        self.is_exit_ctypes_bool = ctypes.c_bool(False)
        self.default_results_folder = Path(__file__).parent / "results"
        self.default_results_folder.mkdir(exist_ok=True, parents=True)

        self.has_data_attr: list[str] = []

    def run_prog(self) -> None:
        try:
            print("--------------------------------------------------")
            print("Gravity Simulator - Command Line Interface")
            print("Exit program anytime by hitting Ctrl + C")

            print()
            print("Applied options:")
            applied_option = False
            if self.is_save_plots:
                print("Save plots: True")
                applied_option = True
            if hasattr(self, "c_lib_path"):
                print(f"C library path: {self.c_lib_path.absolute()}")
                applied_option = True

            if not applied_option:
                print("No options applied.")

            print()

            ### Initialization ###
            print("Initializing...")
            # Load C library
            print("Loading C library...", end="")
            if hasattr(self, "c_lib_path"):
                self.c_lib = self._load_c_lib(self.c_lib_path)
            else:
                self.c_lib = self._load_c_lib()

            if self.c_lib is None:
                sys.exit("C library not loaded. Exiting program...")
            utils.initialize_c_lib(self.c_lib)

            # Initialize simulator
            self.simulator = Simulator(self.c_lib)
            print("Done.")
            print("--------------------------------------------------")

            ### Start main loop ###
            self.main_loop()

        except KeyboardInterrupt:
            self.is_exit_ctypes_bool.value = True
            print("KeyboardInterrupt detected. Exiting program...")
            sys.exit(1)

    def main_loop(self) -> None:
        while True:
            action = self._user_interface_before_simulation()

            if action == 1:
                self._delete_previous_simulation_data()
                user_cancelled_simulation = self._launch_simulation()
                if user_cancelled_simulation:
                    continue
                is_exit_program = self._user_interface_after_simulation()
                if is_exit_program:
                    break
            elif action == 2:
                self._delete_previous_simulation_data()
                self._read_results()
                is_exit_program = self._user_interface_after_simulation()
                if is_exit_program:
                    break
            elif action == 3:
                print("Exiting the program...")
                break
            else:
                raise ValueError("Invalid action")

    @staticmethod
    def _read_command_line_arg() -> argparse.Namespace:
        parser = argparse.ArgumentParser(description="Gravity Simulator")
        parser.add_argument(
            "--save_plots",
            "-s",
            action="store_true",
            help="Save plots directly to results folder",
        )
        parser.add_argument(
            "--c_lib_path",
            "-c",
            help="C library path",
        )
        args = parser.parse_args()

        return args

    @staticmethod
    def _load_c_lib(c_lib_path: Optional[Path] = None) -> Optional[ctypes.CDLL]:
        c_lib = None
        try:
            if c_lib_path is None:
                c_lib = utils.load_c_lib()
            else:
                c_lib = utils.load_c_lib(c_lib_path)
        except FileNotFoundError as err_msg:
            print(f"FileNotFoundError: {err_msg}")
            print(
                "To fix this error, please provide the correct path to the C library, "
                + "with the --c_lib_path option. "
                + 'E.g. python3 gravity_sim --c_lib_path "/path/to/c_lib.so".'
            )
        except OSError as err_msg:
            print(f"OSError: {err_msg}")

        if c_lib is None:
            raise ValueError("C library not loaded")

        return c_lib

    @staticmethod
    def _user_interface_before_simulation() -> int:
        msg = (
            "Select an action:\n"
            + "1. Launch simulation\n"
            + "2. Read simulation data\n"
            + "3. Exit\n"
            + "Enter action (Number): "
        )
        return_value = GravitySimulatorCLI.get_int(msg, larger_than=0, smaller_than=4)
        print()
        return return_value

    def _user_interface_after_simulation(self) -> bool:
        msg = (
            "Select an action:\n"
            + "1. Plot 2D trajectory (xy plane)\n"
            + "2. Plot 3D trajectory\n"
            + "3. Animate 2D trajectory (gif)\n"
            + "4. Animate 3D trajectory (gif)\n"
            + "5. Plot relative energy error\n"
            + "6. Plot relative angular momentum error\n"
            + "7. Plot dt\n"
            + "8. Plot eccentricity\n"
            + "9. Plot inclination\n"
            + "10. Read data size\n"
            + "11. Trim data\n"
            + "12. Save simulation data\n"
            + "13. Restart program\n"
            + "14. Exit\n"
            + "Enter action (Number): "
        )

        while True:
            action = self.get_int(msg, larger_than=0, smaller_than=15)
            print()

            match action:
                case 1:
                    self._plot_2d_trajectory_wrapper()
                case 2:
                    self._plot_3d_trajectory_wrapper()
                case 3:
                    self._animate_2d_trajectory_wrapper()
                case 4:
                    self._animate_3d_trajectory_wrapper()
                case 5:
                    self._plot_rel_energy_error_wrapper()
                case 6:
                    self._plot_rel_angular_momentum_error_wrapper()
                case 7:
                    self._plot_dt_wrapper()
                case 8:
                    self._plot_eccentricity_wrapper()
                case 9:
                    self._plot_inclination_wrapper()
                case 10:
                    print(f"There are {self.data_size_} lines of data.")
                    print()
                case 11:
                    print(f"There are {self.data_size_} lines of data.")
                    self._trim_data()
                case 12:
                    self._save_results()
                case 13:
                    return False
                case 14:
                    print("Exiting the program...")
                    return True

    def _launch_simulation(self) -> bool:
        """Launch simulation

        Returns
        -------
        user_cancelled_simulation : bool
            True if user cancelled the simulation, False otherwise
        """
        while True:
            (
                gravitational_system,
                integrator_params,
                acceleration_params,
                storing_params,
                settings,
                tf,
                system_name,
            ) = self._get_user_simulation_input()
            self._print_user_simulation_input(
                integrator_params, storing_params, tf, system_name
            )
            if self.get_bool("Proceed? (y/n): "):
                print()
                break
            else:
                if self.get_bool("Back to main menu? (y/n): "):
                    print()
                    return True

        gravitational_system.load(system_name)
        self.simulator.launch_simulation(
            gravitational_system,
            integrator_params,
            acceleration_params,
            storing_params,
            settings,
            tf,
            self.is_exit_ctypes_bool,
        )
        self.gravitational_system = self.simulator.gravitational_system
        self.gravitational_system.name = system_name
        self.sol_state_ = self.simulator.sol_state_
        self.sol_time_ = self.simulator.sol_time_
        self.sol_dt_ = self.simulator.sol_dt_
        self.data_size_ = self.simulator.data_size_

        self.has_data_attr.append("gravitational_system")
        self.has_data_attr.append("sol_state_")
        self.has_data_attr.append("sol_time_")
        self.has_data_attr.append("sol_dt_")
        self.has_data_attr.append("data_size_")

        return False

    def _delete_previous_simulation_data(self) -> None:
        attr_list_copy = self.has_data_attr.copy()
        for attr in attr_list_copy:
            if hasattr(self, attr):
                delattr(self, attr)
            self.has_data_attr.remove(attr)

    def _get_user_simulation_input(
        self,
    ) -> tuple[GravitationalSystem, dict, dict, dict, dict, float, str]:
        """Get user input for simulation

        Returns
        -------
        gravitational_system : GravitationalSystem
        integrator_params : dict
        acceleration_params : dict
        storing_params : dict
        settings : dict
        tf : float
        system_name : str
        """
        print("------------------------------------")
        gravitational_system = GravitationalSystem()
        integrator_params: dict[str, str | float | bool] = {
            "initial_dt": 0.0,
            "whfast_kepler_tol": 1e-12,
            "whfast_kepler_max_iter": 500,
            "whfast_kepler_auto_remove": False,
            "whfast_auto_remove_tol": 1e-8,
        }
        acceleration_params: dict[str, str | int | float] = {
            "method": "pairwise",
            "opening_angle": 0.0,
            "softening_length": 0.0,
            "order": 0,
        }
        storing_params: dict[str, str | int] = {
            "method": "default",
        }
        settings: dict[str, bool | int] = {
            "make_copy_params": False,
            "make_copy_system": False,
            "verbose": 2,
            "disable_progress_bar": False,
        }

        ### Select system ###
        msg = "Select a system:\n"
        for i, system in enumerate(self.C_LIB_BUILT_IN_SYSTEMS):
            msg += f"{i + 1}. {system}\n"
        msg += "Enter system number (int): "
        system_number = self.get_int(
            msg, larger_than=0, smaller_than=len(self.C_LIB_BUILT_IN_SYSTEMS) + 1
        )
        system_name = self.C_LIB_BUILT_IN_SYSTEMS[system_number - 1]
        print()

        ### Recommended settings ###
        try_recommended_settings = self.get_bool("Try recommended settings? (y/n): ")
        if try_recommended_settings:
            recommended_settings: list = (
                self.simulator.RECOMMENDED_SETTINGS_BUILT_IN_SYSTEMS[system_name]
            )
            integrator_params["integrator"] = "ias15"
            tf = recommended_settings[0]
            self.tf_units = recommended_settings[1]
            if self.tf_units == "years":
                tf *= self.simulator.DAYS_PER_YEAR
            integrator_params["dt"] = 0.0
            self.dt_units = "days"
            integrator_params["tolerance"] = recommended_settings[2]
            storing_params["storing_freq"] = recommended_settings[3]

        if not try_recommended_settings:
            print()
            ### Integrator parameters ###
            msg = "Choose an integrator:\n"
            for i, integrator in enumerate(self.simulator.AVAILABLE_INTEGRATORS):
                msg += f"{i + 1}. {integrator}\n"
            msg += "Enter integrator number (int): "
            integrator_number = self.get_int(
                msg,
                larger_than=0,
                smaller_than=len(self.simulator.AVAILABLE_INTEGRATORS) + 1,
            )
            integrator_params["integrator"] = self.simulator.AVAILABLE_INTEGRATORS[
                integrator_number - 1
            ]
            print()

            ### Input tf ###
            # tf = 0 is allowed as user may want to plot the
            # initial position of the system
            while True:
                user_input_tf = input("Enter tf (days/year) (e.g. 200y or 100d): ")
                if matches := re.search(
                    r"([0-9]*\.?[0-9]*)(?:\.|\W*)*(day|year|d|y)?",
                    user_input_tf,
                    re.IGNORECASE,
                ):
                    if not matches.group(1):
                        print("Invalid input. Please try again.")
                        print()
                        continue

                    try:
                        tf = float(matches.group(1))
                        if tf < 0:
                            raise ValueError
                    except ValueError:
                        print("Invalid input. Please try again.")
                        print()
                        continue

                    break
            if matches.group(2) in ["year", "y"]:
                self.tf_units = "years"
                tf *= self.simulator.DAYS_PER_YEAR
            else:
                self.tf_units = "days"

            self.has_data_attr.append("tf_units")
            print()

            ### Input dt / tolerance ###
            if (
                integrator_params["integrator"]
                in self.simulator.FIXED_STEP_SIZE_INTEGRATORS
            ):
                while True:
                    user_input_dt = input("Enter dt (days/year) (e.g. 1d): ")
                    if matches := re.search(
                        r"([0-9]*\.?[0-9]*)(?:\.|\W*)*(day|year|d|y)?\s*",
                        user_input_dt,
                        re.IGNORECASE,
                    ):
                        if not matches.group(1):
                            print("Invalid input. Please try again.")
                            print()
                            continue

                        try:
                            dt = float(matches.group(1))
                            if dt <= 0:
                                raise ValueError
                        except ValueError:
                            print("Invalid input. Please try again.")
                            print()
                            continue

                        break

                if matches.group(2) in ["year", "y"]:
                    self.dt_units = "years"
                    dt *= self.simulator.DAYS_PER_YEAR
                else:
                    self.dt_units = "days"

                self.has_data_attr.append("dt_units")

                integrator_params["dt"] = dt
                integrator_params["tolerance"] = 0.0

            elif (
                integrator_params["integrator"]
                in self.simulator.ADAPTIVE_STEP_SIZE_INTEGRATORS
            ):
                tolerance = self.get_float("Enter tolerance: ", larger_than=0.0)
                integrator_params["dt"] = 0.0
                self.dt_units = "days"
                integrator_params["tolerance"] = tolerance

            else:
                raise ValueError("Invalid integrator")

            print()

            ### Get storing frequency ###
            storing_freq = self.get_int(
                "Enter storing frequency (int): ", larger_than=0
            )
            storing_params["storing_freq"] = storing_freq

        return (
            gravitational_system,
            integrator_params,
            acceleration_params,
            storing_params,
            settings,
            tf,
            system_name,
        )

    def _print_user_simulation_input(
        self,
        integrator_params: dict,
        storing_params: dict,
        tf: float,
        system_name: str,
    ) -> None:
        print("------------------------------------")
        print("Input:")
        print(f"System: {system_name}")
        print(f"Integrator: {integrator_params['integrator']}")
        if self.tf_units == "years":
            print(f"tf: {tf / self.simulator.DAYS_PER_YEAR} years")
        else:
            print(f"tf: {tf} days")
        if integrator_params["dt"] != 0.0:
            if self.dt_units == "years":
                print(
                    f"dt: {integrator_params['dt'] / self.simulator.DAYS_PER_YEAR} years"
                )
            else:
                print(f"dt: {integrator_params['dt']} days")
        else:
            print(f"tolerance: {integrator_params['tolerance']}")
        print(f"Storing frequency: {storing_params['storing_freq']}")

        ### Estimate number of data points for fixed step size integrators ###
        if integrator_params["dt"] != 0.0:
            n_steps = int(tf / integrator_params["dt"])
            sol_size = int(n_steps / storing_params["storing_freq"] + 1)
            print(f"Estimated number of data points: {sol_size}")

        print()

    @staticmethod
    def get_bool(msg: str) -> bool:
        """Prompt user for boolean input

        Parameters
        ----------
        msg : str
            Message to display to user

        Returns
        -------
        user_input : bool

        Notes
        -----
        This function has following side effects:
        - Print msg
        - Print "Invalid input. Please try again." if user enters an invalid input
        """
        while True:
            if matches := re.search(r"^\s*(yes|no|y|n)\s*$", input(msg), re.IGNORECASE):
                if matches.group(1).lower() in ["y", "yes"]:
                    return True
                elif matches.group(1).lower() in ["n", "no"]:
                    return False

            print("Invalid input. Please try again.\n")

    @staticmethod
    def get_int(
        msg: str,
        larger_than: Optional[int] = None,
        smaller_than: Optional[int] = None,
        allow_cancel: Optional[bool] = False,
    ) -> int:
        """Prompt user for integer input

        Parameters
        ----------
        msg : str
            Message to display to user
        larger_than : Optional[int], optional
            Value must be larger than this, by default None
        smaller_than : Optional[int], optional
            Value must be smaller than this, by default None
        allow_cancel : Optional[bool], optional
            Allow user to enter "cancel" to cancel input, by default False

        Returns
        -------
        user_input : int

        Raises
        ------
        ValueError
            If user enters "cancel" and allow_cancel is enabled

        Notes
        -----
        This function has following side effects:
        - Print msg
        - Print "Value too small! Please try again." if user enters a value smaller than larger_than
        - Print "Value too big! Please try again." if user enters a value larger than smaller_than
        """
        while True:
            try:
                user_input = input(msg)
                if allow_cancel and user_input.strip().lower() == "cancel":
                    break

                user_input_int = int(user_input)
                if larger_than is not None and user_input_int <= larger_than:
                    print("Value too small! Please try again.")
                    print()
                    continue

                if smaller_than is not None and user_input_int >= smaller_than:
                    print("Value too big! Please try again.")
                    print()
                    continue

                return user_input_int

            except ValueError:
                print("Invalid input. Please try again.")
                print()

        raise ValueError("User canceled input")

    @staticmethod
    def get_float(
        msg: str,
        larger_than: Optional[float] = None,
        smaller_than: Optional[float] = None,
        allow_cancel: Optional[bool] = False,
    ) -> float:
        """Prompt user for float input

        Parameters
        ----------
        msg : str
            Message to display to user
        larger_than : Optional[float], optional
            Value must be larger than this, by default None
        smaller_than : Optional[float], optional
            Value must be smaller than this, by default None
        allow_cancel : Optional[bool], optional
            Allow user to enter "cancel" to cancel input, by default False

        Returns
        -------
        user_input : float

        Raises
        ------
        ValueError
            If user enters "cancel" and allow_cancel is enabled

        Notes
        -----
        This function has following side effects:
        - Print msg
        - Print "Value too small! Please try again." if user enters a value smaller than larger_than
        - Print "Value too big! Please try again." if user enters a value larger than smaller_than
        - Print "Invalid input. Please try again." if user enters an invalid input
        """
        while True:
            try:
                value = input(msg)
                if allow_cancel and value.strip().lower() == "cancel":
                    break

                user_input = float(value)
                if larger_than is not None and user_input <= larger_than:
                    print("Value too small! Please try again.")
                    print()
                    continue

                if smaller_than is not None and user_input >= smaller_than:
                    print("Value too big! Please try again.")
                    print()
                    continue

                return user_input

            except ValueError:
                print("Invalid input. Please try again.")
                print()

        raise ValueError("User canceled input")

    def _plot_2d_trajectory_wrapper(self) -> None:
        # Get kwargs
        kwargs: dict[str, Any] = {}
        kwargs["sol_state"] = self.sol_state_
        if self.is_save_plots:
            i = 0
            while True:
                save_fig_path = (
                    self.default_results_folder / f"2d_trajectory_{i:05d}.pdf"
                )
                if save_fig_path.exists():
                    i += 1
                else:
                    break

            kwargs["save_fig"] = True
            kwargs["save_fig_path"] = save_fig_path
            print(f"Saving plot to {save_fig_path}...", end="")
        else:
            print("Plotting 2D trajectory (xy plane)...(Please check the window)")

        if "gravitational_system" in self.has_data_attr:
            if self.gravitational_system.name in self.SOLAR_LIKE_SYSTEM_OBJECTS_PAIRS:
                labels = self.SOLAR_LIKE_SYSTEM_OBJECTS_PAIRS[
                    self.gravitational_system.name
                ]
                colors: list[str] = []
                for label in labels:
                    if label in plotting.SOLAR_SYSTEM_COLORS:
                        color = plotting.SOLAR_SYSTEM_COLORS[label]
                        if color is not None:
                            colors.append(color)

                kwargs["colors"] = colors
                kwargs["labels"] = labels
                kwargs["legend"] = True

        # Plot 2D trajectory
        plotting.plot_2d_trajectory(
            **kwargs,
        )
        if self.is_save_plots:
            print("Done!")

        print("------------------------------------")

    def _plot_3d_trajectory_wrapper(self) -> None:
        # Get kwargs
        kwargs: dict[str, Any] = {}
        kwargs["sol_state"] = self.sol_state_
        if self.is_save_plots:
            i = 0
            while True:
                save_fig_path = (
                    self.default_results_folder / f"3d_trajectory_{i:05d}.pdf"
                )
                if save_fig_path.exists():
                    i += 1
                else:
                    break

            kwargs["save_fig"] = True
            kwargs["save_fig_path"] = save_fig_path
            print(f"Saving plot to {save_fig_path}...")
        else:
            print("Plotting 3D trajectory...(Please check the window)")

        if "gravitational_system" in self.has_data_attr:
            if self.gravitational_system.name in self.SOLAR_LIKE_SYSTEM_OBJECTS_PAIRS:
                labels: list[str] = self.SOLAR_LIKE_SYSTEM_OBJECTS_PAIRS[
                    self.gravitational_system.name
                ]
                colors: list[str] = []
                for label in labels:
                    if label in plotting.SOLAR_SYSTEM_COLORS:
                        color = plotting.SOLAR_SYSTEM_COLORS[label]
                        if color is not None:
                            colors.append(color)

                kwargs["colors"] = colors
                kwargs["labels"] = labels
                kwargs["legend"] = True

        # Plot 3D trajectory
        plotting.plot_3d_trajectory(
            **kwargs,
        )
        if self.is_save_plots:
            print("Done!")
        print("------------------------------------")

    def _animation_get_user_input(
        self,
    ) -> tuple[float, int, str, float, int, bool, bool, bool]:
        while True:
            ### FPS ###
            fps = self.get_float("Enter FPS: ", larger_than=0.0)

            ### Plotting frequency ###
            while True:
                desired_time = self.get_float(
                    "Enter desired length for the gif (in seconds): ",
                    larger_than=0.0,
                )
                if desired_time <= self.data_size_ / fps:
                    break
                else:
                    print(
                        f"For FPS = {fps:.1f}, the maximum length is {(self.data_size_ / fps):.1f} s!"
                    )
                    print()

            print()

            plotting_freq = int(self.data_size_ / (desired_time * fps))
            print(f"Plotting frequency: {plotting_freq}")
            frame_size = int(self.data_size_ / plotting_freq) + 1
            print(f"Estimated time length: {(frame_size / fps):.1f} s")
            print()

            ### File name ###
            file_name = input("Enter file name: ").strip()
            if not file_name.endswith(".gif"):
                file_name += ".gif"

            ### DPI ###
            dpi = self.get_float(
                "Enter dots per inch (dpi) (recommended: 200): ",
                larger_than=0.0,
            )

            ### Trajectory length ###
            traj_len = self.get_int(
                "Enter the number of points for the trail (Enter -1 for full trajectory): ",
                larger_than=-2,
            )
            if traj_len > self.data_size_:
                print(
                    "Trail length is larger than the data size. Adjusted to full trajectory."
                )
                traj_len = -1

            ### Dynamic axes limits ###
            is_dynamic_axes = self.get_bool("Use dynamic axes limit? (y/n): ")

            ### Fixed step size ###
            is_maintain_fixed_dt = self.get_bool(
                "Use fixed step size for the animation? (y/n): "
            )

            ### Print summary ###
            print()
            print(f"FPS = {fps:.1f}")
            print(f"Plotting frequency: {plotting_freq}")
            print(f"Estimated frame size: about {frame_size} frames")
            print(f"Estimated time length: {(frame_size / fps):.1f} s")
            print(f"File name: {file_name}")
            print(f"dpi: {dpi:.1f}")
            if traj_len != -1:
                print(f"Trail length: {traj_len}")
            else:
                print("Trail length: full trajectory")
            print(f"Dynamic axes limits: {is_dynamic_axes}")
            print(f"Fixed step size: {is_maintain_fixed_dt}")

            ### User confirmation ###
            is_cancel = False
            if self.get_bool("Proceed? (y/n)"):
                print()
                break
            else:
                print()

            if self.get_bool("Return to menu? (y/n)"):
                is_cancel = True
                print()
                break
            else:
                print()

        return (
            fps,
            plotting_freq,
            file_name,
            dpi,
            traj_len,
            is_dynamic_axes,
            is_maintain_fixed_dt,
            is_cancel,
        )

    def _animate_2d_trajectory_wrapper(self) -> None:
        (
            fps,
            plotting_freq,
            file_name,
            dpi,
            traj_len,
            is_dynamic_axes,
            is_maintain_fixed_dt,
            is_cancel,
        ) = self._animation_get_user_input()
        if is_cancel:
            return

        # Get kwargs
        kwargs: dict[str, Any] = {
            "file_path": self.default_results_folder / file_name,
            "sol_state": self.sol_state_,
            "fps": fps,
            "plotting_freq": plotting_freq,
            "dpi": dpi,
            "is_dynamic_axes": is_dynamic_axes,
            "is_maintain_fixed_dt": is_maintain_fixed_dt,
            "sol_time": self.sol_time_,
            "traj_len": traj_len,
        }

        if "gravitational_system" in self.has_data_attr:
            if self.gravitational_system.name in self.SOLAR_LIKE_SYSTEM_OBJECTS_PAIRS:
                labels = self.SOLAR_LIKE_SYSTEM_OBJECTS_PAIRS[
                    self.gravitational_system.name
                ]
                colors: list[str] = []
                for label in labels:
                    if label in plotting.SOLAR_SYSTEM_COLORS:
                        color = plotting.SOLAR_SYSTEM_COLORS[label]
                        if color is not None:
                            colors.append(color)

                kwargs["colors"] = colors
                kwargs["labels"] = labels
                kwargs["legend"] = True

        print("Animating 2D trajectory (xy plane) in .gif...")
        plotting.animate_2d_traj_gif(**kwargs)
        print(f'Output completed! Please check "{kwargs["file_path"]}"')
        print("------------------------------------")

    def _animate_3d_trajectory_wrapper(self) -> None:
        (
            fps,
            plotting_freq,
            file_name,
            dpi,
            traj_len,
            is_dynamic_axes,
            is_maintain_fixed_dt,
            is_cancel,
        ) = self._animation_get_user_input()
        if is_cancel:
            return

        # Get kwargs
        kwargs: dict[str, Any] = {
            "file_path": self.default_results_folder / file_name,
            "sol_state": self.sol_state_,
            "fps": fps,
            "plotting_freq": plotting_freq,
            "dpi": dpi,
            "is_dynamic_axes": is_dynamic_axes,
            "is_maintain_fixed_dt": is_maintain_fixed_dt,
            "sol_time": self.sol_time_,
            "traj_len": traj_len,
        }

        if "gravitational_system" in self.has_data_attr:
            if self.gravitational_system.name in self.SOLAR_LIKE_SYSTEM_OBJECTS_PAIRS:
                labels = self.SOLAR_LIKE_SYSTEM_OBJECTS_PAIRS[
                    self.gravitational_system.name
                ]
                colors: list[str] = []
                for label in labels:
                    if label in plotting.SOLAR_SYSTEM_COLORS:
                        color = plotting.SOLAR_SYSTEM_COLORS[label]
                        if color is not None:
                            colors.append(color)

                kwargs["colors"] = colors
                kwargs["labels"] = labels
                kwargs["legend"] = True

        print("Animating 3D trajectory (xy plane) in .gif...")
        plotting.animate_3d_traj_gif(**kwargs)
        print(f'Output completed! Please check "{kwargs["file_path"]}"')
        print("------------------------------------")

    def _plot_dt_wrapper(self) -> None:
        if "sol_dt_" not in self.has_data_attr:
            print("Error: dt data not available.")
            print()
            return

        # Get kwargs
        kwargs: dict[str, Any] = {}
        if self.dt_units == "years":
            kwargs["quantity"] = self.sol_dt_ / self.simulator.DAYS_PER_YEAR
        else:
            kwargs["quantity"] = self.sol_dt_

        if self.tf_units == "years":
            kwargs["sol_time"] = self.sol_time_ / self.simulator.DAYS_PER_YEAR
        else:
            kwargs["sol_time"] = self.sol_time_

        kwargs["title"] = "dt against time"
        kwargs["xlabel"] = f"Time ({self.tf_units})"
        kwargs["ylabel"] = f"dt ({self.dt_units})"

        if self.is_save_plots:
            i = 0
            while True:
                save_fig_path = self.default_results_folder / f"dt_plot_{i:05d}.pdf"
                if save_fig_path.exists():
                    i += 1
                else:
                    break
            kwargs["save_fig"] = True
            kwargs["save_fig_path"] = save_fig_path

            print(f"Saving plot to {save_fig_path}...", end="")
        else:
            print("Plotting dt...(Please check the window)")

        # Plot dt
        plotting.plot_quantity_against_time(
            **kwargs,
        )

        if self.is_save_plots:
            print("Done!")
        print("------------------------------------")

    def _plot_rel_energy_error_wrapper(self) -> None:
        if "sol_energy_" not in self.has_data_attr:
            if "gravitational_system" in self.has_data_attr:
                self.sol_energy_ = self.simulator.compute_energy(
                    self.gravitational_system.objects_count,
                    self.gravitational_system.m,
                    self.gravitational_system.G,
                    self.sol_state_,
                    self.is_exit_ctypes_bool,
                )
                self.has_data_attr.append("sol_energy_")
            else:
                "Error: unable to compute energy without gravitational system data."
                return

        if self.sol_energy_[0] == 0.0:
            print(
                "Error: relative energy error cannot be computed as the initial energy is 0."
            )
            print()
            return

        # Get kwargs
        kwargs: dict[str, Any] = {}
        kwargs["quantity"] = np.abs(
            (self.sol_energy_ - self.sol_energy_[0]) / self.sol_energy_[0]
        )

        if self.tf_units == "years":
            kwargs["sol_time"] = self.sol_time_ / self.simulator.DAYS_PER_YEAR
        else:
            kwargs["sol_time"] = self.sol_time_

        kwargs["is_log_y"] = True
        kwargs["title"] = "Relative energy error against time"
        kwargs["xlabel"] = f"Time ({self.tf_units})"
        kwargs["ylabel"] = "$|(E(t)-E_0)/E_0|$"

        if self.is_save_plots:
            i = 0
            while True:
                save_fig_path = (
                    self.default_results_folder / f"rel_energy_error_{i:05d}.pdf"
                )
                if save_fig_path.exists():
                    i += 1
                else:
                    break

            kwargs["save_fig"] = True
            kwargs["save_fig_path"] = save_fig_path
            print(f"Saving plot to {save_fig_path}...")

        else:
            print("Plotting relative energy error...(Please check the window)")

        # Plot relative energy error
        plotting.plot_quantity_against_time(
            **kwargs,
        )

        if self.is_save_plots:
            print("Done!")
        print("------------------------------------")

    def _plot_rel_angular_momentum_error_wrapper(self) -> None:
        if "sol_angular_momentum_" not in self.has_data_attr:
            if "gravitational_system" in self.has_data_attr:
                self.sol_angular_momentum_ = self.simulator.compute_angular_momentum(
                    self.gravitational_system.objects_count,
                    self.gravitational_system.m,
                    self.sol_state_,
                    self.is_exit_ctypes_bool,
                )
                self.has_data_attr.append("sol_angular_momentum_")
            else:
                "Error: unable to compute angular momentum without gravitational system data."
                return

        if self.sol_angular_momentum_[0] == 0.0:
            print(
                "Error: relative angular momentum error cannot be computed as the initial angular momentum is 0."
            )
            print()
            return

        # Get kwargs
        kwargs: dict[str, Any] = {}
        kwargs["quantity"] = np.abs(
            (self.sol_angular_momentum_ - self.sol_angular_momentum_[0])
            / self.sol_angular_momentum_[0]
        )

        if self.tf_units == "years":
            kwargs["sol_time"] = self.sol_time_ / self.simulator.DAYS_PER_YEAR
        else:
            kwargs["sol_time"] = self.sol_time_

        kwargs["is_log_y"] = True
        kwargs["title"] = "Relative angular momentum error against time"
        kwargs["xlabel"] = f"Time ({self.tf_units})"
        kwargs["ylabel"] = "$|(L(t)-L_0)/L_0|$"

        if self.is_save_plots:
            i = 0
            while True:
                save_fig_path = (
                    self.default_results_folder
                    / f"rel_angular_momentum_error_{i:05d}.pdf"
                )
                if save_fig_path.exists():
                    i += 1
                else:
                    break

            kwargs["save_fig"] = True
            kwargs["save_fig_path"] = save_fig_path
            print(f"Saving plot to {save_fig_path}...")

        else:
            print(
                "Plotting relative angular momentum error...(Please check the window)"
            )

        # Plot relative angular momentum error
        plotting.plot_quantity_against_time(
            **kwargs,
        )

        if self.is_save_plots:
            print("Done!")
        print("------------------------------------")

    def _plot_eccentricity_wrapper(self) -> None:
        if "sol_eccentricity_" not in self.has_data_attr:
            if "gravitational_system" in self.has_data_attr:
                self.sol_eccentricity_ = self.simulator.compute_eccentricity(
                    self.gravitational_system.objects_count,
                    self.gravitational_system.m,
                    self.gravitational_system.G,
                    self.sol_state_,
                )
                self.has_data_attr.append("sol_eccentricity_")
            else:
                print(
                    "Error: unable to compute eccentricity without gravitational system data."
                )
                return

        # Get kwargs
        kwargs: dict[str, Any] = {}
        kwargs["eccentricity_or_inclination"] = self.sol_eccentricity_
        if self.tf_units == "years":
            kwargs["sol_time"] = self.sol_time_ / self.simulator.DAYS_PER_YEAR
        else:
            kwargs["sol_time"] = self.sol_time_

        kwargs["title"] = "Eccentricity against time"
        kwargs["xlabel"] = f"Time ({self.tf_units})"
        kwargs["ylabel"] = "Eccentricity"

        if "gravitational_system" in self.has_data_attr:
            if self.gravitational_system.name in self.SOLAR_LIKE_SYSTEM_OBJECTS_PAIRS:
                labels: list[str] = self.SOLAR_LIKE_SYSTEM_OBJECTS_PAIRS[
                    self.gravitational_system.name
                ][1:]  # Exclude the first object)
                colors: list[str] = []
                for label in labels:
                    if label in plotting.SOLAR_SYSTEM_COLORS:
                        color = plotting.SOLAR_SYSTEM_COLORS[label]
                        if color is not None:
                            colors.append(color)

                kwargs["colors"] = colors
                kwargs["labels"] = labels
                kwargs["legend"] = True

        if self.is_save_plots:
            i = 0
            while True:
                save_fig_path = (
                    self.default_results_folder / f"plot_eccentricity_{i:05d}.pdf"
                )
                if save_fig_path.exists():
                    i += 1
                else:
                    break

            kwargs["save_fig"] = True
            kwargs["save_fig_path"] = save_fig_path
            print(f"Saving plot to {save_fig_path}...")

        else:
            print("Plotting eccentricity...(Please check the window)")

        # Plot eccentricity
        plotting.plot_eccentricity_or_inclination(
            **kwargs,
        )

        if self.is_save_plots:
            print("Done!")
        print("------------------------------------")

    def _plot_inclination_wrapper(self) -> None:
        if "sol_inclination_" not in self.has_data_attr:
            if "gravitational_system" in self.has_data_attr:
                self.sol_inclination_ = self.simulator.compute_inclination(
                    self.gravitational_system.objects_count,
                    self.sol_state_,
                )
                self.has_data_attr.append("sol_inclination_")
            else:
                print(
                    "Error: unable to compute inclination without gravitational system data."
                )
                return

        # Get kwargs
        kwargs: dict[str, Any] = {}
        kwargs["eccentricity_or_inclination"] = self.sol_inclination_
        if self.tf_units == "years":
            kwargs["sol_time"] = self.sol_time_ / self.simulator.DAYS_PER_YEAR
        else:
            kwargs["sol_time"] = self.sol_time_

        kwargs["title"] = "Inclination against time"
        kwargs["xlabel"] = f"Time ({self.tf_units})"
        kwargs["ylabel"] = "Inclination"

        if "gravitational_system" in self.has_data_attr:
            if self.gravitational_system.name in self.SOLAR_LIKE_SYSTEM_OBJECTS_PAIRS:
                labels: list[str] = self.SOLAR_LIKE_SYSTEM_OBJECTS_PAIRS[
                    self.gravitational_system.name
                ][1:]  # Exclude the first object)
                colors: list[str] = []
                for label in labels:
                    if label in plotting.SOLAR_SYSTEM_COLORS:
                        color = plotting.SOLAR_SYSTEM_COLORS[label]
                        if color is not None:
                            colors.append(color)

                kwargs["colors"] = colors
                kwargs["labels"] = labels
                kwargs["legend"] = True

        if self.is_save_plots:
            i = 0
            while True:
                save_fig_path = (
                    self.default_results_folder / f"plot_inclination_{i:05d}.pdf"
                )
                if save_fig_path.exists():
                    i += 1
                else:
                    break

            kwargs["save_fig"] = True
            kwargs["save_fig_path"] = save_fig_path
            print(f"Saving plot to {save_fig_path}...")

        else:
            print("Plotting inclination...(Please check the window)")

        # Plot eccentricity
        plotting.plot_eccentricity_or_inclination(
            **kwargs,
        )

        if self.is_save_plots:
            print("Done!")
        print("------------------------------------")

    def _trim_data(self) -> None:
        # Get user input
        try:
            while True:
                print()
                desired_trim_size = self.get_int(
                    'Enter desired data size (Enter "cancel" to cancel): ',
                    larger_than=0,
                    smaller_than=self.data_size_,
                    allow_cancel=True,
                )

                trim_freq = int(
                    (self.data_size_ + desired_trim_size - 1) / desired_trim_size
                )
                trim_size = int(self.data_size_ / trim_freq)
                if self.get_bool(
                    f"The trimmed data size would be {trim_size}. Continue? (y/n): "
                ):
                    print()
                    break

        # User entered "cancel"
        except ValueError:
            return

        # Trim data
        for attr in self.has_data_attr:
            data = getattr(self, attr)
            if isinstance(data, np.ndarray) or isinstance(data, int):
                setattr(self, attr, utils.trim_data(data, trim_freq))

        print(f"Done! Trimmed data size: {self.data_size_}")
        print("------------------------------------")

    def _save_results(self) -> None:
        results_attr = [
            "sol_state_",
            "sol_time_",
            "sol_dt_",
        ]
        for attr in results_attr:
            if not hasattr(self, attr):
                print(f'Error: "{attr}" data not available.')
                print()
                return

        if "sol_energy_" not in self.has_data_attr:
            if "gravitational_system" in self.has_data_attr:
                if self.get_bool("Energy data not available. Compute energy? (y/n): "):
                    self.sol_energy_ = self.simulator.compute_energy(
                        self.gravitational_system.objects_count,
                        self.gravitational_system.m,
                        self.gravitational_system.G,
                        self.sol_state_,
                        self.is_exit_ctypes_bool,
                    )
                    self.has_data_attr.append("sol_energy_")
                else:
                    if not self.get_bool("Proceed saving without energy data? (y/n): "):
                        return
            else:
                if not self.get_bool(
                    "Energy data not available. Proceed saving without energy? (y/n): "
                ):
                    return

        if hasattr(self, "sol_energy_"):
            energy_data = self.sol_energy_
        else:
            energy_data = np.zeros(self.data_size_)

        results_file_path = self.default_results_folder / (
            str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")) + "_result.csv"
        )

        utils.save_results_csv(
            results_file_path,
            self.sol_state_,
            self.sol_time_,
            self.sol_dt_,
            energy_data,
        )

        print(f'Done! Results saved to "{results_file_path}".')
        print("------------------------------------")

    def _read_results(self) -> None:
        while True:
            results_file_path = Path(
                input("Enter the path / file name to the results file: ")
                .strip('"')
                .strip("'")
            )
            if results_file_path.exists():
                break
            elif (self.default_results_folder / results_file_path).exists():
                results_file_path = self.default_results_folder / results_file_path
                break

            print("File not found! Please try again.")
            print()

        print("Reading results...")
        results_dict = utils.read_results_csv(results_file_path)
        self.sol_time_ = results_dict["time"]
        self.sol_dt_ = results_dict["dt"]
        self.sol_energy_ = results_dict["energy"]
        self.sol_state_ = results_dict["state"]
        self.data_size_ = len(self.sol_state_)

        print("Done!")
        print("------------------------------------")
