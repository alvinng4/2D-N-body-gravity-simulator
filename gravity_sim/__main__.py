"""
Terminal user interface for N-body gravity simulator
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

import numpy as np

import common
from gravitational_system import GravitationalSystem
from simulator import Simulator
import plotting


class GravitySimulatorCLI:
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
        "whfast": "WHFast",
    }
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

    def __init__(self):
        # --------------------Read command line arguments--------------------
        self._read_command_line_arg()

        # Use c library to perform simulation
        self.is_c_lib = self.args.numpy
        if self.is_c_lib:
            self.c_lib = common.load_c_lib()

        # --------------------Initialize attributes--------------------
        self.is_exit_ctypes_bool = ctypes.c_bool(False)
        self.simulator = Simulator(self.c_lib, self.is_exit_ctypes_bool)
        self.tolerance = None
        self.dt = None
        self.solar_like_systems = [
            "sun_earth_moon",
            "solar_system",
            "solar_system_plus",
        ]
        self.recommended_settings = {
            # "template": ["tf", "tf unit", "tolerance", "store_every_n"],
            "circular_binary_orbit": [50, "days", 1e-9, 1],
            "eccentric_binary_orbit": [2.6, "years", 1e-9, 1],
            "3d_helix": [20, "days", 1e-9, 1],
            "sun_earth_moon": [1, "years", 1e-9, 1],
            "figure-8": [20, "days", 1e-9, 1],
            "pyth-3-body": [70, "days", 1e-9, 1],
            "solar_system": [200, "years", 1e-9, 1],
            "solar_system_plus": [250, "years", 1e-9, 1],
        }

    def run_prog(self):
        # Catch KeyboardInterrupt
        try:
            # Restart program once the loop is finished.
            while True:
                print("\nGravity simulator")
                print("Exit the program anytime by hitting Ctrl + C\n")

                self._user_interface_before_simulation()
                if self.is_simulate:
                    self.computed_energy = False
                    self.computed_angular_momentum = False
                    self.computed_eccentricity = False
                    self.computed_inclination = False
                    self._launch_simulation()

                else:
                    self.computed_energy = True
                    self.computed_angular_momentum = False
                    self.computed_eccentricity = False
                    self.computed_inclination = False
                    self._read_simulation_data()

                if self.tf_unit == "years":
                    self.sol_time_in_tf_unit = (
                        self.simulator.sol_time / self.DAYS_PER_YEAR
                    )
                else:
                    self.sol_time_in_tf_unit = self.simulator.sol_time

                self._user_interface_after_simulation()

        except KeyboardInterrupt:
            self.is_exit_ctypes_bool.value = True
            sys.exit("\nKeyboard Interrupt detected (Cltr + C). Exiting the program...")

    def _user_interface_before_simulation(self):
        msg = (
            "Select an action:\n"
            + "1. Launch simulation\n"
            + "2. Read simulation data\n"
            + "3. Exit\n"
            + "Enter action (Number): "
        )
        action = common.get_int(msg, larger_than=0, smaller_than=4)
        print()

        match action:
            case 1:
                self.is_simulate = True
            case 2:
                self.is_simulate = False
            case 3:
                sys.exit(0)

    def _user_interface_after_simulation(self):
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
            + "13. Compare relative energy error\n"
            + "14. Restart program\n"
            + "15. Exit\n"
            + "Enter action (Number): "
        )

        while True:
            action = common.get_int(msg, larger_than=0, smaller_than=16)
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
                    self._plot_rel_energy_wrapper()
                case 6:
                    self._plot_rel_angular_momentum_wrapper()
                case 7:
                    self._plot_dt_wrapper()
                case 8:
                    self._plot_eccentricity_wrapper()
                case 9:
                    self._plot_inclination_wrapper()
                case 10:
                    print(f"There are {self.simulator.data_size} lines of data.")
                    print()
                case 11:
                    print(f"There are {self.simulator.data_size} lines of data.")
                    self.trim_data()
                case 12:
                    if not self.computed_energy:
                        self.simulator.compute_energy()
                        self.computed_energy = True
                    self._save_results()
                case 13:
                    if not self.computed_energy:
                        self.simulator.compute_energy()
                        self.computed_energy = True
                    plotting.plot_compare_rel_energy(self)
                case 14:
                    break
                case 15:
                    print("Exiting the program...")
                    sys.exit(0)

    def _launch_simulation(self):
        while True:
            self._get_user_simulation_input()
            self._print_user_simulation_input()
            if common.get_bool("Proceed?"):
                print("")
                break

        self.gravitational_system = GravitationalSystem()
        self.gravitational_system.load(self.system)

        self.simulator.launch_simulation(
            self.gravitational_system,
            self.integrator,
            self.tf,
            dt=self.dt,
            tolerance=self.tolerance,
            store_every_n=self.store_every_n,
        )

    def _get_user_simulation_input(self):
        # --------------------Check and Input systems--------------------
        while True:
            self.available_systems = GravitationalSystem.DEFAULT_SYSTEMS.copy()
            file_path = Path(__file__).parent / "customized_systems.csv"
            with open(file_path, "a"):  # Create file if not exist
                pass
            with open(file_path, "r+") as file:
                reader = csv.reader(file)
                for row in reader:
                    self.available_systems.append(row[0])

            print("Available systems:")
            for i, system in enumerate(self.available_systems):
                print(f"{i + 1}. {system}")
            self.system = input("Enter system (Number or name): ")
            # Temporary string
            temp_str = ""
            for i in range(len(self.available_systems)):
                temp_str += "|" + self.available_systems[i]
            if matches := re.search(
                rf"^\s*([1-9][0-9]*|{temp_str})\s*$",
                self.system,
                re.IGNORECASE,
            ):
                if matches.group(1):
                    try:
                        int(matches.group(1))
                    except ValueError:
                        if matches.group(1).lower() in self.available_systems:
                            self.system = matches.group(1).lower()
                            break
                    else:
                        if (int(matches.group(1)) - 1) in range(
                            len(self.available_systems)
                        ):
                            self.system = self.available_systems[
                                int(matches.group(1)) - 1
                            ]
                            break

            print("\nInvalid input. Please try again.\n")

        # --------------------Customize system--------------------
        if self.system == "custom":
            msg = (
                "\n" +
                "Note: For command line interface, only Cartesian coordinates are supported!\n" +
                "If you want to customize a system with orbital elements, please use our API.\n" +
                "\n" +
                "Customizing system..."
            )
            print(msg)
            custom_sys = GravitationalSystem()
            while True:
                system_name = input("Enter the name of the system: ")
                if matches := re.search(
                    r"^\s*(\S+)\s*$",
                    system_name,
                ):
                    if "," in matches.group(1):
                        print(
                            "Invalid input. Please do not use comma (,) inside the name."
                        )
                    else:
                        system_name = matches.group(1)
                        break

            self.system = system_name
            custom_sys.name = system_name

            objects_count = common.get_int("Number of objects: ", larger_than=0)
            print()

            print("Note: The default unit is M_sun, AU and day, G=0.00029591220828411.")
            while True:
                try:
                    custom_G = input(
                        "Enter G (Simply press enter if you want to use default value): "
                    )
                    if custom_G == "":
                        break
                    else:
                        custom_sys.G = float(custom_G)
                        break
                except ValueError:
                    print("Invalid input. Please try again.")
                    print()

            masses = []
            for i in range(objects_count):
                masses.append(
                    common.get_float(f"Please enter the mass for object {i + 1}: ")
                )

            x = np.zeros(objects_count * 3)
            for i in range(objects_count):
                for j in range(3):
                    match j:
                        case 0:
                            variable = "x"
                        case 1:
                            variable = "y"
                        case 2:
                            variable = "z"
                    x[i * 3 + j] = common.get_float(
                        f"Please enter {variable} for object {i + 1}: "
                    )

            v = np.zeros(objects_count * 3)
            for i in range(objects_count):
                for j in range(3):
                    match j:
                        case 0:
                            variable = "x"
                        case 1:
                            variable = "y"
                        case 2:
                            variable = "z"
                    v[i * 3 + j] = common.get_float(
                        f"Please enter v{variable} for object {i + 1}: "
                    )

            for i in range(objects_count):
                custom_sys.add(
                    x[i * 3 : i * 3 + 3],
                    v[i * 3 : i * 3 + 3],
                    masses[i],
                )
            custom_sys.save()

        # --------------------Recommended settings for systems--------------------
        elif self.system in GravitationalSystem.DEFAULT_SYSTEMS:
            print("")
            if common.get_bool(
                "Do you want to use the recommended settings for this system?"
            ):
                print("")
                self.integrator = "ias15"
                (
                    self.tf,
                    self.tf_unit,
                    self.tolerance,
                    self.store_every_n,
                ) = self.recommended_settings[self.system]
                if self.tf_unit == "years":
                    self.tf *= self.DAYS_PER_YEAR

                return None

        print()

        # If user did not choose recommended settings:
        # --------------------Input integrators--------------------
        while True:
            print("Available integrators: ")
            for i, integrator in enumerate(Simulator.AVAILABLE_INTEGRATORS):
                print(
                    f"{i + 1}. {GravitySimulatorCLI.AVAILABLE_INTEGRATORS_TO_PRINTABLE_NAMES[integrator]}"
                )
            self.integrator = input("Enter integrator (Number or name): ")
            # Temporary string
            temp_str = ""
            for i in range(len(Simulator.AVAILABLE_INTEGRATORS)):
                temp_str += "|" + Simulator.AVAILABLE_INTEGRATORS[i]
            if matches := re.search(
                rf"^\s*([1-9][0-9]*|{temp_str})\s*$",
                self.integrator,
                re.IGNORECASE,
            ):
                if matches.group(1):
                    try:
                        int(matches.group(1))
                    except ValueError:
                        if matches.group(1).lower() in Simulator.AVAILABLE_INTEGRATORS:
                            self.integrator = matches.group(1).lower()
                            break
                    else:
                        if (int(matches.group(1)) - 1) in range(
                            len(Simulator.AVAILABLE_INTEGRATORS)
                        ):
                            self.integrator = Simulator.AVAILABLE_INTEGRATORS[
                                int(matches.group(1)) - 1
                            ]
                            break

            print()
            print("Invalid input. Please try again.")
            print()

        print()
        # --------------------Input tf--------------------
        # tf = 0 is allowed as user may want to plot the
        # initial position of the system
        while True:
            self.tf = input("Enter tf (days/year) (e.g. 200y): ")
            if matches := re.search(
                r"([0-9]*\.?[0-9]*)(?:\.|\W*)*(day|year|d|y)?", self.tf, re.IGNORECASE
            ):
                if not matches.group(1):
                    print("Invalid input. Please try again.")
                    print()
                    continue

                try:
                    self.tf = float(matches.group(1))
                    if self.tf < 0:
                        raise ValueError
                except ValueError:
                    print("Invalid input. Please try again.")
                    print()
                    continue

                if matches.group(2) not in ["year", "y"]:
                    self.tf_unit = "days"
                    break
                else:
                    self.tf_unit = "years"
                    self.tf *= self.DAYS_PER_YEAR
                    break

        # --------------------Input dt--------------------
        if self.integrator in self.simulator.FIXED_STEP_SIZE_INTEGRATORS:
            while True:
                self.dt = input("Enter dt (days/year) (e.g. 1d): ")
                if matches := re.search(
                    r"([0-9]*\.?[0-9]*)(?:\.|\W*)*(day|year|d|y)?\s*",
                    self.dt,
                    re.IGNORECASE,
                ):
                    if not matches.group(1):
                        print("Invalid input. Please try again.")
                        print()
                        continue

                    try:
                        self.dt = float(matches.group(1))
                        if self.dt <= 0:
                            raise ValueError
                    except ValueError:
                        print("Invalid input. Please try again.")
                        print()
                        continue

                    if matches.group(2) not in ["year", "y"]:
                        self.dt_unit = "days"
                        break
                    else:
                        self.dt_unit = "years"
                        self.dt *= self.DAYS_PER_YEAR
                        break

        elif self.integrator in self.simulator.ADAPTIVE_STEP_SIZE_INTEGRATORS:
            self.tolerance = common.get_float("Enter tolerance: ", larger_than=0)

        # --------------------Input store_every_n--------------------
        self.store_every_n = common.get_int("Store every nth point (int): ", 0)
        print()

        # Exit
        return None

    def _print_user_simulation_input(self):
        print(f"System: {self.system}")
        print(
            f"Integrator: {GravitySimulatorCLI.AVAILABLE_INTEGRATORS_TO_PRINTABLE_NAMES[self.integrator]}"
        )
        if self.tf_unit == "years":
            print(f"tf: {self.tf / self.DAYS_PER_YEAR:g} years")
        else:
            print(f"tf: {self.tf} days")
        if self.integrator in self.simulator.FIXED_STEP_SIZE_INTEGRATORS:
            if self.dt_unit == "years":
                print(f"dt: {self.dt / self.DAYS_PER_YEAR} years")
            else:
                print(f"dt: {self.dt} days")
        elif self.integrator in self.simulator.ADAPTIVE_STEP_SIZE_INTEGRATORS:
            print(f"Tolerance: {self.tolerance}")
        print(f"Use c_lib: {self.is_c_lib}")
        print(f"Store every nth point: {self.store_every_n}")

        if self.integrator in self.simulator.FIXED_STEP_SIZE_INTEGRATORS:
            npts = math.ceil(self.tf / self.dt)
            if self.store_every_n != 1:
                store_npts = math.floor(npts / self.store_every_n)
            else:
                store_npts = npts
            store_npts += 1  # + 1 for t0

            print(f"Estimated number of points to be stored: {store_npts}")

        print("")

    def _read_command_line_arg(self):
        parser = argparse.ArgumentParser(description="N-body gravity simulator")
        parser.add_argument(
            "--numpy",
            "-n",
            action="store_false",
            help="disable c_lib and use numpy",
        )
        self.args = parser.parse_args()

    def trim_data(self):
        # --------------------Get user input--------------------
        while True:
            print()
            desired_trim_size = common.get_int(
                'Enter desired data size (Enter "cancel" to cancel): ',
                larger_than=1,
                smaller_than=self.simulator.data_size,
                allow_cancel=True,
            )
            if desired_trim_size is None:  # User entered "cancel"
                print()
                return None

            else:
                store_every_n = math.ceil(self.simulator.data_size / desired_trim_size)
                trim_size = math.ceil(self.simulator.data_size / store_every_n)
                if common.get_bool(
                    f"The trimmed data size would be {trim_size}. Continue?"
                ):
                    print()
                    break

        # --------------------Trim data--------------------
        data_size = self.simulator.data_size
        self.simulator.sol_time = common.trim_data(
            store_every_n, data_size, self.simulator.sol_time
        )
        self.simulator.sol_dt = common.trim_data(
            store_every_n, data_size, self.simulator.sol_dt
        )
        self.simulator.sol_state = common.trim_data(
            store_every_n, data_size, self.simulator.sol_state
        )
        if self.computed_energy:
            self.simulator.energy = common.trim_data(
                store_every_n, data_size, self.simulator.energy
            )

        print(f"Trimmed data size = {self.simulator.data_size}")
        print()

        if self.tf_unit == "years":
            self.sol_time_in_tf_unit = self.simulator.sol_time / self.DAYS_PER_YEAR
        else:
            self.sol_time_in_tf_unit = self.simulator.sol_time

    def _save_results(self):
        """
        Save the results in a csv file
        Unit: Solar masses, AU, day
        Format: time, G, dt, total energy, x1, y1, z1, x2, y2, z2, ... vx1, vy1, vz1, vx2, vy2, vz2, ...
        """
        if not self.computed_energy:
            if common.get_bool(
                "WARNING: Energy has not been computed. The energy data will be stored as zeros. Proceed?"
            ):
                self.simulator.energy = np.zeros(self.simulator.data_size)
            else:
                print()
                return None

        # Estimate file size
        num_entries = 3  # Time, energy and dt data
        num_entries += (
            self.simulator.objects_count * 7
        )  # velocity * 3, position * 3, mass
        file_size = (
            num_entries * self.simulator.data_size * 18
        )  # 18 is an approximated empirical value obtained from testing
        file_size /= 1000 * 1000  # Convert to MB

        if 1 < file_size < 1000:
            if not common.get_bool(
                f"File size is estimated to be {file_size:.1f} MB. Continue?"
            ):
                print()
                return None
        elif 1000 <= file_size:
            if not common.get_bool(
                f"File size is estimated to be {(file_size / 1000):.1f} GB. Continue?"
            ):
                print()
                return None

        print()

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
            integrator_name = (
                GravitySimulatorCLI.AVAILABLE_INTEGRATORS_TO_PRINTABLE_NAMES[
                    self.simulator.integrator
                ]
            )
        except KeyError:
            integrator_name = None

        common.save_results(
            file_path,
            self.simulator.system_name,
            integrator_name,
            self.simulator.objects_count,
            self.simulator.G,
            self.simulator.tf,
            self.simulator.dt,
            self.simulator.tolerance,
            self.simulator.data_size,
            self.store_every_n,
            self.simulator.run_time,
            self.simulator.m,
            self.simulator.sol_state,
            self.simulator.sol_time,
            self.simulator.sol_dt,
            self.simulator.energy,
        )
        print("")

    def _read_simulation_data(self):
        self.gravitational_system = GravitationalSystem()
        read_folder_path = Path(__file__).parent / "results"

        while True:
            read_file_path = input(
                "Enter absolute path to the file, or the complete file name if it is inside gravity_plot/results: "
            ).strip()
            if not read_file_path.endswith(".csv"):
                read_file_path += ".csv"
            read_file_path = read_folder_path / read_file_path

            if read_file_path.is_file():
                break
            else:
                print("File not found! Please try again.")

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
            _,
        ) = common.read_results(file_path=read_file_path)
        self.gravitational_system.name = self.simulator.system_name
        self.gravitational_system.objects_names = [
            None for _ in range(self.simulator.objects_count)
        ]
        while True:
            self.tf_unit = input("Enter tf unit for plotting (d/yr): ")
            if matches := re.search(r"(day|year|d|y)", self.tf_unit, re.IGNORECASE):
                if matches.group(1) not in ["year", "y"]:
                    self.tf_unit = "days"
                else:
                    self.tf_unit = "years"

                if common.get_bool(f"Unit for tf is {self.tf_unit}. Proceed?"):
                    print()
                    break

            print("Invalid input. Please try again.")
            print()

    def _plot_2d_trajectory_wrapper(self):
        print("Plotting 2D trajectory (xy plane)...(Please check the window)")

        colors = []
        labels = []
        legend = False

        for objects_name in self.gravitational_system.objects_names:
            try:
                colors.append(self.SOLAR_SYSTEM_COLORS[objects_name])
                legend = True  # Show legend if one of the name is recognized
            except:
                colors.append(None)

            labels.append(objects_name)

        plotting.plot_2d_trajectory(
            self.simulator.objects_count,
            self.simulator.sol_state,
            colors=colors,
            labels=labels,
            legend=legend,
        )
        print()

    def _plot_3d_trajectory_wrapper(self):
        print("Plotting 3D trajectory...(Please check the window)")

        colors = []
        labels = []
        legend = False

        for objects_name in self.gravitational_system.objects_names:
            try:
                colors.append(self.SOLAR_SYSTEM_COLORS[objects_name])
                legend = True  # Show legend if one of the name is recognized
            except:
                colors.append(None)

            labels.append(objects_name)

        plotting.plot_3d_trajectory(
            self.simulator.objects_count,
            self.simulator.sol_state,
            colors=colors,
            labels=labels,
            legend=legend,
        )

        print()

    def _animation_get_user_input(self: int):
        while True:
            fps = common.get_float("Enter FPS: ", larger_than=0)

            while True:
                desired_time = common.get_float(
                    "Enter desired time length for the gif (in seconds): ",
                    larger_than=0,
                )
                if desired_time <= self.simulator.data_size / fps:
                    break
                else:
                    print(
                        f"For FPS = {fps:.1f}, the maximum length is {(self.simulator.data_size / fps):.1f} s!"
                    )
                    print()

            print()

            plot_every_nth_point = math.floor(
                self.simulator.data_size / (desired_time * fps)
            )
            print(f"Plot every nth point: {plot_every_nth_point}")
            frame_size = math.floor(self.simulator.data_size / plot_every_nth_point) + 1
            print(f"Estimated time length: {(frame_size / fps):.1f} s")
            print()

            file_name = input("Enter file name: ").strip()

            dpi = common.get_float(
                "Enter dots per inch (dpi) (recommended: 200): ",
                larger_than=0,
            )
            traj_len = common.get_int(
                "Enter the number of points for the trail (Enter -1 for full trajectory): ",
                larger_than=-2,
                smaller_than=self.simulator.data_size,
            )
            is_dynamic_axes = common.get_bool("Use dynamic axes limit?")
            is_maintain_fixed_dt = common.get_bool("Use fixed step size for the animation?")

            print()

            print(f"FPS = {fps:.1f}")
            print(f"Plot every nth point: {plot_every_nth_point}")
            print(f"Estimated frame size: about {frame_size} frames")
            print(f"Estimated time length: {(frame_size / fps):.1f} s")
            print(f"File name: {file_name}")
            print(f"dpi: {dpi:.1f}")
            print(f"Trail length: {traj_len}")
            print(f"Dynamic axes limits: {is_dynamic_axes}")
            print(f"Fixed step size: {is_maintain_fixed_dt}")

            is_cancel = False
            if common.get_bool("Proceed?"):
                print()
                break
            else:
                print()

            if common.get_bool("Return to menu?"):
                is_cancel = True
                print()
                break
            else:
                print()

        return (
            fps,
            plot_every_nth_point,
            file_name,
            dpi,
            traj_len,
            is_dynamic_axes,
            is_maintain_fixed_dt,
            is_cancel,
        )

    def _animate_2d_trajectory_wrapper(self):
        (
            fps,
            plot_every_nth_point,
            file_name,
            dpi,
            traj_len,
            is_dynamic_axes,
            is_maintain_fixed_dt,
            is_cancel,
        ) = self._animation_get_user_input()
        if not is_cancel:
            print("Animating 2D trajectory (xy plane) in .gif...")

        try:
            colors = [self.SOLAR_SYSTEM_COLORS[objects_name] for objects_name in self.gravitational_system.objects_names]
            labels = [objects_name for objects_name in self.gravitational_system.objects_names]
            legend = True
        except KeyError:
            colors = None
            labels = None
            legend = False

        plotting.animate_2d_traj_gif(
            self.simulator.objects_count,
            self.simulator.sol_state,
            fps,
            plot_every_nth_point,
            dpi,
            traj_len=traj_len,
            is_dynamic_axes=is_dynamic_axes,
            is_maintain_fixed_dt=is_maintain_fixed_dt,
            colors=colors,
            labels=labels,
            legend=legend,
            file_name=file_name,
            sol_time=self.simulator.sol_time,
        )
        print()

    def _animate_3d_trajectory_wrapper(self):
        (
            fps,
            plot_every_nth_point,
            file_name,
            dpi,
            traj_len,
            is_dynamic_axes,
            is_maintain_fixed_dt,
            is_cancel,
        ) = self._animation_get_user_input()
        if not is_cancel:
            print("Animating 3D trajectory in .gif...")

        try:
            colors = [self.SOLAR_SYSTEM_COLORS[objects_name] for objects_name in self.gravitational_system.objects_names]
            labels = [objects_name for objects_name in self.gravitational_system.objects_names]
            legend = True
        except KeyError:
            colors = None
            labels = None
            legend = False

        plotting.animate_3d_traj_gif(
            self.simulator.objects_count,
            self.simulator.sol_state,
            fps,
            plot_every_nth_point,
            dpi,
            traj_len=traj_len,
            is_dynamic_axes=is_dynamic_axes,
            is_maintain_fixed_dt=is_maintain_fixed_dt,
            colors=colors,
            labels=labels,
            legend=legend,
            file_name=file_name,
            sol_time=self.simulator.sol_time,
        )
        print()

    def _plot_rel_energy_wrapper(self):
        if not self.computed_energy:
            self.simulator.compute_energy()
            self.computed_energy = True
        print("Plotting relative energy error...(Please check the window)")
        plotting.plot_rel_energy(
            self.simulator.energy,
            self.sol_time_in_tf_unit,
            xlabel=f"Time ({self.tf_unit})",
        )
        print()

    def _plot_rel_angular_momentum_wrapper(self):
        if not self.computed_angular_momentum:
            self.simulator.compute_angular_momentum()
            self.computed_angular_momentum = True

        print("Plotting relative angular momentum error...(Please check the window)")
        plotting.plot_rel_angular_momentum(
            self.simulator.angular_momentum,
            self.sol_time_in_tf_unit,
            xlabel=f"Time ({self.tf_unit})",
        )
        print()

    def _plot_dt_wrapper(self):
        print("Plotting dt...(Please check the window)")
        plotting.plot_dt(
            self.simulator.sol_dt,
            self.sol_time_in_tf_unit,
            xlabel=f"Time ({self.tf_unit})",
            ylabel="dt (days)",
        )
        print()

    def _plot_eccentricity_wrapper(self) -> None:
        if not self.computed_eccentricity:
            self.simulator.compute_eccentricity()
            self.computed_eccentricity = True

        try:
            colors = [self.SOLAR_SYSTEM_COLORS[objects_name] for objects_name in self.gravitational_system.objects_names[1:]]
            labels = [objects_name for objects_name in self.gravitational_system.objects_names[1:]]
            legend = True
        except KeyError:
            colors = None
            labels = None
            legend = False

        print("Plotting eccentricity...(Please check the window)")
        plotting.plot_eccentricity(
            self.simulator.eccentricity,
            self.sol_time_in_tf_unit,
            colors=colors,
            labels=labels,
            legend=legend,
            xlabel=f"Time ({self.tf_unit})",
        )
        print()

    def _plot_inclination_wrapper(self) -> None:
        if not self.computed_inclination:
            self.simulator.compute_inclination()
            self.computed_inclination = True

        try:
            colors = [self.SOLAR_SYSTEM_COLORS[objects_name] for objects_name in self.gravitational_system.objects_names[1:]]
            labels = [objects_name for objects_name in self.gravitational_system.objects_names[1:]]
            legend = True
        except KeyError:
            colors = None
            labels = None
            legend = False

        print("Plotting inclination...(Please check the window)")
        plotting.plot_inclination(
            self.simulator.inclination,
            self.sol_time_in_tf_unit,
            colors=colors,
            labels=labels,
            legend=legend,
            xlabel=f"Time ({self.tf_unit})",
        )
        print()


if __name__ == "__main__":
    grav_plot = GravitySimulatorCLI()
    grav_plot.run_prog()
