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

import numpy as np

from common import get_bool
from common import get_int
from common import get_float
from simulator import Simulator
from plotter import Plotter
from progress_bar import Progress_bar


class GravitySimulator:
    DAYS_PER_YEAR = 365.242189

    def __init__(self):
        # --------------------Read command line arguments--------------------
        self._read_command_line_arg()

        # Use c library to perform simulation
        self.is_c_lib = self.args.numpy
        self.c_lib = None
        if self.is_c_lib:
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
            except:
                if get_bool("Loading c_lib failed. Run with numpy?"):
                    self.is_c_lib = False
                else:
                    print(
                        "If you want to run with C library, try to compile the "
                        + "C files provided in the src folders."
                    )
                    print("Exiting the program...")
                    sys.exit(0)

        # --------------------Initialize attributes--------------------
        self.is_exit_ctypes_bool = ctypes.c_bool(False)
        self.simulator = Simulator(self.c_lib, self.is_exit_ctypes_bool)
        self.tolerance = None
        self.dt = None
        self.default_systems = [
            "circular_binary_orbit",
            "eccentric_binary_orbit",
            "3d_helix",
            "sun_earth_moon",
            "figure-8",
            "pyth-3-body",
            "solar_system",
            "solar_system_plus",
            "custom",
        ]
        self.available_integrators = [
            "euler",
            "euler_cromer",
            "rk4",
            "leapfrog",
            "rkf45",
            "dopri",
            "dverk",
            "rkf78",
            "ias15",
        ]
        self.available_integrators_to_printable_names = {
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
        self.solar_like_systems = [
            "sun_earth_moon",
            "solar_system",
            "solar_system_plus",
        ]
        self.solar_like_systems_colors = {
            "Sun": "orange",
            "Mercury": "slategrey",
            "Venus": "wheat",
            "Earth": "skyblue",
            "Mars": "red",
            "Jupiter": "brown",
            "Saturn": "gold",
            "Uranus": "paleturquoise",
            "Neptune": "blue",
            "Moon": "grey",
            "Pluto": None,
            "Ceres": None,
            "Vesta": None,
        }
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
                    self._launch_simulation()

                else:
                    self.computed_energy = True
                    self.computed_angular_momentum = False
                    self._read_simulation_data()

                self.data_size = len(self.simulator.sol_time)

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
        action = get_int(msg, larger_than=0, smaller_than=4)
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
            + "8. Read data size\n"
            + "9. Trim data\n"
            + "10. Save simulation data\n"
            + "11. Compare relative energy error\n"
            + "12. Restart program\n"
            + "13. Exit\n"
            + "Enter action (Number): "
        )

        while True:
            action = get_int(msg, larger_than=0, smaller_than=14)
            print()

            match action:
                case 1:
                    Plotter.plot_2d_trajectory(self)
                case 2:
                    Plotter.plot_3d_trajectory(self)
                case 3:
                    (
                        fps,
                        plot_every_nth_point,
                        file_name,
                        dpi,
                        is_dynamic_axes,
                        is_custom_axes,
                        axes_lim,
                        is_cancel,
                        is_maintain_fixed_dt,
                    ) = Plotter.animation_user_interface(2, self)
                    if not is_cancel:
                        Plotter.animation_2d_traj_gif(
                            self,
                            fps,
                            plot_every_nth_point,
                            file_name,
                            dpi,
                            is_dynamic_axes,
                            is_custom_axes,
                            axes_lim,
                            is_maintain_fixed_dt,
                        )
                case 4:
                    (
                        fps,
                        plot_every_nth_point,
                        file_name,
                        dpi,
                        is_dynamic_axes,
                        is_custom_axes,
                        axes_lim,
                        is_cancel,
                        is_maintain_fixed_dt,
                    ) = Plotter.animation_user_interface(3, self)
                    if not is_cancel:
                        Plotter.animation_3d_traj_gif(
                            self,
                            fps,
                            plot_every_nth_point,
                            file_name,
                            dpi,
                            is_dynamic_axes,
                            is_custom_axes,
                            axes_lim,
                            is_maintain_fixed_dt,
                        )
                case 5:
                    if not self.computed_energy:
                        self.simulator.compute_energy()
                        self.computed_energy = True
                    Plotter.plot_rel_energy(self)
                case 6:
                    if not self.computed_angular_momentum:
                        self.simulator.compute_angular_momentum()
                        self.computed_angular_momentum = True
                    Plotter.plot_rel_angular_momentum(self)
                case 7:
                    Plotter.plot_dt(self)
                case 8:
                    print(f"There are {self.data_size} lines of data.")
                    print()
                case 9:
                    print(f"There are {self.data_size} lines of data.")
                    self.trim_data()
                case 10:
                    if not self.computed_energy:
                        self.simulator.compute_energy()
                        self.computed_energy = True
                    self._save_result()
                case 11:
                    if not self.computed_energy:
                        self.simulator.compute_energy()
                        self.computed_energy = True
                    Plotter.plot_compare_rel_energy(self)
                case 12:
                    break
                case 13:
                    print("Exiting the program...")
                    sys.exit(0)

    def _launch_simulation(self):
        while True:
            self._get_user_simulation_input()
            self._print_user_simulation_input()
            if get_bool("Proceed?"):
                print("")
                break

        self.simulator.store_every_n = self.store_every_n
        self.simulator.system = self.system
        self.simulator.integrator = self.integrator
        self.simulator.tf = self.tf
        self.simulator.unit = self.tf_unit
        self.simulator.tolerance = self.tolerance
        self.simulator.dt = self.dt

        if (not self.is_c_lib) or (self.simulator.system not in self.default_systems):
            (
                self.simulator.x,
                self.simulator.v,
                self.simulator.m,
                self.simulator.objects_count,
                self.simulator.G,
            ) = self.simulator.initialize_system_numpy(self.simulator.system)

        if self.simulator.system not in self.default_systems:
            self.simulator.simulation(is_custom_sys=True)
        else:
            self.simulator.simulation(is_custom_sys=False)

    def _get_user_simulation_input(self):
        # --------------------Check and Input systems--------------------
        while True:
            self.available_systems = self.default_systems.copy()
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
            print("\nCustomizing system...")
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
                    elif matches.group(1) in self.available_systems:
                        print("System name already exist! Please try another one.")
                    else:
                        system_name = matches.group(1)
                        break

            self.system = system_name
            objects_count = get_int("Number of objects: ", larger_than=0)
            print()

            print("Note: The default unit is M_sun, AU and day, G=0.00029591220828411.")
            masses = []
            for i in range(objects_count):
                masses.append(get_float(f"Please enter the mass for object {i + 1}: "))

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
                    x[i * 3 + j] = get_float(
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
                    v[i * 3 + j] = get_float(
                        f"Please enter v{variable} for object {i + 1}: "
                    )

            file_path = Path(__file__).parent / "customized_systems.csv"
            with open(file_path, "a", newline="") as file:
                writer = csv.writer(file)
                writer.writerow(
                    [system_name, objects_count] + masses + x.tolist() + v.tolist()
                )

        # --------------------Recommended settings for systems--------------------
        elif self.system in self.default_systems:
            print("")
            if get_bool("Do you want to use the recommended settings for this system?"):
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
            for i, integrator in enumerate(self.available_integrators):
                print(
                    f"{i + 1}. {self.available_integrators_to_printable_names[integrator]}"
                )
            self.integrator = input("Enter integrator (Number or name): ")
            # Temporary string
            temp_str = ""
            for i in range(len(self.available_integrators)):
                temp_str += "|" + self.available_integrators[i]
            if matches := re.search(
                rf"^\s*([1-9]{temp_str})\s*$",
                self.integrator,
                re.IGNORECASE,
            ):
                if matches.group(1):
                    try:
                        int(matches.group(1))
                    except ValueError:
                        if matches.group(1).lower() in self.available_integrators:
                            self.integrator = matches.group(1).lower()
                            break
                    else:
                        if (int(matches.group(1)) - 1) in range(
                            len(self.available_integrators)
                        ):
                            self.integrator = self.available_integrators[
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
            self.tf = input("Enter tf (d/yr): ")
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
        if self.integrator in ["euler", "euler_cromer", "rk4", "leapfrog"]:
            while True:
                self.dt = input("Enter dt (d/yr): ")
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

        elif self.integrator in ["rkf45", "dopri", "dverk", "rkf78", "ias15"]:
            self.tolerance = get_float("Enter tolerance: ", larger_than=0)

        # --------------------Input store_every_n--------------------
        self.store_every_n = get_int("Store every nth point (int): ", 0)
        print()

        # Exit
        return None

    def _print_user_simulation_input(self):
        print(f"System: {self.system}")
        print(
            f"Integrator: {self.available_integrators_to_printable_names[self.integrator]}"
        )
        if self.tf_unit == "years":
            print(f"tf: {self.tf / self.DAYS_PER_YEAR:g} years")
        else:
            print(f"tf: {self.tf} days")
        if self.integrator in ["euler", "euler_cromer", "rk4", "leapfrog"]:
            if self.dt_unit == "years":
                print(f"dt: {self.dt / self.DAYS_PER_YEAR} years")
            else:
                print(f"dt: {self.dt} days")
        elif self.integrator in ["rkf45", "dopri", "dverk", "rkf78", "ias15"]:
            print(f"Tolerance: {self.tolerance}")
        print(f"Use c_lib: {self.is_c_lib}")
        print(f"Store every nth point: {self.store_every_n}")

        if self.integrator in ["euler", "euler_cromer", "rk4", "leapfrog"]:
            npts = int(np.floor((self.tf / self.dt))) + 1  # + 1 for t0

            store_npts = npts
            if self.store_every_n != 1:
                store_npts = (
                    math.floor((npts - 1) / self.store_every_n) + 1
                )  # + 1 for t0

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
            desired_trim_size = get_int(
                'Enter desired data size (Enter "cancel" to cancel): ',
                larger_than=1,
                smaller_than=self.data_size,
                allow_cancel=True,
            )
            if desired_trim_size is None:  # User entered "cancel"
                is_trim_data = False
                print()
                break

            else:
                divide_factor = math.ceil(self.data_size / desired_trim_size)
                trim_size = math.floor(self.data_size / divide_factor) + 1
                if get_bool(f"The trimmed data size would be {trim_size}. Continue?"):
                    is_trim_data = True
                    print()
                    break

        # --------------------Trim data--------------------

        # Note: trim size is calculated in the the user input section:
        # divide_factor = math.ceil(self.data_size / desired_trim_size)
        # trim_size = math.floor(self.data_size / divide_factor) + 1

        if is_trim_data == True:
            trimmed_sol_time = np.zeros(trim_size)
            trimmed_sol_dt = np.zeros(trim_size)
            trimmed_sol_state = np.zeros(
                (trim_size, self.simulator.objects_count * 3 * 2)
            )
            if self.computed_energy == True:
                trimmed_energy = np.zeros(trim_size)

            j = 0
            for i in range(self.data_size):
                if i % divide_factor == 0:
                    trimmed_sol_time[j] = self.simulator.sol_time[i]
                    trimmed_sol_dt[j] = self.simulator.sol_dt[i]
                    trimmed_sol_state[j] = self.simulator.sol_state[i]
                    if self.computed_energy == True:
                        trimmed_energy[j] = self.simulator.energy[i]
                    j += 1

            if trimmed_sol_time[-1] != self.simulator.sol_time[-1]:
                trimmed_sol_time[-1] = self.simulator.sol_time[-1]
                trimmed_sol_dt[-1] = self.simulator.sol_dt[-1]
                trimmed_sol_state[-1] = self.simulator.sol_state[-1]
                if self.computed_energy == True:
                    trimmed_energy[-1] = self.simulator.energy[-1]

            self.data_size = len(trimmed_sol_time)
            print(f"Trimmed data size = {self.data_size}")
            print()

            self.simulator.sol_time = trimmed_sol_time
            self.simulator.sol_dt = trimmed_sol_dt
            self.simulator.sol_state = trimmed_sol_state
            if self.computed_energy == True:
                self.simulator.energy = trimmed_energy

            if self.tf_unit == "years":
                self.sol_time_in_tf_unit = self.simulator.sol_time / self.DAYS_PER_YEAR
            else:
                self.sol_time_in_tf_unit = self.simulator.sol_time

    def _save_result(self):
        """
        Save the result in a csv file
        Unit: Solar masses, AU, day
        Format: time, dt, total energy, x1, y1, z1, x2, y2, z2, ... vx1, vy1, vz1, vx2, vy2, vz2, ...
        """
        if not self.computed_energy:
            if get_bool(
                "WARNING: Energy has not been computed. The energy data will be stored as zeros. Proceed?"
            ):
                self.simulator.energy = np.zeros(self.data_size)
            else:
                print()
                return None

        # Estimate file size
        num_entries = 3  # Time, energy and dt data
        num_entries += (
            self.simulator.objects_count * 7
        )  # velocity * 3, position * 3, mass
        file_size = (
            num_entries * self.data_size * 18
        )  # 18 is an approximated empirical value obtained from testing
        file_size /= 1000 * 1000  # Convert to MB

        if 1 < file_size < 1000:
            if not get_bool(
                f"File size is estimated to be {file_size:.1f} MB. Continue?"
            ):
                print()
                return None
        elif 1000 <= file_size:
            if not get_bool(
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

        # Storing metadata
        with open(file_path, "w", newline="") as file:
            writer = csv.writer(file, quoting=csv.QUOTE_NONE)
            writer.writerow(
                [
                    f"# Data saved on (YYYY-MM-DD): {str(datetime.datetime.now().strftime('%Y-%m-%d'))}"
                ]
            )
            writer.writerow([f"# System Name: {self.simulator.system}"])

            try:
                integrator_name = self.available_integrators_to_printable_names[
                    self.simulator.integrator
                ]
            except KeyError:
                integrator_name = None

            writer.writerow([f"# Integrator: {integrator_name}"])
            writer.writerow([f"# Number of objects: {self.simulator.objects_count}"])
            writer.writerow([f"# Simulation time (days): {self.simulator.tf}"])
            writer.writerow([f"# dt (days): {self.simulator.dt}"])
            writer.writerow([f"# Tolerance: {self.simulator.tolerance}"])
            writer.writerow([f"# Data size: {self.data_size}"])
            writer.writerow([f"# Store every nth point: {self.store_every_n}"])
            writer.writerow([f"# Run time (s): {self.simulator.run_time}"])
            masses_str = " ".join(map(str, self.simulator.m))
            writer.writerow([f"# masses: {masses_str}"])

        progress_bar = Progress_bar()
        with progress_bar:
            with open(file_path, "a", newline="") as file:
                writer = csv.writer(file)
                for count in progress_bar.track(range(self.data_size)):
                    row = np.insert(
                        self.simulator.sol_state[count],
                        0,
                        self.simulator.energy[count],
                    )
                    row = np.insert(row, 0, self.simulator.sol_dt[count])
                    row = np.insert(row, 0, self.simulator.sol_time[count])
                    writer.writerow(row.tolist())
        print(f"Storing completed. Please check {file_path}")

        print("")

    def _read_simulation_data(self):
        self.simulator.system = None
        self.simulator.integrator = None
        self.simulator.dt = None
        self.simulator.tolerance = None
        self.store_every_n = None
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

        # Read data
        start = timeit.default_timer()
        progress_bar = Progress_bar()

        # Read metadata
        has_proper_metadata = [
            False,
            False,
            False,
        ]  # to check if objects_count, simulation time and data size can be obtained properly
        with open(read_file_path, "r") as file:
            reader = csv.reader(row for row in file if row.startswith("#"))

            for row in reader:
                if row[0].startswith("# System Name: "):
                    self.simulator.system = row[0].replace("# System Name: ", "")

                if row[0].startswith("# Integrator: "):
                    integrator = row[0].replace("# Integrator: ", "")

                    try:
                        self.simulator.integrator = list(
                            self.available_integrators_to_printable_names.keys()
                        )[
                            list(
                                self.available_integrators_to_printable_names.values()
                            ).index(integrator)
                        ]
                    except KeyError:
                        pass

                if row[0].startswith("# Number of objects: "):
                    try:
                        self.simulator.objects_count = int(
                            row[0].replace("# Number of objects: ", "")
                        )
                        has_proper_metadata[0] = True
                    except ValueError:
                        pass

                if row[0].startswith("# Simulation time (days): "):
                    try:
                        self.simulator.tf = float(
                            row[0].replace("# Simulation time (days): ", "")
                        )
                        has_proper_metadata[1] = True
                    except ValueError:
                        pass

                if row[0].startswith("# dt (seconds): "):
                    try:
                        self.simulator.dt = float(
                            row[0].replace("# dt (seconds): ", "")
                        )
                    except ValueError:
                        pass

                if row[0].startswith("# Tolerance: "):
                    try:
                        self.simulator.tolerance = float(
                            row[0].replace("# Tolerance: ", "")
                        )
                    except ValueError:
                        pass

                if row[0].startswith("# Data size: "):
                    try:
                        self.simulator.data_size = int(
                            row[0].replace("# Data size: ", "")
                        )
                        has_proper_metadata[2] = True
                    except ValueError:
                        pass

                if row[0].startswith("# Store every nth point: "):
                    try:
                        self.store_every_n = int(
                            row[0].replace("# Store every nth point: ", "")
                        )
                    except ValueError:
                        pass

                if row[0].startswith("# Run time (s): "):
                    try:
                        self.simulator.run_time = float(
                            row[0].strip("# Run time (s): ")
                        )
                    except ValueError:
                        pass

                if row[0].startswith("# masses: "):
                    try:
                        masses_str = row[0].strip("# masses: ").strip()
                        self.simulator.m = np.array(masses_str.split(" "), dtype=float)
                    except ValueError:
                        pass

        try:
            with open(read_file_path, "r") as file:
                reader = csv.reader(row for row in file if not row.startswith("#"))

                # Get object_count and file size
                if not has_proper_metadata[0]:
                    self.simulator.objects_count = round(
                        (len(next(reader)) - 3) / (3 * 2)
                    )

                if not has_proper_metadata[2]:
                    self.simulator.data_size = sum(1 for _ in file) + 1

                with progress_bar:
                    task = progress_bar.add_task("", total=self.simulator.data_size)

                    # Allocate memory
                    self.simulator.sol_time = np.zeros(50000)
                    self.simulator.sol_dt = np.zeros(50000)
                    self.simulator.energy = np.zeros(50000)
                    self.simulator.sol_state = np.zeros(
                        (50000, self.simulator.objects_count * 3 * 2)
                    )

                    print()
                    print("Reading data...")

                    i = 0
                    file.seek(0)
                    for row in reader:
                        self.simulator.sol_time[i] = row[0]
                        self.simulator.sol_dt[i] = row[1]
                        self.simulator.energy[i] = row[2]
                        for j in range(self.simulator.objects_count):
                            self.simulator.sol_state[i][j * 3 + 0] = row[3 + j * 3 + 0]
                            self.simulator.sol_state[i][j * 3 + 1] = row[3 + j * 3 + 1]
                            self.simulator.sol_state[i][j * 3 + 2] = row[3 + j * 3 + 2]
                            self.simulator.sol_state[i][
                                self.simulator.objects_count * 3 + j * 3 + 0
                            ] = row[3 + self.simulator.objects_count * 3 + j * 3 + 0]
                            self.simulator.sol_state[i][
                                self.simulator.objects_count * 3 + j * 3 + 1
                            ] = row[3 + self.simulator.objects_count * 3 + j * 3 + 1]
                            self.simulator.sol_state[i][
                                self.simulator.objects_count * 3 + j * 3 + 2
                            ] = row[3 + self.simulator.objects_count * 3 + j * 3 + 2]

                        i += 1
                        progress_bar.update(task, completed=i)

                        # Extending memory buffer
                        if i % 50000 == 0:
                            self.simulator.sol_time = np.concatenate(
                                (self.simulator.sol_time, np.zeros(50000))
                            )
                            self.simulator.sol_dt = np.concatenate(
                                (self.simulator.sol_dt, np.zeros(50000))
                            )
                            self.simulator.energy = np.concatenate(
                                (self.simulator.energy, np.zeros(50000))
                            )
                            self.simulator.sol_state = np.concatenate(
                                (
                                    self.simulator.sol_state,
                                    np.zeros(
                                        (50000, self.simulator.objects_count * 3 * 2)
                                    ),
                                )
                            )

                    self.simulator.sol_time = self.simulator.sol_time[:i]
                    self.simulator.sol_dt = self.simulator.sol_dt[:i]
                    self.simulator.energy = self.simulator.energy[:i]
                    self.simulator.sol_state = self.simulator.sol_state[:i]

                if not has_proper_metadata[1]:
                    self.simulator.tf = self.simulator.sol_time[-1]

                stop = timeit.default_timer()
                print("Reading completed.\n")
                print(f"Run time: {(stop - start):.3f} s")
                print(f"Data size: {self.simulator.data_size}")
                print("")

                while True:
                    self.tf_unit = input("Enter tf unit for plotting (d/yr): ")
                    if matches := re.search(
                        r"(day|year|d|y)", self.tf_unit, re.IGNORECASE
                    ):
                        if matches.group(1) not in ["year", "y"]:
                            self.tf_unit = "days"
                        else:
                            self.tf_unit = "years"

                        if get_bool(f"Unit for tf is {self.tf_unit}. Proceed?"):
                            print()
                            break

                    print("Invalid input. Please try again.")
                    print()

        except FileNotFoundError:
            sys.exit("Error: file is not found. Exiting the program")


if __name__ == "__main__":
    grav_plot = GravitySimulator()
    grav_plot.run_prog()
