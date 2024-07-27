import csv
import datetime
from pathlib import Path
import re
import warnings

import numpy as np

from progress_bar import Progress_bar

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


def get_bool(msg: str) -> bool:
    """
    Get boolean input from user

    Args: msg (str): Message to display to user
    Print: "{msg} (y/n): "
    Returns: bool
    """
    while True:
        if matches := re.search(
            r"^\s*(yes|no|y|n)\s*$", input(f"{msg} (y/n): "), re.IGNORECASE
        ):
            if matches.group(1).lower() in ["y", "yes"]:
                return True
            elif matches.group(1).lower() in ["n", "no"]:
                return False

        print("Invalid input. Please try again.\n")


def get_int(
    msg: str,
    larger_than: int = None,
    smaller_than: int = None,
    allow_cancel: bool = False,
) -> int:
    """
    Get integer input from user

    Args: msg (str): Message to display to user
          larger_than (int): Value must be larger than this
          smaller_than (int): Value must be smaller than this
          allow_cancel (bool): Allow user to enter "cancel" to cancel input
    Print: "{msg}"
    Returns: int, None if allow_cancel is True and user enters "cancel"
    """
    while True:
        try:
            value = input(f"{msg}")
            if allow_cancel and value.strip().lower() == "cancel":
                return None
            else:
                value = int(value)

            if (larger_than is None or value > larger_than) and (
                smaller_than is None or value < smaller_than
            ):
                return value

            if larger_than is not None and value <= larger_than:
                print("Value too small! Please try again.")
                print()

            elif smaller_than is not None and value >= smaller_than:
                print("Value too big! Please try again.")
                print()

        except ValueError:
            print("Invalid input. Please try again.")
            print()


def get_float(
    msg: str,
    larger_than: float = None,
    smaller_than: float = None,
    allow_cancel: bool = False,
) -> float:
    """
    Get integer input from user

    Args: msg (str): Message to display to user
          larger_than (float): Value must be larger than this
          smaller_than (float): Value must be smaller than this
          allow_cancel (bool): Allow user to enter "cancel" to cancel input
    Print: "{msg}"
    Returns: float, None if allow_cancel is True and user enters "cancel"
    """
    while True:
        try:
            value = input(f"{msg}")
            if allow_cancel and value.strip().lower() == "cancel":
                return None
            else:
                value = float(value)

            if (larger_than is None or value > larger_than) and (
                smaller_than is None or value < smaller_than
            ):
                return value

            if larger_than is not None and value <= larger_than:
                print("Value too small! Please try again.")
                print()

            elif smaller_than is not None and value >= smaller_than:
                print("Value too big! Please try again.")
                print()

        except ValueError:
            print("Invalid input. Please try again.")
            print()


def acceleration(objects_count, x, m, G):
    """
    Calculate acceleration by a = - GM/r^3 vec{r} using a matrix approach

    Old version of code using nested for loops:
    # Allocate memory
    a = np.zeros((objects_count, 3))

    # Calculations
    for j in range(objects_count):
        for k in range(j + 1, objects_count):
            R = x[j] - x[k]
            temp_value = G * R / np.linalg.norm(R) ** 3
            a[j] += -temp_value * m[k]
            a[k] += temp_value * m[j]

    return a
    """
    # Compute the displacement vector
    r_ij = x[np.newaxis, :, :] - x[:, np.newaxis, :]

    # Compute the distance squared
    r_squared = np.sum(r_ij**2, axis=2)
    softening = 1e-12  # To prevent r = 0
    r_squared[r_squared < softening] = softening

    # Compute 1/r^3
    inv_r_cubed = r_squared ** (-1.5)

    # Set diagonal elements to 0 to avoid self-interaction
    np.fill_diagonal(inv_r_cubed, 0.0)

    # Compute the acceleration
    a = G * np.sum(
        r_ij * inv_r_cubed[:, :, np.newaxis] * m[np.newaxis, :, np.newaxis], axis=1
    )

    return a


def save_results(
    file_path: str,
    system_name: str,
    integrator_name: str,
    objects_count: int,
    G: float,
    tf: float,
    dt: float,
    tolerance: float,
    data_size: int,
    store_every_n: int,
    run_time: float,
    masses: np.ndarray,
    sol_state: np.ndarray = None,
    sol_time: np.ndarray = None,
    sol_dt: np.ndarray = None,
    energy: np.ndarray = None,
    only_metadata: bool = False,
    no_print: bool = False,
    no_progress_bar: bool = False,
):
    # Storing metadata
    with open(file_path, "w", newline="") as file:
        writer = csv.writer(file, quoting=csv.QUOTE_NONE)
        writer.writerow(
            [
                f"# Data saved on (YYYY-MM-DD): {str(datetime.datetime.now().strftime('%Y-%m-%d'))}"
            ]
        )
        writer.writerow([f"# System Name: {system_name}"])
        writer.writerow([f"# Integrator: {integrator_name}"])
        writer.writerow([f"# Number of objects: {objects_count}"])
        writer.writerow([f"# Gravitational constant: {G}"])
        writer.writerow([f"# Simulation time (days): {tf}"])
        writer.writerow([f"# dt (days): {dt}"])
        writer.writerow([f"# Tolerance: {tolerance}"])
        writer.writerow([f"# Data size: {data_size}"])
        writer.writerow([f"# Store every nth point: {store_every_n}"])
        writer.writerow([f"# Run time (s): {run_time}"])
        masses_str = " ".join(map(str, masses))
        writer.writerow([f"# masses: {masses_str}"])

    if not only_metadata:
        if not no_progress_bar:
            progress_bar = Progress_bar()
            with progress_bar:
                with open(file_path, "a", newline="") as file:
                    writer = csv.writer(file)
                    for count in progress_bar.track(range(data_size)):
                        row = np.insert(
                            sol_state[count],
                            0,
                            energy[count],
                        )
                        row = np.insert(row, 0, sol_dt[count])
                        row = np.insert(row, 0, sol_time[count])
                        writer.writerow(row.tolist())
        else:
            with open(file_path, "a", newline="") as file:
                writer = csv.writer(file)
                for count in range(data_size):
                    row = np.insert(
                        sol_state[count],
                        0,
                        energy[count],
                    )
                    row = np.insert(row, 0, sol_dt[count])
                    row = np.insert(row, 0, sol_time[count])
                    writer.writerow(row.tolist())

    if not no_print:
        print(f"Storing completed. Please check {file_path}")


def read_results(
    file_path: str,
    start: int = 0,
    end: int = -1,
    step: int = 1,
    memory_buffer_size: int = 50000,
    no_print: bool = False,
    no_progress_bar: bool = False,
):
    """
    Read the results from a csv file

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

    Returns
    -------
    sol_state : np.ndarray
        Solution state
    sol_time : np.ndarray
        Solution time
    sol_dt : np.ndarray
        Solution dt
    energy : np.ndarray
        Total energy of each states
    system_name : str
        Name of the system
    integrator : str
        Name of the integrator
    objects_count : int
        Number of objects
    G : float
        Gravitational constant
    tf : float
        Simulation time
    dt : float
        Time step size
    tolerance : float
        Tolerance for adaptive step size integrators
    store_every_n : int
        Store every nth point
    run_time : float
        Run time of the simulation
    m : np.ndarray
        Masses of the objects
    store_count : int
        New data size
    """

    if start < 0 or end < -1 or step < 1:
        raise ValueError("Invalid start, end or step values.")

    system_name = None
    integrator = None
    dt = None
    tolerance = None
    store_every_n = None

    if not Path(file_path).is_file():
        raise FileNotFoundError(f"File {file_path} not found.")

    progress_bar = Progress_bar()

    # Read metadata
    metadata = []

    # to check if objects_count, simulation time and data size can be obtained properly
    has_proper_metadata = [False, False, False]

    with open(file_path, "r") as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0].startswith("#"):
                metadata.append(row)
            else:
                break

    for row in metadata:
        if row[0].startswith("# System Name: "):
            system_name = row[0].replace("# System Name: ", "")

        elif row[0].startswith("# Integrator: "):
            integrator = row[0].replace("# Integrator: ", "")

            try:
                integrator = list(AVAILABLE_INTEGRATORS_TO_PRINTABLE_NAMES.keys())[
                    list(AVAILABLE_INTEGRATORS_TO_PRINTABLE_NAMES.values()).index(
                        integrator
                    )
                ]
            except KeyError:
                pass

        elif row[0].startswith("# Number of objects: "):
            try:
                objects_count = int(row[0].replace("# Number of objects: ", ""))
                has_proper_metadata[0] = True
            except ValueError:
                pass

        elif row[0].startswith("# Gravitational constant: "):
            try:
                G = float(row[0].replace("# Gravitational constant: ", ""))
            except ValueError:
                pass

        elif row[0].startswith("# Simulation time (days): "):
            try:
                tf = float(row[0].replace("# Simulation time (days): ", ""))
                has_proper_metadata[1] = True
            except ValueError:
                pass

        elif row[0].startswith("# dt (days): "):
            try:
                dt = float(row[0].replace("# dt (days): ", ""))
            except ValueError:
                pass

        elif row[0].startswith("# Tolerance: "):
            try:
                tolerance = float(row[0].replace("# Tolerance: ", ""))
            except ValueError:
                pass

        elif row[0].startswith("# Data size: "):
            try:
                data_size = int(row[0].replace("# Data size: ", ""))
                has_proper_metadata[2] = True
            except ValueError:
                pass

        elif row[0].startswith("# Store every nth point: "):
            try:
                store_every_n = int(row[0].replace("# Store every nth point: ", ""))
            except ValueError:
                pass

        elif row[0].startswith("# Run time (s): "):
            try:
                run_time = float(row[0].strip("# Run time (s): "))
            except ValueError:
                pass

        elif row[0].startswith("# masses: "):
            try:
                masses_str = row[0].strip("# masses: ").strip()
                m = np.array(masses_str.split(" "), dtype=float)
            except ValueError:
                pass

    if not has_proper_metadata[2]:
        if not no_print:
            print(
                "Original data size in metadata is not found. Disabling progress bar."
            )
        no_progress_bar = True

    # Read data
    with open(file_path, "r") as file:
        reader = csv.reader(file)

        for row in reader:
            if row[0].startswith("#"):
                continue

            # Get object_count
            if not has_proper_metadata[0]:
                objects_count = (len(row) - 3) // (3 * 2)

            break

        file.seek(0)
        with progress_bar:
            if not no_progress_bar:
                task = progress_bar.add_task("", total=data_size)

            # Allocate memory
            sol_time = np.zeros(memory_buffer_size)
            sol_dt = np.zeros(memory_buffer_size)
            energy = np.zeros(memory_buffer_size)
            sol_state = np.zeros((memory_buffer_size, objects_count * 3 * 2))

            if not no_print:
                print()
                print("Reading data...")

            count = 0
            store_count = 0
            file.seek(0)
            for row in reader:
                if row[0].startswith("#"):
                    continue

                if count < start:
                    count += 1
                    continue

                if end != -1 and (count - start) >= end:
                    break

                if (count - start) % step != 0:
                    count += 1
                    continue

                sol_time[store_count] = row[0]
                sol_dt[store_count] = row[1]
                energy[store_count] = row[2]
                for j in range(objects_count):
                    sol_state[store_count][j * 3 + 0] = row[3 + j * 3 + 0]
                    sol_state[store_count][j * 3 + 1] = row[3 + j * 3 + 1]
                    sol_state[store_count][j * 3 + 2] = row[3 + j * 3 + 2]
                    sol_state[store_count][objects_count * 3 + j * 3 + 0] = row[
                        3 + objects_count * 3 + j * 3 + 0
                    ]
                    sol_state[store_count][objects_count * 3 + j * 3 + 1] = row[
                        3 + objects_count * 3 + j * 3 + 1
                    ]
                    sol_state[store_count][objects_count * 3 + j * 3 + 2] = row[
                        3 + objects_count * 3 + j * 3 + 2
                    ]

                count += 1
                store_count += 1
                if not no_progress_bar:
                    progress_bar.update(task, completed=store_count)

                # Extending memory buffer
                if store_count % memory_buffer_size == 0:
                    sol_time = np.concatenate((sol_time, np.zeros(memory_buffer_size)))
                    sol_dt = np.concatenate((sol_dt, np.zeros(memory_buffer_size)))
                    energy = np.concatenate((energy, np.zeros(memory_buffer_size)))
                    sol_state = np.concatenate(
                        (
                            sol_state,
                            np.zeros((memory_buffer_size, objects_count * 3 * 2)),
                        )
                    )

            sol_time = sol_time[:store_count]
            sol_dt = sol_dt[:store_count]
            energy = energy[:store_count]
            sol_state = sol_state[:store_count]
            if len(sol_time) > 0:
                tf = sol_time[-1]
            else:
                warnings.warn("No data is found.")

        if not no_print:
            print("Reading completed.\n")
            print(f"Read data size: {store_count}")
            print()

    return (
        sol_state,
        sol_time,
        sol_dt,
        energy,
        system_name,
        integrator,
        objects_count,
        G,
        tf,
        dt,
        tolerance,
        store_every_n,
        run_time,
        m,
        store_count,
    )
