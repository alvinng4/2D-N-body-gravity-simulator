"""
Utility functions for gravity simulator

Author: Ching Yin Ng
"""

import csv
import ctypes
import platform
import sys
import time
import timeit
from pathlib import Path
from typing import Optional

import numpy as np
import rich.progress

# Increase field limit for reading CSV files
new_field_lim = sys.maxsize
while True:
    try:
        csv.field_size_limit(new_field_lim)
        break
    except OverflowError:
        new_field_lim = new_field_lim // 10


class Progress_bar(rich.progress.Progress):
    def __init__(self):
        super().__init__(
            rich.progress.BarColumn(),
            rich.progress.TextColumn("[green]{task.percentage:>3.0f}%"),
            rich.progress.TextColumn("•"),
            rich.progress.TimeElapsedColumn(),
            rich.progress.TextColumn("•"),
            rich.progress.TimeRemainingColumn(),
        )


class Progress_bar_with_data_size(rich.progress.Progress):
    def __init__(self):
        super().__init__(
            rich.progress.BarColumn(),
            rich.progress.TextColumn("[green]{task.percentage:>3.0f}%"),
            rich.progress.TextColumn("•"),
            rich.progress.TimeElapsedColumn(),
            rich.progress.TextColumn("•"),
            rich.progress.TimeRemainingColumn(),
            "• [magenta]Data size: {task.fields[store_count]}",
        )


def progress_bar_c_lib_function(
    npts: int,
    count: ctypes.c_int | ctypes.c_int64,
    is_exit_ctypes_bool: ctypes.c_bool,
):
    """Progress bar for function with C library

    Parameters
    ----------
    npts : int
        Total number of points
    count : ctypes.c_int | ctypes.c_int64
        Current count
    is_exit_ctypes_bool : ctypes.c_bool
        Flag to indicate if the function should be terminated
    """

    progress_bar = Progress_bar()
    with progress_bar:
        task = progress_bar.add_task("", total=npts, count=count.value)

        while not is_exit_ctypes_bool.value:
            # Update progress bar
            progress_bar.update(task, completed=count.value)
            time.sleep(0.1)

            if count.value >= npts:
                break

        progress_bar.update(task, completed=npts)


def progress_bar_c_lib_simulation(
    tf: float,
    t: ctypes.c_double,
    store_count: ctypes.c_int64,
    is_exit_ctypes_bool: ctypes.c_bool,
) -> None:
    """Progress bar for simulation with C library

    Parameters
    ----------
    tf : float
        Total simulation time
    t : ctypes.c_double
        Current simulation time
    store_count : ctypes.c_int64
        Number of data points stored
    is_exit_ctypes_bool : ctypes.c_bool
        Flag to indicate if the simulation should be terminated
    """

    progress_bar = Progress_bar_with_data_size()
    with progress_bar:
        task = progress_bar.add_task("", total=tf, store_count=store_count.value)

        while not is_exit_ctypes_bool.value:
            # Update progress bar
            progress_bar.update(task, completed=t.value, store_count=store_count.value)
            time.sleep(0.1)

            if t.value >= tf:
                break

        progress_bar.update(task, completed=tf, store_count=store_count.value)


def load_c_lib(c_lib_path: Optional[Path] = None) -> ctypes.CDLL:
    """Load the C dynamic-link library

    Returns
    -------
    c_lib : ctypes.CDLL
        C dynamic-link library object

    Raises
    ------
    OSError
        If the platform is not supported
    FileNotFoundError
        If the C library is not found at the path
    """
    if c_lib_path is None:
        if platform.system() == "Windows":
            c_lib_path = Path(__file__).parent.parent / "src" / "c_lib.dll"
        elif platform.system() == "Darwin":
            c_lib_path = Path(__file__).parent.parent / "src" / "c_lib.dylib"
        elif platform.system() == "Linux":
            c_lib_path = Path(__file__).parent.parent / "src" / "c_lib.so"
        else:
            raise Exception(
                f'Platform "{platform.system()}" not supported. Supported platforms'
                + ": Windows, macOS, Linux."
                + "You may bypass this error by providing the path to the C library"
            )

    if not Path(c_lib_path).exists():
        raise FileNotFoundError(f'C library not found at path: "{c_lib_path}"')

    return ctypes.cdll.LoadLibrary(str(c_lib_path))


def initialize_c_lib(c_lib: ctypes.CDLL) -> None:
    """Initialize C library

    Parameters
    ----------
    c_lib : ctypes.CDLL
        C dynamic-link library object

    Raises
    ------
    ValueError
        If any C functions that should be available are not found in the C library
    """
    for function in ["launch_simulation"]:
        if not hasattr(c_lib, function):
            raise ValueError(f"Function {function} not found in the C library.")

    c_lib.launch_simulation_python.restype = ctypes.c_int


def trim_data(
    data: int | np.ndarray,
    trim_freq: int,
) -> int | np.ndarray:
    """Trim data to reduce data size

    Parameters
    ----------
    data : int | np.ndarray
        Data to be trimmed
    trim_freq : int
        Frequency to trim data

    Returns
    -------
    trimmed_data : int | np.ndarray
        Trimmed data
    """
    if isinstance(data, int):
        return data // trim_freq
    elif isinstance(data, np.ndarray):
        return data[::trim_freq]
    else:
        raise TypeError("Data type not supported.")


def save_results_csv(
    file_path: str | Path,
    sol_state_: np.ndarray,
    sol_time_: np.ndarray,
    sol_dt_: np.ndarray,
    sol_energy_: np.ndarray,
    disable_progress_bar: bool = False,
) -> None:
    """Save simulation results to a CSV file

    Notes
    -----
    Unit: Solar masses, AU, day
    Format: time, dt, total energy, x1, y1, z1, x2, y2, z2, ... vx1, vy1, vz1, vx2, vy2, vz2, ...

    Parameters
    ----------
    file_path : str | Path
        Path to save the CSV file
    sol_state_ : np.ndarray
        State of the system at each time step
    sol_time_ : np.ndarray
        Time at each time step
    sol_dt_ : np.ndarray
        Time step at each time step
    sol_energy_ : np.ndarray
        Energy of the system at each time step
    disable_progress_bar : bool, optional
        Disable progress bar, by default False
    """
    if not isinstance(file_path, Path):
        file_path = Path(file_path)

    data_size = len(sol_state_)
    if not disable_progress_bar:
        print("Saving results to CSV file...")
        start = timeit.default_timer()
        progress_bar = Progress_bar()
        with progress_bar:
            with file_path.open("w", newline="") as file:
                writer = csv.writer(file)
                for count in progress_bar.track(range(data_size)):
                    row = np.insert(
                        sol_state_[count],
                        0,
                        sol_energy_[count],
                    )
                    row = np.insert(row, 0, sol_dt_[count])
                    row = np.insert(row, 0, sol_time_[count])
                    writer.writerow(row.tolist())
        end = timeit.default_timer()
        print(f"Run time: {end - start:.2f} s")
    else:
        with file_path.open("w", newline="") as file:
            writer = csv.writer(file)
            for count in range(data_size):
                row = np.insert(
                    sol_state_[count],
                    0,
                    sol_energy_[count],
                )
                row = np.insert(row, 0, sol_dt_[count])
                row = np.insert(row, 0, sol_time_[count])
                writer.writerow(row.tolist())


def read_results_csv(
    file_path: str | Path,
) -> dict[str, np.ndarray]:
    """Read simulation results from a CSV file

    Notes
    -----
    Unit: Solar masses, AU, day
    Format: time, dt, total energy, x1, y1, z1, x2, y2, z2, ... vx1, vy1, vz1, vx2, vy2, vz2, ...

    Parameters
    ----------
    file_path : str | Path
        Path to the CSV file

    Returns
    -------
    results : dict[str, np.ndarray]
        Simulation results, with field names: time, dt, energy, state
    """
    if not isinstance(file_path, Path):
        file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    time = []
    dt = []
    energy = []
    state = []

    with file_path.open("r") as file:
        reader = csv.reader(file)
        for row in reader:
            time.append(float(row[0]))
            dt.append(float(row[1]))
            energy.append(float(row[2]))
            state.append(np.array(row[3:], dtype=float))

    results = {
        "time": np.array(time),
        "dt": np.array(dt),
        "energy": np.array(energy),
        "state": np.array(state),
    }

    return results


def keplerian_to_cartesian(
    semi_major_axis: float,
    eccentricity: float,
    inclination: float,
    argument_of_periapsis: float,
    longitude_of_ascending_node: float,
    true_anomaly: float,
    total_mass: float,
    G: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert keplerian elements to cartesian coordinates

    Parameters
    ----------
    semi_major_axis : float
    eccentricity : float
    inclination : float
    argument_of_periapsis : float
    longitude_of_ascending_node : float
    true_anomaly : float
    total_mass : float
    G : float

    Reference
    ---------
    Moving Planets Around: An Introduction to N-Body
    Simulations Applied to Exoplanetary Systems, Chapter 2
    """

    cos_inc = np.cos(inclination)
    sin_inc = np.sin(inclination)

    cos_arg_periapsis = np.cos(argument_of_periapsis)
    sin_arg_periapsis = np.sin(argument_of_periapsis)

    cos_long_asc_node = np.cos(longitude_of_ascending_node)
    sin_long_asc_node = np.sin(longitude_of_ascending_node)

    cos_true_anomaly = np.cos(true_anomaly)
    sin_true_anomaly = np.sin(true_anomaly)

    # ecc_unit_vec is the unit vector pointing towards periapsis
    ecc_unit_vec = np.zeros(3)
    ecc_unit_vec[0] = (
        cos_long_asc_node * cos_arg_periapsis
        - sin_long_asc_node * sin_arg_periapsis * cos_inc
    )
    ecc_unit_vec[1] = (
        sin_long_asc_node * cos_arg_periapsis
        + cos_long_asc_node * sin_arg_periapsis * cos_inc
    )
    ecc_unit_vec[2] = sin_arg_periapsis * sin_inc

    # q_unit_vec is the unit vector that is perpendicular to ecc_unit_vec and orbital angular momentum vector
    q_unit_vec = np.zeros(3)
    q_unit_vec[0] = (
        -cos_long_asc_node * sin_arg_periapsis
        - sin_long_asc_node * cos_arg_periapsis * cos_inc
    )
    q_unit_vec[1] = (
        -sin_long_asc_node * sin_arg_periapsis
        + cos_long_asc_node * cos_arg_periapsis * cos_inc
    )
    q_unit_vec[2] = cos_arg_periapsis * sin_inc

    # Calculate the position vector
    x = (
        semi_major_axis
        * (1.0 - eccentricity**2)
        / (1.0 + eccentricity * cos_true_anomaly)
        * (cos_true_anomaly * ecc_unit_vec + sin_true_anomaly * q_unit_vec)
    )
    v = np.sqrt(G * total_mass / (semi_major_axis * (1.0 - eccentricity**2))) * (
        -sin_true_anomaly * ecc_unit_vec
        + (eccentricity + cos_true_anomaly) * q_unit_vec
    )

    if np.isnan(x).any() or np.isnan(v).any():
        raise ValueError("Invalid values. Please check the input values.")

    return x, v
