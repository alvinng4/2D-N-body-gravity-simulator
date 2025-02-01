"""
Simulator module for gravity simulator.

Author: Ching Yin Ng
"""

import copy
import ctypes
import threading
import time
import timeit
from pathlib import Path
from queue import Queue
from typing import Optional

import numpy as np

from . import utils
from .gravitational_system import GravitationalSystem


class Simulator:
    DAYS_PER_YEAR = 365.242189
    AVAILABLE_ACCELERATION_METHODS = ["pairwise", "massless", "barnes_hut"]
    AVAILABLE_STORING_METHODS = ["default", "flush", "disabled"]
    AVAILABLE_INTEGRATORS = [
        "euler",
        "euler_cromer",
        "rk4",
        "leapfrog",
        "rkf45",
        "dopri",
        "dverk",
        "rkf78",
        "ias15",
        "whfast",
    ]
    FIXED_STEP_SIZE_INTEGRATORS = ["euler", "euler_cromer", "rk4", "leapfrog", "whfast"]
    ADAPTIVE_STEP_SIZE_INTEGRATORS = ["rkf45", "dopri", "dverk", "rkf78", "ias15"]
    # Recommended settings for built-in systems with IAS15 integrator
    RECOMMENDED_SETTINGS_BUILT_IN_SYSTEMS = {
        # "template": ["tf", "tf unit", "tolerance", "storing_freq"],
        "circular_binary_orbit": [50.0, "days", 1e-9, 1],
        "eccentric_binary_orbit": [2.6, "years", 1e-9, 1],
        "3d_helix": [20.0, "days", 1e-9, 1],
        "sun_earth_moon": [1.0, "years", 1e-9, 1],
        "figure-8": [20.0, "days", 1e-9, 1],
        "pyth-3-body": [70.0, "days", 1e-9, 1],
        "solar_system": [200.0, "years", 1e-9, 1],
        "solar_system_plus": [250.0, "years", 1e-9, 1],
    }

    def __init__(self, c_lib: ctypes.CDLL) -> None:
        self.c_lib = c_lib

    def launch_simulation(
        self,
        gravitational_system: GravitationalSystem,
        integrator_params: dict,
        acceleration_params: dict,
        storing_params: dict,
        settings: dict,
        tf: float,
        is_exit_ctypes_bool: ctypes.c_bool,
    ) -> None:
        """Launch simulation

        Parameters
        ----------
        gravitational_system : GravitationalSystem
        integrator_params : dict
        acceleration_params : dict
        storing_params : dict
        settings : dict
        tf : float
        is_exit_ctypes_bool : ctypes.c_bool
            Flag to indicate if the simulation should be terminated

        Notes
        -----
        - This function would not check the validity of the input parameters.
        """
        if settings["make_copy_params"]:
            integrator_params = integrator_params.copy()
            acceleration_params = acceleration_params.copy()
            storing_params = storing_params.copy()
            settings = settings.copy()

        if settings["make_copy_system"]:
            gravitational_system = copy.deepcopy(gravitational_system)

        self.gravitational_system = gravitational_system
        self.integrator_params = integrator_params
        self.acceleration_params = acceleration_params
        self.storing_params = storing_params
        self.settings = settings

        if "flush_path" in storing_params:
            Path(storing_params["flush_path"]).parent.mkdir(parents=True, exist_ok=True)
            if Path(storing_params["flush_path"]).is_file():
                Path(storing_params["flush_path"]).unlink()
            flush_path_ctypes = storing_params["flush_path"].encode("utf-8")
        else:
            flush_path_ctypes = None

        sol_state_ctypes = ctypes.POINTER(ctypes.c_double)()
        sol_time_ctypes = ctypes.POINTER(ctypes.c_double)()
        sol_dt_ctypes = ctypes.POINTER(ctypes.c_double)()
        sol_size_ctypes = ctypes.c_int64()
        t_ctypes = ctypes.c_double()
        simulation_last_dt_ctypes = ctypes.c_double()
        run_time_ctypes = ctypes.c_double()

        if not settings["disable_progress_bar"]:
            progress_bar_thread = threading.Thread(
                target=utils.progress_bar_c_lib_simulation,
                args=(tf, t_ctypes, sol_size_ctypes, is_exit_ctypes_bool),
            )

        def simulation_wrapper(c_lib_launch_simulation_python, return_queue, *args):
            return_queue.put(c_lib_launch_simulation_python(*args))

        queue: Queue = Queue()
        simulation_thread = threading.Thread(
            target=simulation_wrapper,
            args=(
                self.c_lib.launch_simulation_python,
                queue,
                gravitational_system.x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                gravitational_system.v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                gravitational_system.m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.c_int(gravitational_system.objects_count),
                ctypes.c_double(gravitational_system.G),
                integrator_params["integrator"].encode("utf-8"),
                ctypes.c_double(integrator_params["dt"]),
                ctypes.c_double(integrator_params["tolerance"]),
                ctypes.c_double(integrator_params["initial_dt"]),
                ctypes.c_double(integrator_params["whfast_kepler_tol"]),
                ctypes.c_int(integrator_params["whfast_kepler_max_iter"]),
                ctypes.c_bool(integrator_params["whfast_kepler_auto_remove"]),
                ctypes.c_double(integrator_params["whfast_kepler_auto_remove_tol"]),
                acceleration_params["method"].encode("utf-8"),
                ctypes.c_double(acceleration_params["opening_angle"]),
                ctypes.c_double(acceleration_params["softening_length"]),
                ctypes.c_int(acceleration_params["order"]),
                storing_params["method"].encode("utf-8"),
                flush_path_ctypes,
                ctypes.c_int(storing_params["storing_freq"]),
                ctypes.byref(sol_state_ctypes),
                ctypes.byref(sol_time_ctypes),
                ctypes.byref(sol_dt_ctypes),
                ctypes.byref(sol_size_ctypes),
                ctypes.byref(t_ctypes),
                ctypes.byref(simulation_last_dt_ctypes),
                ctypes.byref(run_time_ctypes),
                ctypes.c_int(settings["verbose"]),
                ctypes.byref(is_exit_ctypes_bool),
                ctypes.c_double(tf),
            ),
        )

        ### Begin simulation ###
        if settings["verbose"] > 1:
            print("Simulation in progress...")
        simulation_thread.start()
        if not settings["disable_progress_bar"]:
            progress_bar_thread.start()

        while simulation_thread.is_alive():
            time.sleep(0.05)

        simulation_thread.join()
        if not settings["disable_progress_bar"]:
            t_ctypes.value = tf
            progress_bar_thread.join()
        return_code = queue.get()

        if return_code != 0:
            raise RuntimeError("Simulation failed.")

        ### End simulation ###
        self.run_time_ = run_time_ctypes.value
        if settings["verbose"] > 1:
            print(f"Simulation completed! Run time: {self.run_time_:.3f} s")
            print()

        if storing_params["method"] == "default":
            self.data_size_ = sol_size_ctypes.value

            self.sol_state_ = np.ctypeslib.as_array(
                sol_state_ctypes,
                shape=(self.data_size_, self.gravitational_system.objects_count * 6),
            ).copy()
            self.sol_time_ = np.ctypeslib.as_array(
                sol_time_ctypes, shape=(self.data_size_,)
            ).copy()
            self.sol_dt_ = np.ctypeslib.as_array(
                sol_dt_ctypes, shape=(self.data_size_,)
            ).copy()

            self.c_lib.free_memory_real(sol_state_ctypes)
            self.c_lib.free_memory_real(sol_time_ctypes)
            self.c_lib.free_memory_real(sol_dt_ctypes)

    def compute_energy(
        self,
        objects_count: int,
        m: np.ndarray,
        G: float,
        sol_state: np.ndarray,
        is_exit_ctypes_bool: Optional[ctypes.c_bool] = None,
    ) -> np.ndarray:
        """Compute energy of the system

        Parameters
        ----------
        objects_count : int
            Number of objects in the system
        m : np.ndarray
            Masses of the objects
        G : float
            Gravitational constant
        sol_state : np.ndarray
            Solution state of the system
        is_exit_ctypes_bool : ctypes.c_bool, optional
            Flag to indicate if the function should be terminated, by default None
        Returns
        -------
        np.ndarray
            Energy of the system at each time step
        """
        if is_exit_ctypes_bool is None:
            is_exit_ctypes_bool = ctypes.c_bool(False)

        print("Computing energy...")
        npts = len(sol_state)
        energy = np.zeros(npts)

        start = timeit.default_timer()
        count = ctypes.c_int(0)
        progress_bar_thread = threading.Thread(
            target=utils.progress_bar_c_lib_function,
            args=(npts, count, is_exit_ctypes_bool),
        )
        progress_bar_thread.start()
        compute_energy_thread = threading.Thread(
            target=self.c_lib.compute_energy_python,
            args=(
                ctypes.c_int(objects_count),
                m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.c_double(G),
                ctypes.c_int(npts),
                ctypes.byref(count),
                energy.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                sol_state.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.byref(is_exit_ctypes_bool),
            ),
        )
        compute_energy_thread.start()
        compute_energy_thread.join()

        # Close progress bar thread
        count.value = npts
        progress_bar_thread.join()

        stop = timeit.default_timer()
        print(f"Run time: {(stop - start):.3f} s")
        print("")

        return energy

    def compute_linear_momentum(
        self,
        objects_count: int,
        m: np.ndarray,
        sol_state: np.ndarray,
        is_exit_ctypes_bool: Optional[ctypes.c_bool] = None,
    ) -> np.ndarray:
        """Compute linear momentum of the system

        Parameters
        ----------
        objects_count : int
            Number of objects in the system
        m : np.ndarray
            Masses of the objects
        sol_state : np.ndarray
            Solution state of the system
        is_exit_ctypes_bool : ctypes.c_bool, optional
            Flag to indicate if the function should be terminated, by default None
        Returns
        -------
        np.ndarray
            Linear momentum of the system at each time step
        """
        if is_exit_ctypes_bool is None:
            is_exit_ctypes_bool = ctypes.c_bool(False)

        print("Computing linear momentum...")
        npts = len(sol_state)
        linear_momentum = np.zeros(npts)

        start = timeit.default_timer()
        count = ctypes.c_int(0)
        progress_bar_thread = threading.Thread(
            target=utils.progress_bar_c_lib_function,
            args=(npts, count, is_exit_ctypes_bool),
        )
        progress_bar_thread.start()
        compute_linear_momentum_thread = threading.Thread(
            target=self.c_lib.compute_linear_momentum_python,
            args=(
                ctypes.c_int(objects_count),
                m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.c_int(npts),
                ctypes.byref(count),
                linear_momentum.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                sol_state.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.byref(is_exit_ctypes_bool),
            ),
        )
        compute_linear_momentum_thread.start()
        compute_linear_momentum_thread.join()

        # Close progress bar thread
        count.value = npts
        progress_bar_thread.join()

        stop = timeit.default_timer()
        print(f"Run time: {(stop - start):.3f} s")
        print("")

        return linear_momentum

    def compute_angular_momentum(
        self,
        objects_count: int,
        m: np.ndarray,
        sol_state: np.ndarray,
        is_exit_ctypes_bool: Optional[ctypes.c_bool] = None,
    ) -> np.ndarray:
        """Compute angular momentum of the system

        Parameters
        ----------
        objects_count : int
            Number of objects in the system
        m : np.ndarray
            Masses of the objects
        sol_state : np.ndarray
            Solution state of the system
        is_exit_ctypes_bool : ctypes.c_bool, optional
            Flag to indicate if the function should be terminated, by default None
        Returns
        -------
        np.ndarray
            Angular momentum of the system at each time step
        """
        if is_exit_ctypes_bool is None:
            is_exit_ctypes_bool = ctypes.c_bool(False)

        print("Computing angular momentum...")
        npts = len(sol_state)
        angular_momentum = np.zeros(npts)

        start = timeit.default_timer()
        count = ctypes.c_int(0)
        progress_bar_thread = threading.Thread(
            target=utils.progress_bar_c_lib_function,
            args=(npts, count, is_exit_ctypes_bool),
        )
        progress_bar_thread.start()
        compute_angular_momentum_thread = threading.Thread(
            target=self.c_lib.compute_angular_momentum_python,
            args=(
                ctypes.c_int(objects_count),
                m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.c_int(npts),
                ctypes.byref(count),
                angular_momentum.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                sol_state.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.byref(is_exit_ctypes_bool),
            ),
        )
        compute_angular_momentum_thread.start()
        compute_angular_momentum_thread.join()

        # Close progress bar thread
        count.value = npts
        progress_bar_thread.join()

        stop = timeit.default_timer()
        print(f"Run time: {(stop - start):.3f} s")
        print("")

        return angular_momentum

    @staticmethod
    def compute_eccentricity(
        objects_count: int,
        m: np.ndarray,
        G: float,
        sol_state: np.ndarray,
    ) -> np.ndarray:
        """Compute the eccentricity using the sol_state array,
        assuming that the first object is the central object

        Parameters
        ----------
        objects_count : int
            Number of objects in the system
        m : np.ndarray
            Masses of the objects
        G : float
            Gravitational constant
        sol_state : np.ndarray
            Solution state of the system

        Returns
        -------
        np.ndarray
            Eccentricity of the system at each time step

        Notes
        -----
        - The function assumes that the first object is the central object.
        - C library function is not used here since this can be done with
          purely numpy vectorized operations. Nevertheless, we may consider
          implementing this in C library in the future.
        """
        print("Computing eccentricity (Assuming the first body is the central star)...")
        eccentricity = np.zeros(len(sol_state))

        start = timeit.default_timer()
        x = (
            sol_state[:, 3 : (objects_count * 3)]
            .reshape(-1, (objects_count - 1), 3)
            .copy()
        )
        v = (
            sol_state[:, (objects_count + 1) * 3 :]
            .reshape(-1, (objects_count - 1), 3)
            .copy()
        )

        x = x - sol_state[:, :3].reshape(-1, 1, 3)
        v = v - sol_state[:, (objects_count) * 3 : (objects_count + 1) * 3].reshape(
            -1, 1, 3
        )

        eccentricity = (
            np.cross(v, np.cross(x, v)) / (G * (m[0] + m[1:]))[:, np.newaxis]
            - x / np.linalg.norm(x, axis=2)[:, :, np.newaxis]
        )
        eccentricity = np.linalg.norm(eccentricity, axis=2)

        stop = timeit.default_timer()
        print(f"Run time: {(stop - start):.3f} s")
        print("")

        return eccentricity

    @staticmethod
    def compute_inclination(
        objects_count: int,
        sol_state: np.ndarray,
    ) -> np.ndarray:
        """Compute the inclination using the sol_state array,
        assuming that the first object is the central object

        Parameters
        ----------
        objects_count : int
            Number of objects in the system
        m : np.ndarray
            Masses of the objects
        G : float
            Gravitational constant
        sol_state : np.ndarray
            Solution state of the system

        Returns
        -------
        np.ndarray
            Inclination of the system at each time step

        Notes
        -----
        - The function assumes that the first object is the central object.
        - C library function is not used here since this can be done with
            purely numpy vectorized operations. Nevertheless, we may consider
            implementing this in C library in the future.
        """
        print("Computing inclination (Assuming the first body is the central star)...")
        inclination = np.zeros(len(sol_state))

        start = timeit.default_timer()
        x = (
            sol_state[:, 3 : (objects_count * 3)]
            .reshape(-1, (objects_count - 1), 3)
            .copy()
        )
        v = (
            sol_state[:, (objects_count + 1) * 3 :]
            .reshape(-1, (objects_count - 1), 3)
            .copy()
        )

        x = x - sol_state[:, :3].reshape(-1, 1, 3)
        v = v - sol_state[:, (objects_count) * 3 : (objects_count + 1) * 3].reshape(
            -1, 1, 3
        )

        unit_angular_momentum_vector = (
            np.cross(x, v) / np.linalg.norm(np.cross(x, v), axis=2)[:, :, np.newaxis]
        )
        unit_z = np.array([0, 0, 1])

        inclination = np.arccos(np.sum(unit_angular_momentum_vector * unit_z, axis=2))

        stop = timeit.default_timer()
        print(f"Run time: {(stop - start):.3f} s")
        print("")

        return inclination
