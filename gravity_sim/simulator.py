import csv
import ctypes
import datetime
from pathlib import Path
import threading
import timeit
import warnings

import numpy as np

import common
from integrator_simple import SimpleIntegrator
from integrator_rk_embedded import RKEmbedded
from integrator_ias15 import IAS15
from integrator_whfast import WHFast


class Simulator:
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
    FIXED_STEP_SIZE_INTEGRATORS = ["euler", "euler_cromer", "rk4", "leapfrog", "whfast"]
    ADAPTIVE_STEP_SIZE_INTEGRATORS = ["rkf45", "dopri", "dverk", "rkf78", "ias15"]

    def __init__(self, c_lib=None, is_exit_ctypes_bool=None):
        self.c_lib = c_lib
        if c_lib is not None:
            self.integration_mode = "c_lib"
        else:
            self.integration_mode = "numpy"

        if is_exit_ctypes_bool is None:
            self.is_exit_ctypes_bool = ctypes.c_bool(False)
        else:
            self.is_exit_ctypes_bool = is_exit_ctypes_bool

        self.store_every_n = None
        self.run_time = None

    def launch_simulation(
        self,
        gravitational_system,
        integrator: str,
        tf: float,
        dt: float = None,
        tolerance: float = None,
        store_every_n: int = 1,
        acceleration_method: str = "pairwise",
        storing_method: str = "default",
        flush_path: str = None,
        no_progress_bar: bool = False,
        no_print: bool = False,
        softening_length: float = 0.0,
        barnes_hut_theta: float = 0.5,
        **kwargs,
    ) -> None:
        """
        Launch simulation

        Parameters
        ----------
        gravitational_system : GravitationalSystem
        integrator : str
        tf : float
            Integration time
        dt : float, optional
            Time step
        tolerance : float, optional
            Tolerance for adaptive time step integrators
        store_every_n : int, optional
            Store every nth time step
        acceleration_method : str, optional
            Acceleration method -- "pairwise", "massless"
                pairwise: Pairwise acceleration
                massless: System with massless particles
        storing_method : str, optional
            Storing method -- "default", "flush", "no_store"
                default: Store output into memory
                flush : Store output into data file directly
                no_store : Do not store output
        flush_path : str, optional
            Path to store flushed results
        no_progress_bar : bool, optional
            Disable progress bar
        no_print : bool, optional
            Disable print statements
        softening_length : float, optional
            Softening length for acceleration acceleration
        barnes_hut_theta : float, optional
            Theta parameter for Barnes-Hut algorithm
        kwargs : dict
            Additional keyword arguments
        """
        self.system_name = gravitational_system.name
        self.x = gravitational_system.x.copy()
        self.v = gravitational_system.v.copy()
        self.m = gravitational_system.m.copy()
        self.objects_count = gravitational_system.objects_count
        self.G = gravitational_system.G

        self.integrator = integrator.strip().lower()
        if self.integrator not in self.AVAILABLE_INTEGRATORS:
            raise ValueError(
                f"Invalid integrator. Must be one of {self.AVAILABLE_INTEGRATORS}."
            )
        self.tf = tf
        self.store_every_n = int(store_every_n)
        self.acceleration_method = acceleration_method.strip().lower()
        self.flush_path = flush_path
        self.no_progress_bar = no_progress_bar
        self.no_print = no_print
        self.softening_length = softening_length
        self.barnes_hut_theta = barnes_hut_theta
        self.kwargs = kwargs

        if (dt is None) and (self.integrator in self.FIXED_STEP_SIZE_INTEGRATORS):
            raise ValueError("dt must be provided for fixed step size integrators.")
        elif (tolerance is None) and (self.integrator in self.ADAPTIVE_STEP_SIZE_INTEGRATORS):
            raise ValueError("tolerance must be provided for adaptive step size integrators.")

        if dt is not None:
            self.dt = dt
        else:
            self.dt = 0.0
            
        self.tolerance = tolerance

        if self.acceleration_method == "barnes-hut" and (self.m == 0.0).any():
            warnings.warn(
                "barnes-hut: Massless particles detected, adding m=1e-30 to massless particles."
            )
            self.m[self.m == 0.0] = 1e-30

        match storing_method:
            case "default":
                self._storing_method = 0
            case "flush":
                self._storing_method = 1
            case "no_store":
                self._storing_method = 2

        # Flushing           
        if self._storing_method == 1:
            if self.flush_path is None:
                file_path = Path(__file__).parent / "results"
                file_path.mkdir(parents=True, exist_ok=True)
                flush_path = (
                    file_path
                    / f"temp_data_{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))}.csv"
                )
        else:
            flush_path = None

        if (self.integration_mode == "numpy"):
            if self.acceleration_method != "pairwise":
                warnings.warn(
                    "Only pairwise acceleration is available for NumPy. "
                    + 'Setting acceleration method to "pairwise".'
                )

            if self.softening_length != 0.0:
                warnings.warn(
                    "Softening length is not available for NumPy. "
                    + "Ignoring the provided softening length."
                )

        if not self.no_print:
            print("Simulating the system...")
        start = timeit.default_timer()

        if self.integration_mode == "c_lib":
            match self.integrator:
                case "euler" | "euler_cromer" | "rk4" | "leapfrog":
                    integrator = SimpleIntegrator(
                        self.store_every_n,
                        self.c_lib,
                        self.is_exit_ctypes_bool,
                    )

                    (
                        self.sol_state,
                        self.sol_time,
                        self.sol_dt,
                        store_count,
                    ) = integrator.simulation_c_lib(
                        self.integrator,
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.dt,
                        self.tf,
                        self.acceleration_method,
                        self.softening_length,
                        self.barnes_hut_theta,
                        self._storing_method,
                        str(flush_path),
                        self.no_progress_bar,
                    )

                case "rkf45" | "dopri" | "dverk" | "rkf78":
                    integrator = RKEmbedded(
                        self.integrator,
                        self.store_every_n,
                        self.c_lib,
                        self.is_exit_ctypes_bool,
                    )

                    (
                        self.sol_state,
                        self.sol_time,
                        self.sol_dt,
                        store_count,
                        self.dt,
                    ) = integrator.simulation_c_lib(
                        integrator.order,
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.dt,
                        self.tf,
                        self.tolerance,
                        self.tolerance,
                        self.acceleration_method,
                        self.softening_length,
                        self.barnes_hut_theta,
                        self._storing_method,
                        str(flush_path),
                        self.no_progress_bar,
                    )

                case "ias15":
                    integrator = IAS15(
                        self.store_every_n, self.c_lib, self.is_exit_ctypes_bool
                    )
                    (
                        self.sol_state,
                        self.sol_time,
                        self.sol_dt,
                        store_count,
                        self.dt,
                    ) = integrator.simulation_c_lib(
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.dt,
                        self.tf,
                        self.tolerance,
                        self.acceleration_method,
                        self.softening_length,
                        self.barnes_hut_theta,
                        self._storing_method,
                        str(flush_path),
                        self.no_progress_bar,
                    )
                case "whfast":
                    integrator = WHFast(
                        self.store_every_n,
                        self.c_lib,
                        self.is_exit_ctypes_bool,
                    )

                    (
                        self.sol_state,
                        self.sol_time,
                        self.sol_dt,
                        store_count,
                        self.x,
                        self.v,
                        self.m,
                        self.objects_count,
                    ) = integrator.simulation_c_lib(
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.dt,
                        self.tf,
                        self.acceleration_method,
                        self.softening_length,
                        self.barnes_hut_theta,
                        self._storing_method,
                        str(flush_path),
                        self.no_progress_bar,
                        **self.kwargs,
                    )

        else:
            match self.integrator:
                case "euler" | "euler_cromer" | "rk4" | "leapfrog":
                    integrator = SimpleIntegrator(
                        self.store_every_n,
                        self.c_lib,
                        self.is_exit_ctypes_bool,
                    )

                    (
                        self.sol_state,
                        self.sol_time,
                        self.sol_dt,
                    ) = integrator.simulation_numpy(
                        self.integrator,
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.dt,
                        self.tf,
                        self._storing_method,
                        str(flush_path),
                        self.no_progress_bar,
                    )

                case "rkf45" | "dopri" | "dverk" | "rkf78":
                    integrator = RKEmbedded(
                        self.integrator,
                        self.store_every_n,
                        self.c_lib,
                        self.is_exit_ctypes_bool,
                    )

                    (
                        self.sol_state,
                        self.sol_time,
                        self.sol_dt,
                        self.dt,
                    ) = integrator.simulation_numpy(
                        integrator.order,
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.dt,
                        self.tf,
                        self.tolerance,
                        self.tolerance,
                        self._storing_method,
                        str(flush_path),
                        self.no_progress_bar,
                    )

                case "ias15":
                    integrator = IAS15(
                        self.store_every_n,
                        self.c_lib,
                        self.is_exit_ctypes_bool,
                    )

                    (
                        self.sol_state,
                        self.sol_time,
                        self.sol_dt,
                        self.dt,
                    ) = integrator.simulation_numpy(
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.dt,
                        self.tf,
                        self.tolerance,
                        self._storing_method,
                        str(flush_path),
                        self.no_progress_bar,
                    )
                case "whfast":
                    integrator = WHFast(
                        self.store_every_n,
                        self.c_lib,
                        self.is_exit_ctypes_bool,
                    )

                    (
                        self.sol_state,
                        self.sol_time,
                        self.sol_dt,
                    ) = integrator.simulation_numpy(
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.dt,
                        self.tf,
                        self._storing_method,
                        str(flush_path),
                        self.no_progress_bar,
                    )

        stop = timeit.default_timer()
        self.run_time = stop - start

        if not self.no_print:
            print(f"Run time: {self.run_time:.3f} s")
            print("")

        # Flush
        if self._storing_method == 1:
            if self.flush_path is None:
                results_path = (
                    Path(__file__).parent
                    / "results"
                    / (
                        str(datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
                        + "_result.csv"
                    )
                )
            else:
                results_path = Path(flush_path)

            try:
                integrator_name = self.AVAILABLE_INTEGRATORS_TO_PRINTABLE_NAMES[
                    self.integrator
                ]
            except KeyError:
                integrator_name = None

            self.combine_metadata_with_flushed_data(
                results_file_path=str(results_path),
                flushed_file_path=str(flush_path),
                system_name=self.system_name,
                integrator_name=integrator_name,
                objects_count=self.objects_count,
                G=self.G,
                tf=self.tf,
                dt=self.dt,
                tolerance=self.tolerance,
                data_size=store_count,
                store_every_n=self.store_every_n,
                run_time=self.run_time,
                masses=self.m,
                no_print=self.no_print,
            )
            flush_path.unlink()

    def resume_simulation(self, tf: float) -> None:
        """
        Resume simulation

        Restart the simulation from the last saved state

        Parameters
        ----------
        tf : float
            Integration time
        """
        self.tf = tf
        if self.integrator == "whfast":
            warnings.warn(
                "resume_simulation() is unstable for WHFast, and are likely to " + 
                "induce errors."
            )

        # Flushing           
        if self._storing_method == 1:
            if self.flush_path is None:
                file_path = Path(__file__).parent / "results"
                file_path.mkdir(parents=True, exist_ok=True)
                flush_path = (
                    file_path
                    / f"temp_data_{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))}.csv"
                )
        else:
            flush_path = None

        if not self.no_print:
            print("Simulating the system...")
        start = timeit.default_timer()

        if self.integration_mode == "c_lib":
            match self.integrator:
                case "euler" | "euler_cromer" | "rk4" | "leapfrog":
                    integrator = SimpleIntegrator(
                        self.store_every_n,
                        self.c_lib,
                        self.is_exit_ctypes_bool,
                    )

                    (
                        self.sol_state,
                        self.sol_time,
                        self.sol_dt,
                        store_count,
                    ) = integrator.simulation_c_lib(
                        self.integrator,
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.dt,
                        self.tf,
                        self.acceleration_method,
                        self.softening_length,
                        self.barnes_hut_theta,
                        self._storing_method,
                        str(flush_path),
                        self.no_progress_bar,
                    )

                case "rkf45" | "dopri" | "dverk" | "rkf78":
                    integrator = RKEmbedded(
                        self.integrator,
                        self.store_every_n,
                        self.c_lib,
                        self.is_exit_ctypes_bool,
                    )

                    (
                        self.sol_state,
                        self.sol_time,
                        self.sol_dt,
                        store_count,
                        self.dt,
                    ) = integrator.simulation_c_lib(
                        integrator.order,
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.dt,
                        self.tf,
                        self.tolerance,
                        self.tolerance,
                        self.acceleration_method,
                        self.softening_length,
                        self.barnes_hut_theta,
                        self._storing_method,
                        str(flush_path),
                        self.no_progress_bar,
                    )

                case "ias15":
                    integrator = IAS15(
                        self.store_every_n, self.c_lib, self.is_exit_ctypes_bool
                    )
                    (
                        self.sol_state,
                        self.sol_time,
                        self.sol_dt,
                        store_count,
                        self.dt,
                    ) = integrator.simulation_c_lib(
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.dt,
                        self.tf,
                        self.tolerance,
                        self.acceleration_method,
                        self.softening_length,
                        self.barnes_hut_theta,
                        self._storing_method,
                        str(flush_path),
                        self.no_progress_bar,
                    )
                case "whfast":
                    integrator = WHFast(
                        self.store_every_n,
                        self.c_lib,
                        self.is_exit_ctypes_bool,
                    )

                    (
                        self.sol_state,
                        self.sol_time,
                        self.sol_dt,
                        store_count,
                        self.x,
                        self.v,
                        self.m,
                        self.objects_count,
                    ) = integrator.simulation_c_lib(
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.dt,
                        self.tf,
                        self.acceleration_method,
                        self.softening_length,
                        self.barnes_hut_theta,
                        self._storing_method,
                        str(flush_path),
                        self.no_progress_bar,
                        **self.kwargs,
                    )

        else:
            match self.integrator:
                case "euler" | "euler_cromer" | "rk4" | "leapfrog":
                    integrator = SimpleIntegrator(
                        self.store_every_n,
                        self.c_lib,
                        self.is_exit_ctypes_bool,
                    )

                    (
                        self.sol_state,
                        self.sol_time,
                        self.sol_dt,
                    ) = integrator.simulation_numpy(
                        self.integrator,
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.dt,
                        self.tf,
                        self._storing_method,
                        str(flush_path),
                        self.no_progress_bar,
                    )

                case "rkf45" | "dopri" | "dverk" | "rkf78":
                    integrator = RKEmbedded(
                        self.integrator,
                        self.store_every_n,
                        self.c_lib,
                        self.is_exit_ctypes_bool,
                    )

                    (
                        self.sol_state,
                        self.sol_time,
                        self.sol_dt,
                        self.dt,
                    ) = integrator.simulation_numpy(
                        integrator.order,
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.dt,
                        self.tf,
                        self.tolerance,
                        self.tolerance,
                        self._storing_method,
                        str(flush_path),
                        self.no_progress_bar,
                    )

                case "ias15":
                    integrator = IAS15(
                        self.store_every_n,
                        self.c_lib,
                        self.is_exit_ctypes_bool,
                    )

                    (
                        self.sol_state,
                        self.sol_time,
                        self.sol_dt,
                        self.dt,
                    ) = integrator.simulation_numpy(
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.dt,
                        self.tf,
                        self.tolerance,
                        self._storing_method,
                        str(flush_path),
                        self.no_progress_bar,
                    )
                case "whfast":
                    integrator = WHFast(
                        self.store_every_n,
                        self.c_lib,
                        self.is_exit_ctypes_bool,
                    )

                    (
                        self.sol_state,
                        self.sol_time,
                        self.sol_dt,
                    ) = integrator.simulation_numpy(
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.dt,
                        self.tf,
                        self._storing_method,
                        str(flush_path),
                        self.no_progress_bar,
                    )

        stop = timeit.default_timer()
        self.run_time = stop - start

        if not self.no_print:
            print(f"Run time: {self.run_time:.3f} s")
            print("")

        # Flush
        if self._storing_method == 1:
            if self.flush_path is None:
                results_path = (
                    Path(__file__).parent
                    / "results"
                    / (
                        str(datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
                        + "_result.csv"
                    )
                )
            else:
                results_path = Path(flush_path)

            try:
                integrator_name = self.AVAILABLE_INTEGRATORS_TO_PRINTABLE_NAMES[
                    self.integrator
                ]
            except KeyError:
                integrator_name = None

            self.combine_metadata_with_flushed_data(
                results_file_path=str(results_path),
                flushed_file_path=str(flush_path),
                system_name=self.system_name,
                integrator_name=integrator_name,
                objects_count=self.objects_count,
                G=self.G,
                tf=self.tf,
                dt=self.dt,
                tolerance=self.tolerance,
                data_size=store_count,
                store_every_n=self.store_every_n,
                run_time=self.run_time,
                masses=self.m,
                no_print=self.no_print,
            )
            flush_path.unlink()

    def compute_energy(self):
        """
        Compute the total energy using the sol_state array
        """
        print("Computing energy...")
        npts = len(self.sol_state)
        self.energy = np.zeros(npts)

        start = timeit.default_timer()
        if self.integration_mode == "c_lib":
            count = ctypes.c_int(0)
            progress_bar_thread = threading.Thread(
                target=common.progress_bar_c_lib_fixed_step_size,
                args=(npts, count, self.is_exit_ctypes_bool),
            )
            progress_bar_thread.start()
            compute_energy_thread = threading.Thread(
                target=self.c_lib.compute_energy,
                args=(
                    ctypes.c_int(self.objects_count),
                    ctypes.c_int(npts),
                    ctypes.byref(count),
                    self.energy.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    self.sol_state.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    self.m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    ctypes.c_double(self.G),
                    ctypes.byref(self.is_exit_ctypes_bool),
                ),
            )
            compute_energy_thread.start()
            compute_energy_thread.join()

            # Close progress bar thread
            count.value = npts
            progress_bar_thread.join()

        else:
            progress_bar = common.Progress_bar()
            with progress_bar:
                for count in progress_bar.track(range(npts), description=""):
                    x = self.sol_state[count]
                    for i in range(self.objects_count):
                        # KE
                        self.energy[count] += (
                            0.5
                            * self.m[i]
                            * np.linalg.norm(
                                x[
                                    (self.objects_count + i)
                                    * 3 : (self.objects_count + 1 + i)
                                    * 3
                                ]
                            )
                            ** 2
                        )
                        # PE
                        for j in range(i + 1, self.objects_count):
                            self.energy[count] -= (
                                self.G
                                * self.m[i]
                                * self.m[j]
                                / np.linalg.norm(
                                    x[i * 3 : (i + 1) * 3] - x[j * 3 : (j + 1) * 3]
                                )
                            )

        stop = timeit.default_timer()
        print(f"Run time: {(stop - start):.3f} s")
        print("")

    def compute_linear_momentum(self):
        """
        Compute the total linear momentum using the sol_state array
        """
        print("Computing linear momentum...")
        npts = len(self.sol_state)
        self.linear_momentum = np.zeros(npts)

        start = timeit.default_timer()
        if self.integration_mode == "c_lib":
            count = ctypes.c_int(0)
            progress_bar_thread = threading.Thread(
                target=common.progress_bar_c_lib_fixed_step_size,
                args=(npts, count, self.is_exit_ctypes_bool),
            )
            progress_bar_thread.start()
            compute_linear_momentum_thread = threading.Thread(
                target=self.c_lib.compute_linear_momentum,
                args=(
                    ctypes.c_int(self.objects_count),
                    ctypes.c_int(npts),
                    ctypes.byref(count),
                    self.linear_momentum.ctypes.data_as(
                        ctypes.POINTER(ctypes.c_double)
                    ),
                    self.sol_state.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    self.m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    ctypes.byref(self.is_exit_ctypes_bool),
                ),
            )
            compute_linear_momentum_thread.start()
            compute_linear_momentum_thread.join()

            # Close progress bar thread
            count.value = npts
            progress_bar_thread.join()

        else:
            progress_bar = common.Progress_bar()
            with progress_bar:
                for count in progress_bar.track(range(npts), description=""):
                    temp_vec = np.zeros(3)
                    for i in range(self.objects_count):
                        temp_vec += (
                            self.m[i]
                            * self.sol_state[count][
                                (self.objects_count + i)
                                * 3 : (self.objects_count + 1 + i)
                                * 3
                            ]
                        )

                    self.linear_momentum[count] = np.linalg.norm(temp_vec)

        stop = timeit.default_timer()
        print(f"Run time: {(stop - start):.3f} s")
        print("")

    def compute_angular_momentum(self):
        """
        Compute the total angular momentum using the sol_state array
        """
        print("Computing angular momentum...")
        npts = len(self.sol_state)
        self.angular_momentum = np.zeros(npts)

        start = timeit.default_timer()
        if self.integration_mode == "c_lib":
            count = ctypes.c_int(0)
            progress_bar_thread = threading.Thread(
                target=common.progress_bar_c_lib_fixed_step_size,
                args=(npts, count, self.is_exit_ctypes_bool),
            )
            progress_bar_thread.start()
            compute_angular_momentum_thread = threading.Thread(
                target=self.c_lib.compute_angular_momentum,
                args=(
                    ctypes.c_int(self.objects_count),
                    ctypes.c_int(npts),
                    ctypes.byref(count),
                    self.angular_momentum.ctypes.data_as(
                        ctypes.POINTER(ctypes.c_double)
                    ),
                    self.sol_state.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    self.m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    ctypes.byref(self.is_exit_ctypes_bool),
                ),
            )
            compute_angular_momentum_thread.start()
            compute_angular_momentum_thread.join()

            # Close progress bar thread
            count.value = npts
            progress_bar_thread.join()

        else:
            progress_bar = common.Progress_bar()
            with progress_bar:
                for count in progress_bar.track(range(npts), description=""):
                    temp_vec = np.zeros(3)
                    for i in range(self.objects_count):
                        temp_vec += self.m[i] * np.cross(
                            self.sol_state[count][(i * 3) : (i + 1) * 3],
                            self.sol_state[count][
                                (self.objects_count + i)
                                * 3 : (self.objects_count + 1 + i)
                                * 3
                            ],
                        )
                    self.angular_momentum[count] = np.linalg.norm(temp_vec)

        stop = timeit.default_timer()
        print(f"Run time: {(stop - start):.3f} s")
        print("")

    def compute_eccentricity(self):
        """
        Compute the eccentricity using the sol_state array,
        assuming that the first object is the central object

        C library is not implemented since this is mostly NumPy
        vector operations, which can be done almost instantly.
        """
        print("Computing eccentricity (Assuming the first body is the central star)...")
        self.eccentricity = np.zeros(self.data_size)

        start = timeit.default_timer()
        x = (
            self.sol_state[:, 3 : (self.objects_count * 3)]
            .reshape(-1, (self.objects_count - 1), 3)
            .copy()
        )
        v = (
            self.sol_state[:, (self.objects_count + 1) * 3 :]
            .reshape(-1, (self.objects_count - 1), 3)
            .copy()
        )

        x = x - self.sol_state[:, :3].reshape(-1, 1, 3)
        v = v - self.sol_state[
            :, (self.objects_count) * 3 : (self.objects_count + 1) * 3
        ].reshape(-1, 1, 3)

        eccentricity = (
            np.cross(v, np.cross(x, v))
            / (self.G * (self.m[0] + self.m[1:]))[:, np.newaxis]
            - x / np.linalg.norm(x, axis=2)[:, :, np.newaxis]
        )
        self.eccentricity = np.linalg.norm(eccentricity, axis=2)

        stop = timeit.default_timer()
        print(f"Run time: {(stop - start):.3f} s")
        print("")

    def compute_inclination(self):
        """
        Compute the inclination using the sol_state array,
        assuming that the first object is the central object

        C library is not implemented since this is mostly NumPy
        vector operations, which can be done almost instantly.
        """
        print("Computing inclination (Assuming the first body is the central star)...")
        self.inclination = np.zeros(self.data_size)

        start = timeit.default_timer()
        x = (
            self.sol_state[:, 3 : (self.objects_count * 3)]
            .reshape(-1, (self.objects_count - 1), 3)
            .copy()
        )
        v = (
            self.sol_state[:, (self.objects_count + 1) * 3 :]
            .reshape(-1, (self.objects_count - 1), 3)
            .copy()
        )

        x = x - self.sol_state[:, :3].reshape(-1, 1, 3)
        v = v - self.sol_state[
            :, (self.objects_count) * 3 : (self.objects_count + 1) * 3
        ].reshape(-1, 1, 3)

        unit_angular_momentum_vector = (
            np.cross(x, v) / np.linalg.norm(np.cross(x, v), axis=2)[:, :, np.newaxis]
        )
        unit_z = np.array([0, 0, 1])

        self.inclination = np.arccos(
            np.sum(unit_angular_momentum_vector * unit_z, axis=2)
        )

        stop = timeit.default_timer()
        print(f"Run time: {(stop - start):.3f} s")
        print("")

    @staticmethod
    def combine_metadata_with_flushed_data(
        results_file_path: str,
        flushed_file_path: str,
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
        no_print: bool = False,
    ):
        common.save_results(
            file_path=results_file_path,
            system_name=system_name,
            integrator_name=integrator_name,
            objects_count=objects_count,
            G=G,
            tf=tf,
            dt=dt,
            tolerance=tolerance,
            data_size=data_size,
            store_every_n=store_every_n,
            run_time=run_time,
            masses=masses,
            only_metadata=True,
            no_print=no_print,
        )

        with open(results_file_path, "a", newline="") as results_file:
            writer = csv.writer(results_file)
            with open(flushed_file_path, "r", newline="") as flushed_file:
                reader = csv.reader(flushed_file)
                for row in reader:
                    writer.writerow(row)

    @property
    def integration_mode(self) -> str:
        return self._integration_mode

    @integration_mode.setter
    def integration_mode(self, integration_mode: str):
        match integration_mode:
            case "numpy":
                self._integration_mode = "numpy"
            case "c_lib":
                if self.c_lib is not None:
                    self._integration_mode = "c_lib"
                else:
                    warnings.warn("C library not available. Defaulting to numpy.")
                    self._integration_mode = "numpy"
            case _:
                raise ValueError(
                    'Invalid integration mode. Must be "numpy" or "c_lib".'
                )

    @property
    def data_size(self) -> int:
        return len(self.sol_time)
