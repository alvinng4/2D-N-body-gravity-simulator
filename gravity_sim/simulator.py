import ctypes
import threading
import timeit
import warnings

import numpy as np

from progress_bar import Progress_bar
from progress_bar import progress_bar_c_lib_fixed_step_size
from integrator_simple import SimpleIntegrator
from integrator_rk_embedded import RKEmbedded
from integrator_ias15 import IAS15


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
    ]

    def __init__(self, c_lib=None, is_exit_ctypes_bool=None):
        self.c_lib = c_lib
        self.is_c_lib = c_lib is not None

        if is_exit_ctypes_bool is None:
            self.is_exit_ctypes_bool = ctypes.c_bool(False)
        else:
            self.is_exit_ctypes_bool = is_exit_ctypes_bool

        self.run_time = None

    def launch_simulation(
        self,
        gravitational_system,
        integrator: str,
        tf: float,
        dt: float = None,
        tolerance: float = None,
        store_every_n: int = 1,
        acceleration: str = "pairwise",
        no_progress_bar: bool = False,
        no_print: bool = False,
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
        acceleration : str, optional
            Acceleration method -- "pairwise", "massless"
                pairwise: Pairwise acceleration
                massless: System with massless particles
        no_progress_bar : bool, optional
            Disable progress bar
        no_print : bool, optional
            Disable print statements
        """
        self.x = gravitational_system.x.copy()
        self.v = gravitational_system.v.copy()
        self.m = gravitational_system.m.copy()
        self.objects_count = gravitational_system.objects_count
        self.G = gravitational_system.G

        self.integrator = integrator
        self.tf = tf
        self.store_every_n = store_every_n

        if (dt is None) and (tolerance is None):
            raise ValueError("Either dt or tolerance must be provided.")

        self.dt = dt
        self.tolerance = tolerance

        if not self.is_c_lib and acceleration != "pairwise":
            warnings.warn(
                "Only pairwise acceleration is available for numpy integrators. "
                + 'Setting acceleration method to "pairwise".'
            )
            acceleration_func = self.c_lib.acceleration_pairwise
        elif self.is_c_lib:
            if acceleration == "pairwise":
                acceleration_func = self.c_lib.acceleration_pairwise
            elif acceleration == "massless":
                acceleration_func = self.c_lib.acceleration_with_massless
            elif acceleration == "barnes_hut":
                raise NotImplementedError

        if not no_print:
            print("Simulating the system...")
            start = timeit.default_timer()

        if self.is_c_lib:
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
                    ) = integrator.simulation_c_lib(
                        self.integrator,
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.dt,
                        self.tf,
                        acceleration_func,
                        no_progress_bar,
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
                    ) = integrator.simulation_c_lib(
                        integrator.order,
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.tf,
                        self.tolerance,
                        self.tolerance,
                        acceleration_func,
                        no_progress_bar,
                    )

                case "ias15":
                    integrator = IAS15(
                        self.store_every_n, self.c_lib, self.is_exit_ctypes_bool
                    )
                    (
                        self.sol_state,
                        self.sol_time,
                        self.sol_dt,
                    ) = integrator.simulation_c_lib(
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.tf,
                        self.tolerance,
                        acceleration_func,
                        no_progress_bar,
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
                        no_progress_bar,
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
                    ) = integrator.simulation_numpy(
                        integrator.order,
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.tf,
                        self.tolerance,
                        self.tolerance,
                        no_progress_bar,
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
                    ) = integrator.simulation_numpy(
                        self.objects_count,
                        self.x,
                        self.v,
                        self.m,
                        self.G,
                        self.tf,
                        self.tolerance,
                        no_progress_bar,
                    )

        if not no_print:
            stop = timeit.default_timer()
            self.run_time = stop - start
            print(f"Run time: {self.run_time:.3f} s")
            print("")

    def compute_energy(self):
        """
        Compute the total energy using the sol_state array
        """
        print("Computing energy...")
        npts = len(self.sol_state)
        self.energy = np.zeros(npts)

        start = timeit.default_timer()
        if self.is_c_lib == True:
            count = ctypes.c_int(0)
            progress_bar_thread = threading.Thread(
                target=progress_bar_c_lib_fixed_step_size,
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

        elif self.is_c_lib == False:
            progress_bar = Progress_bar()
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
        # Check if self.m has data
        if len(self.m) == 0:
            print("Error: Mass data is missing. Returning to menu.")
            print()
            return 1

        print("Computing linear momentum...")
        npts = len(self.sol_state)
        self.linear_momentum = np.zeros(npts)

        start = timeit.default_timer()
        if self.is_c_lib == True:
            count = ctypes.c_int(0)
            progress_bar_thread = threading.Thread(
                target=progress_bar_c_lib_fixed_step_size,
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

        elif self.is_c_lib == False:
            progress_bar = Progress_bar()
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

        return 0

    def compute_angular_momentum(self):
        """
        Compute the total angular momentum using the sol_state array
        """
        # Check if self.m has data
        if len(self.m) == 0:
            print("Error: Mass data is missing. Returning to menu.")
            print()
            return 1

        print("Computing angular momentum...")
        npts = len(self.sol_state)
        self.angular_momentum = np.zeros(npts)

        start = timeit.default_timer()
        if self.is_c_lib == True:
            count = ctypes.c_int(0)
            progress_bar_thread = threading.Thread(
                target=progress_bar_c_lib_fixed_step_size,
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

        elif self.is_c_lib == False:
            progress_bar = Progress_bar()
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

        return 0
