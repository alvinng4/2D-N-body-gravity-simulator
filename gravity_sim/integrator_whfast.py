"""
WHFast integrator

References: 
1. Moving Planets Around: An Introduction to N-Body Simulations
   Applied to Exoplanetary Systems, Chapter 9 and Appendix B
2. Rein & Tamayo 2015, https://arxiv.org/abs/1505.03036
"""

import ctypes
import math
import threading
import time
import sys
from queue import Queue

import numpy as np

from common import acceleration
from progress_bar import progress_bar_c_lib_fixed_step_size
from progress_bar import Progress_bar_with_data_size


class WHFast:
    def __init__(self, store_every_n=1, c_lib=None, is_exit_ctypes_bool=None):
        self.store_every_n = store_every_n
        self.c_lib = c_lib

        if is_exit_ctypes_bool is None:
            self.is_exit_ctypes_bool = ctypes.c_bool(False)
        else:
            self.is_exit_ctypes_bool = is_exit_ctypes_bool

    def simulation_c_lib(
        self,
        objects_count,
        x,
        v,
        m,
        G,
        dt,
        tf,
        acceleration_func,
        flush: bool = False,
        flush_path: str = "",
        no_progress_bar: bool = False,
    ):
        class Solutions(ctypes.Structure):
            _fields_ = [
                ("sol_state", ctypes.POINTER(ctypes.c_double)),
                ("sol_time", ctypes.POINTER(ctypes.c_double)),
                ("sol_dt", ctypes.POINTER(ctypes.c_double)),
            ]

        self.c_lib.whfast.restype = ctypes.c_int

        npts = int(np.floor((tf / dt))) + 1  # + 1 for t0
        if self.store_every_n != 1:
            store_npts = math.floor((npts - 1) / self.store_every_n) + 1  # + 1 for t0
        else:
            store_npts = npts

        store_count = ctypes.c_int(0)

        if not no_progress_bar:
            progress_bar_thread = threading.Thread(
                target=progress_bar_c_lib_fixed_step_size,
                args=(store_npts, store_count, self.is_exit_ctypes_bool),
            )
            progress_bar_thread.start()

        # parameters are double no matter what "real" is defined
        def whfast_wrapper(c_lib_integrator, return_queue, *args):
            return_code = c_lib_integrator(*args)
            return_queue.put(return_code)

        queue = Queue()
        solution = Solutions()
        flush_path = flush_path.encode("utf-8")

        whfast_thread = threading.Thread(
            target=whfast_wrapper,
            args=(
                self.c_lib.euler,
                queue,
                ctypes.c_int(objects_count),
                x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.c_double(G),
                ctypes.c_double(dt),
                acceleration_func,
                ctypes.c_int64(npts),
                ctypes.c_int(store_npts),
                ctypes.c_int(self.store_every_n),
                ctypes.byref(store_count),
                ctypes.c_bool(flush),
                flush_path,
                ctypes.byref(solution),
                ctypes.byref(self.is_exit_ctypes_bool),
            ),
        )

        whfast_thread.start()

        # Keeps the main thread running to catch keyboard interrupt
        # This is added since the main thread is not catching
        # exceptions on Windows
        while whfast_thread.is_alive():
            time.sleep(0.1)

        whfast_thread.join()
        return_code = queue.get()

        # Check return code
        match return_code:
            # Simulation successful
            case 0:
                pass
            # C library failed to allocate memory
            case 1:
                self.is_exit_ctypes_bool.value = True
                sys.exit(1)
            # User sent KeyboardInterrupt
            case 2:
                pass  # Should be caught in run_prog and exit

        # Close the progress_bar_thread
        if not no_progress_bar:
            if store_count.value < (store_npts - 1):
                store_count.value = store_npts - 1
            progress_bar_thread.join()

        if not flush:
            # Convert C arrays to numpy arrays
            return_sol_state = np.ctypeslib.as_array(
                solution.sol_state,
                shape=(store_npts, objects_count * 6),
            ).copy()

            return_sol_time = np.ctypeslib.as_array(
                solution.sol_time, shape=(store_npts,)
            ).copy()
            return_sol_dt = np.ctypeslib.as_array(
                solution.sol_dt, shape=(store_npts,)
            ).copy()

            # Free memory
            self.c_lib.free_memory_real(solution.sol_state)
            self.c_lib.free_memory_real(solution.sol_time)
            self.c_lib.free_memory_real(solution.sol_dt)

            return (
                return_sol_state,
                return_sol_time,
                return_sol_dt,
                store_count.value,
            )

        else:
            return (
                None,
                None,
                None,
                store_count.value,
            )

    def simulation_numpy(
        self,
        objects_count,
        x,
        v,
        m,
        G,
        dt,
        tf,
        flush: bool = False,
        flush_path: str = "",
        no_progress_bar: bool = False,
    ):
        if flush:
            raise NotImplementedError("Flush is not implemented for numpy")

        npts = int(np.floor((tf / dt))) + 1  # + 1 for t0
        if self.store_every_n != 1:
            store_npts = math.floor((npts - 1) / self.store_every_n) + 1  # + 1 for t0
        else:
            store_npts = npts
        self.sol_state = np.zeros((store_npts, objects_count * 3 * 2))
        self.sol_state[0] = np.concatenate(
            (
                np.reshape(x, objects_count * 3),
                np.reshape(v, objects_count * 3),
            )
        )
        self.sol_time = np.linspace(0, dt * (npts - 1), store_npts)
        self.sol_dt = np.full(shape=(store_npts), fill_value=f"{dt}", dtype=float)

        x_err_comp_sum = np.zeros((objects_count, 3))
        v_err_comp_sum = np.zeros((objects_count, 3))

        eta = np.cumsum(m)
        jacobi = np.zeros((objects_count, 6))
        a = np.zeros((objects_count, 3))
        WHFast.cartesian_to_jacobi(objects_count, jacobi, x, v, m, eta)
        WHFast.whfast_acceleration(objects_count, jacobi, x, a, m, eta, G)
        WHFast.whfast_kick(jacobi, a, 0.5 * dt)

        progress_bar = Progress_bar_with_data_size()
        store_count = 0
        with progress_bar:
            if not no_progress_bar:
                task = progress_bar.add_task("", total=store_npts, store_count=1)

            for count in range(npts - 1):
                WHFast.whfast_drift(objects_count, jacobi, m, eta, G, dt)
                WHFast.jacobi_to_cartesian(objects_count, jacobi, x, v, m, eta)
                WHFast.whfast_acceleration(objects_count, jacobi, x, a, m, eta, G)
                WHFast.whfast_kick(jacobi, a, dt)

                # Store solution
                if (count + 1) % self.store_every_n == 0:
                    temp_jacobi = jacobi.copy()
                    WHFast.whfast_kick(temp_jacobi, a, -0.5 * dt)
                    WHFast.jacobi_to_cartesian(objects_count, temp_jacobi, x, v, m, eta)
                    self.sol_state[store_count + 1] = np.concatenate(
                        (
                            np.reshape(x, objects_count * 3),
                            np.reshape(v, objects_count * 3),
                        )
                    )
                    store_count += 1

                if (count + 2) == npts:
                    temp_jacobi = jacobi.copy()
                    WHFast.whfast_kick(temp_jacobi, a, -0.5 * dt)
                    self.sol_state[-1] = np.concatenate(
                        (
                            np.reshape(x, objects_count * 3),
                            np.reshape(v, objects_count * 3),
                        )
                    )

                if not no_progress_bar:
                    progress_bar.update(
                        task, completed=store_count, store_count=store_count + 1
                    )

            if not no_progress_bar:
                progress_bar.update(task, completed=store_npts, store_count=store_npts)

        return self.sol_state, self.sol_time, self.sol_dt

    @staticmethod
    def whfast_kick(jacobi, a, dt) -> None:
        jacobi[:, 3:] += a * dt

    @staticmethod
    def whfast_drift(objects_count, jacobi, m, eta, G, dt) -> None:
        for i in range(1, objects_count):
            gm = G * m[0] * eta[i] / eta[i - 1]
            WHFast.propagate_kepler(jacobi[i], gm, dt)

    @staticmethod
    def whfast_acceleration(objects_count, jacobi, x, a, m, eta, G) -> None:
        """
        Compute the acceleration of the system

        Parameters
        ----------
        objects_count : int
            Number of objects
        jacobi : np.ndarray
            Array of state vector in Jacobi coordinates
        x : np.ndarray
            Array of Cartesian positions
        a : np.ndarray
            Array to store the acceleration
        m : np.ndarray
            Array of masses of the objects
        eta : np.ndarray
            Array of mass where eta[i] is the total mass from
            the 0-th to the i-th object
        G : float
            Gravitational constant
        """
        # fmt: off
        # Empty acceleration array
        a.fill(0.0)

        aux = np.zeros(3)
        for i in range(1, objects_count):
            x_0i = x[i] - x[0]
            a[i] = (
                m[0] * eta[i] / eta[i - 1]
                * (
                    jacobi[i, :3] / np.linalg.norm(jacobi[i, :3]) ** 3
                    - x_0i / np.linalg.norm(x_0i) ** 3
                )
            )

            aux.fill(0.0)
            for j in range(1, i):
                x_ji = x[i] - x[j]
                aux += m[j] * x_ji / np.linalg.norm(x_ji) ** 3
            a[i] -= eta[i] / eta[i - 1] * aux

            aux.fill(0.0)
            for j in range(i + 1, objects_count):
                x_ij = x[j] - x[i]
                aux += m[j] * x_ij / np.linalg.norm(x_ij) ** 3
            a[i] += aux

            aux.fill(0.0)
            for j in range(0, i):
                for k in range(i + 1, objects_count):
                    x_jk = x[k] - x[j]
                    aux += m[j] * m[k] * x_jk / np.linalg.norm(x_jk) ** 3
            a[i] -= aux / eta[i - 1]

        a *= G
        # fmt: on

    @staticmethod
    def cartesian_to_jacobi(
        objects_count: int,
        jacobi: np.ndarray,
        x: np.ndarray,
        v: np.ndarray,
        m: np.ndarray,
        eta: np.ndarray,
    ) -> None:
        """
        Transform Cartesian coordinates to Jacobi coordinates

        Parameters
        ----------
        objects_count : int
            Number of objects
        jacobi : np.ndarray
            The array to store the Jacobi coordinates
        x : np.ndarray
            The Cartesian position array
        v : np.ndarray
            The Cartesian velocity array
        m : np.ndarray
            Array of masses of the objects
        eta : np.ndarray
            Array of mass where eta[i] is the total mass from
            the 0-th to the i-th object
        """
        x_cm = m[0] * x[0]
        v_cm = m[0] * v[0]

        for i in range(1, objects_count):
            jacobi[i, :3] = x[i] - x_cm / eta[i - 1]
            jacobi[i, 3:] = v[i] - v_cm / eta[i - 1]

            x_cm = x_cm * (1 + m[i] / eta[i - 1]) + m[i] * jacobi[i, :3]
            v_cm = v_cm * (1 + m[i] / eta[i - 1]) + m[i] * jacobi[i, 3:]

        jacobi[0, :3] = x_cm / eta[-1]
        jacobi[0, 3:] = v_cm / eta[-1]

    @staticmethod
    def jacobi_to_cartesian(
        objects_count: int,
        jacobi: np.ndarray,
        x: np.ndarray,
        v: np.ndarray,
        m: np.ndarray,
        eta: np.ndarray,
    ) -> None:
        """
        Transform Jacobi coordinates to Cartesian coordinates

        Parameters
        ----------
        objects_count : int
            Number of objects
        jacobi : np.ndarray
            The array of Jacobi coordinates
        x : np.ndarray
            The array to store the Cartesian position
        v : np.ndarray
            The array to store the Cartesian velocity
        m : np.ndarray
            Array of masses of the objects
        eta : np.ndarray
            Array of mass where eta[i] is the total mass from
            the 0-th to the i-th object
        """
        x_cm = eta[-1] * jacobi[0, 0:3]
        v_cm = eta[-1] * jacobi[0, 3:6]

        for i in range(objects_count - 1, 0, -1):
            x_cm = (x_cm - m[i] * jacobi[i, :3]) / eta[i]
            v_cm = (v_cm - m[i] * jacobi[i, 3:]) / eta[i]

            x[i] = jacobi[i, :3] + x_cm
            v[i] = jacobi[i, 3:] + v_cm

            x_cm = eta[i - 1] * x_cm
            v_cm = eta[i - 1] * v_cm

        x[0] = x_cm / m[0]
        v[0] = v_cm / m[0]

    @staticmethod
    def stumpff_functions(z: float) -> tuple[float, float, float, float]:
        """
        Compute the Stumpff functions c0, c1, c2, and c3 for a given argument z.

        Parameters
        ----------
        z : float
            Argument of the Stumpff functions.

        Returns
        -------
        tuple[float, float, float, float]
            Tuple containing the Stumpff functions c0, c1, c2, and c3.
        """
        # # Reduce the argument
        # n = 0
        # while abs(z) > 0.1:
        #     n += 1
        #     z /= 4.0

        # # Compute c3, c2, c1, c0
        # # fmt: off
        # c3 = (1.0 - z / 20.0 * (1.0 - z / 42.0 * (1.0 - z / 72.0 * (1.0 - z / 110.0 \
        #          * (1.0 - z / 156.0 * (1.0 - z / 210.0)))))
        #      ) / 6.0
        # c2 = (1.0 - z / 12.0 * (1.0 - z / 30.0 * (1.0 - z / 56.0 * (1.0 - z / 90.0 \
        #          * (1.0 - z / 132.0 * (1.0 - z / 182.0)))))
        #      ) / 2.0
        # c1 = 1.0 - z * c3
        # c0 = 1.0 - z * c2

        # # Half-angle formulae to recover the actual argument
        # while n > 0:
        #     n -= 1
        #     c3 = (c2 + c0 * c3) / 4.0
        #     c2 = c1 * c1 / 2.0
        #     c1 = c0 * c1
        #     c0 = 2.0 * c0 * c0 - 1.0

        # return c0, c1, c2, c3
        # # fmt: on

        ########################################################
        # Analytical expressions for the Stumpff functions
        if z > 0:
            c3 = (math.sqrt(z) - math.sin(math.sqrt(z))) / (z ** 1.5)
            c2 = (1.0 - math.cos(math.sqrt(z))) / z
            c1 = math.sin(math.sqrt(z)) / math.sqrt(z)
            c0 = math.cos(math.sqrt(z))

        elif z < 0:
            c3 = (math.sqrt(-z) - math.sinh(math.sqrt(-z))) / (z * math.sqrt(-z))
            c2 = (1.0 - math.cosh(math.sqrt(-z))) / z
            c1 = math.sinh(math.sqrt(-z)) / math.sqrt(-z)
            c0 = math.cosh(math.sqrt(-z))

        else:
            c3 = 1.0 / 6.0
            c2 = 0.5
            c1 = 1.0
            c0 = 1.0

        return c0, c1, c2, c3
        ########################################################

    @staticmethod
    def propagate_kepler(
        jacobi_i: np.ndarray, gm: float, dt: float
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Propagate the position and velocity vectors using Kepler's equation.

        Parameters
        ----------
        jacobi : np.ndarray
            State vector in Jacobi coordinates
        gm : float
            Gravitational parameter
        dt : float
            Time step

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            New position and velocity vectors
        """
        # Internal tolerance for solving Kepler's equation
        tol = 1e-12
        
        x = jacobi_i[:3].copy()
        v = jacobi_i[3:].copy()
        x_norm = np.linalg.norm(x)
        v_norm = np.linalg.norm(v)

        alpha = 2.0 * gm / x_norm - v_norm**2

        # Radial velocity
        radial_v = np.dot(x, v) / x_norm

        # Solve Kepler's equation with Newton-Raphson method
        s = dt / x_norm  # initial guess
        for i in range(50):
            # Compute Stumpff functions
            c0, c1, c2, c3 = WHFast.stumpff_functions(alpha * s * s)

            # Evaluate Kepler's equation
            F = (
                x_norm * s * c1
                + x_norm * radial_v * (s * s) * c2
                + gm * (s * s * s) * c3
                - dt
            )

            # Check convergence
            if abs(F) < tol and i > 0:
                break

            # Compute derivative
            dF = x_norm * c0 + x_norm * radial_v * s * c1 + gm * (s * s) * c2

            # Find step size
            ds = -F / dF

            # Advance step
            s += ds

            # Check convergence
            if abs(ds) < tol:
                break

        # The radial distance is equal to the derivative of F
        r = dF

        # Evaluate f and g functions
        f = 1.0 - gm * (s * s) * c2 / x_norm
        g = dt - gm * (s * s * s) * c3

        df = -gm / (r * x_norm) * s * c1
        dg = 1.0 - gm / r * (s * s) * c2

        # Compute position and velocity vectors
        jacobi_i[:3] = f * x + g * v
        jacobi_i[3:] = df * x + dg * v
