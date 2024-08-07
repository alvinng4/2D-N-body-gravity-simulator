"""
Embedded Runge-Kutta Integrators

Some of the codes are adapted from the following book with modifications: 
Moving Planets Around: An Introduction to N-Body 
Simulations Applied to Exoplanetary Systems, Chapter 6
"""

import ctypes
import math
import sys
import threading
import time
from queue import Queue

import numpy as np

import common
from common import acceleration


class RKEmbedded:
    def __init__(
        self,
        integrator: str,
        store_every_n: int = 1,
        c_lib: ctypes.CDLL = None,
        is_exit_ctypes_bool: ctypes.c_bool = None,
    ) -> None:
        self.store_every_n = store_every_n
        self.c_lib = c_lib

        if is_exit_ctypes_bool is None:
            self.is_exit_ctypes_bool = ctypes.c_bool(False)
        else:
            self.is_exit_ctypes_bool = is_exit_ctypes_bool

        match integrator:
            case "rkf45":
                self.order = 45
            case "dopri":
                self.order = 54
            case "dverk":
                self.order = 65
            case "rkf78":
                self.order = 78
            case _:
                raise ValueError("Invalid integrator!")

    def simulation_c_lib(
        self,
        order: int,
        objects_count: int,
        x: np.ndarray,
        v: np.ndarray,
        m: np.ndarray,
        G: float,
        tf: float,
        abs_tolerance: float,
        rel_tolerance: float,
        acceleration_method: str,
        storing_method: int = 0,
        flush_path: str = "",
        no_progress_bar: bool = False,
    ):
        if acceleration_method not in ["pairwise", "massless", "barnes-hut"]:
            raise ValueError("Invalid acceleration method")

        class Solutions(ctypes.Structure):
            _fields_ = [
                ("sol_state", ctypes.POINTER(ctypes.c_double)),
                ("sol_time", ctypes.POINTER(ctypes.c_double)),
                ("sol_dt", ctypes.POINTER(ctypes.c_double)),
            ]

        self.c_lib.rk_embedded.restype = ctypes.c_int

        t = ctypes.c_double(0.0)
        store_count = ctypes.c_int(1)  # 1 for t0

        if not no_progress_bar:
            progress_bar_thread = threading.Thread(
                target=common.progress_bar_c_lib_adaptive_step_size,
                args=(tf, t, store_count, self.is_exit_ctypes_bool),
            )
            progress_bar_thread.start()

        # parameters are double no matter what "real" is defined
        def rk_embedded_wrapper(c_lib_rk_embedded, return_queue, *args):
            return_code = c_lib_rk_embedded(*args)
            return_queue.put(return_code)

        queue = Queue()
        solution = Solutions()

        rk_embedded_thread = threading.Thread(
            target=rk_embedded_wrapper,
            args=(
                self.c_lib.rk_embedded,
                queue,
                ctypes.c_int(order),
                ctypes.c_int(objects_count),
                x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.c_double(G),
                ctypes.byref(t),
                ctypes.c_double(tf),
                ctypes.c_double(abs_tolerance),
                ctypes.c_double(rel_tolerance),
                acceleration_method.encode("utf-8"),
                ctypes.c_int(self.store_every_n),
                ctypes.byref(store_count),
                ctypes.c_int(storing_method),
                flush_path.encode("utf-8"),
                ctypes.byref(solution),
                ctypes.byref(self.is_exit_ctypes_bool),
            ),
        )
        rk_embedded_thread.start()

        # Keeps the main thread running to catch keyboard interrupt
        # This is added since the main thread is not catching
        # exceptions on Windows
        while rk_embedded_thread.is_alive():
            time.sleep(0.05)

        rk_embedded_thread.join()
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
            t.value = tf
            progress_bar_thread.join()

        # Default storing
        if storing_method == 0:
            # Convert C arrays to numpy arrays
            return_sol_state = np.ctypeslib.as_array(
                solution.sol_state,
                shape=(store_count.value, objects_count * 6),
            ).copy()

            return_sol_time = np.ctypeslib.as_array(
                solution.sol_time, shape=(store_count.value,)
            ).copy()
            return_sol_dt = np.ctypeslib.as_array(
                solution.sol_dt, shape=(store_count.value,)
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
        order,
        objects_count,
        x,
        v,
        m,
        G,
        tf,
        abs_tolerance,
        rel_tolerance,
        storing_method: int = 0,
        flush_path: str = "",
        no_progress_bar: bool = False,
    ):
        if storing_method == 1:
            raise NotImplementedError("Flush is not implemented for numpy")
        if storing_method == 2:
            raise NotImplementedError("no_store is not implemented for numpy")

        (
            self.power,
            self.power_test,
            self.coeff,
            self.weights,
            self.weights_test,
        ) = self.rk_embedded_butcher_tableaus_numpy(order)

        a = acceleration(objects_count, x, m, G)
        initial_dt = self.rk_embedded_initial_dt_numpy(
            objects_count,
            self.power,
            x,
            v,
            a,
            m,
            G,
            abs_tolerance,
            rel_tolerance,
        )

        progress_bar = common.Progress_bar_with_data_size()
        with progress_bar:
            return self.simulation(
                progress_bar,
                objects_count,
                x,
                v,
                m,
                G,
                tf,
                initial_dt,
                self.power,
                self.power_test,
                self.coeff,
                self.weights,
                self.weights_test,
                abs_tolerance,
                rel_tolerance,
                self.store_every_n,
                no_progress_bar,
            )

    @staticmethod
    def simulation(
        progress_bar,
        objects_count: int,
        x,
        v,
        m,
        G,
        tf: float,
        initial_dt,
        power,
        power_test,
        coeff,
        weights,
        weights_test,
        abs_tolerance: float,
        rel_tolerance: float,
        store_every_n: int,
        no_progress_bar: bool,
    ):
        """
        Perform simulation using rk_embedded methods

        :return: sol_state, sol_time, sol_dt
        :rtype: numpy.array
        """
        # Initialization
        t = 0.0
        dt = initial_dt
        stages = len(weights)
        min_power = min([power, power_test])
        error_estimation_delta_weights = weights - weights_test

        # Safety factors for step-size control:
        safety_fac_max = 6.0
        safety_fac_min = 0.33
        safety_fac = 0.38 ** (1.0 / (1.0 + min_power))

        # Initialize vk and xk
        vk = np.zeros((stages, objects_count, 3))
        xk = np.zeros((stages, objects_count, 3))

        # Allocate for dense output:
        buffer_size = 50000
        sol_state = np.zeros((buffer_size, objects_count * 2 * 3))
        sol_time = np.zeros(buffer_size)
        sol_dt = np.zeros(buffer_size)

        # Initial values
        sol_state[0] = np.concatenate(
            (np.reshape(x, objects_count * 3), np.reshape(v, objects_count * 3))
        )
        sol_time[0] = t
        sol_dt[0] = 0.0

        # Arrays for compensated summation
        x_err_comp_sum = np.zeros((objects_count, 3))
        v_err_comp_sum = np.zeros((objects_count, 3))
        temp_x_err_comp_sum = np.zeros((objects_count, 3))
        temp_v_err_comp_sum = np.zeros((objects_count, 3))

        # Launch integration:
        count = 1  # 1 for t0
        store_count = 1  # 1 for t0
        if not no_progress_bar:
            task = progress_bar.add_task("", total=tf, store_count=store_count)

        while t < tf:
            # Calculate xk and vk
            vk[0] = acceleration(objects_count, x, m, G)
            xk[0] = np.copy(v)
            for stage in range(1, stages):
                temp_v = np.zeros((objects_count, 3))
                temp_x = np.zeros((objects_count, 3))
                for j in range(stage):
                    temp_v += coeff[stage - 1][j] * vk[j]
                    temp_x += coeff[stage - 1][j] * xk[j]
                vk[stage] = acceleration(objects_count, x + dt * temp_x, m, G)
                xk[stage] = v + dt * temp_v

            # Calculate x_1, v_1 and also delta x, delta v for error estimation
            temp_v = np.zeros((objects_count, 3))
            temp_x = np.zeros((objects_count, 3))
            error_estimation_delta_x = np.zeros((objects_count, 3))
            error_estimation_delta_v = np.zeros((objects_count, 3))
            for stage in range(stages):
                temp_v += weights[stage] * vk[stage]
                temp_x += weights[stage] * xk[stage]
                error_estimation_delta_v += (
                    error_estimation_delta_weights[stage] * vk[stage]
                )
                error_estimation_delta_x += (
                    error_estimation_delta_weights[stage] * xk[stage]
                )

            error_estimation_delta_v *= dt
            error_estimation_delta_x *= dt

            # Calculate x_1 and v_1
            temp_v_err_comp_sum = v_err_comp_sum
            temp_v_err_comp_sum += dt * temp_v
            v_1 = v + temp_v_err_comp_sum
            temp_v_err_comp_sum += v - v_1

            temp_x_err_comp_sum = x_err_comp_sum
            temp_x_err_comp_sum += dt * temp_x
            x_1 = x + temp_x_err_comp_sum
            temp_x_err_comp_sum += x - x_1

            # Error calculation
            tolerance_scale_v = (
                abs_tolerance + np.maximum(np.abs(v), np.abs(v_1)) * rel_tolerance
            )
            tolerance_scale_x = (
                abs_tolerance + np.maximum(np.abs(x), np.abs(x_1)) * rel_tolerance
            )

            # Sum up all the elements of x/tol and v/tol, square and divide by the total number of elements
            sum = np.sum(
                np.square(error_estimation_delta_x / tolerance_scale_x)
            ) + np.sum(np.square(error_estimation_delta_v / tolerance_scale_v))
            error = math.sqrt(sum / (objects_count * 3 * 2))

            if error <= 1 or dt <= tf * 1e-12:
                t += dt
                x = x_1
                v = v_1

                x_err_comp_sum = temp_x_err_comp_sum
                v_err_comp_sum = temp_v_err_comp_sum

                # Store step:
                if (count % store_every_n) == 0:
                    sol_state[store_count] = np.concatenate(
                        (
                            np.reshape(x, objects_count * 3),
                            np.reshape(v, objects_count * 3),
                        )
                    )
                    sol_time[store_count] = t
                    sol_dt[store_count] = dt
                    store_count += 1

                    if not no_progress_bar:
                        progress_bar.update(task, completed=t, store_count=store_count)

                    # Check buffer size and extend if needed:
                    if (store_count == buffer_size) and (t < tf):
                        new_buffer_size = buffer_size * 2
                        sol_state = np.pad(
                            sol_state, ((0, new_buffer_size - buffer_size), (0, 0))
                        )
                        sol_time = np.pad(sol_time, (0, new_buffer_size - buffer_size))
                        sol_dt = np.pad(sol_dt, (0, new_buffer_size - buffer_size))
                        buffer_size = new_buffer_size

                count += 1

            if error < 1e-10:
                error = 1e-10  # Prevent error from being too small

            dt_new = dt * safety_fac / error ** (1.0 / (1.0 + min_power))
            # Prevent dt to be too small or too large relative to the last time step
            if dt_new > safety_fac_max * dt:
                dt = safety_fac_max * dt
            elif dt_new < safety_fac_min * dt:
                dt = safety_fac_min * dt
            else:
                dt = dt_new

            if dt_new / tf < 1e-12:
                dt = tf * 1e-12

            # Correct overshooting:
            if t + dt > tf:
                dt = tf - t

        if not no_progress_bar:
            progress_bar.update(task, completed=tf, store_count=store_count)
        return (
            sol_state[0:store_count],
            sol_time[0:store_count],
            sol_dt[0:store_count],
        )

    @staticmethod
    def rk_embedded_initial_dt_numpy(
        objects_count: int,
        power: int,
        x,
        v,
        a,
        m,
        G,
        abs_tolerance: float,
        rel_tolerance: float,
    ) -> float:
        """
        Calculate the initial time step for embedded rk method

        Modified: Return dt * 1e-2 since this function gives initial dt thats too large
        """
        tolerance_scale_x = abs_tolerance + rel_tolerance * np.abs(x)
        tolerance_scale_v = abs_tolerance + rel_tolerance * np.abs(v)
        sum_0 = np.sum(np.square(x / tolerance_scale_x)) + np.sum(
            np.square(v / tolerance_scale_v)
        )
        sum_1 = np.sum(np.square(v / tolerance_scale_x)) + np.sum(
            np.square(a / tolerance_scale_v)
        )
        d_0 = (sum_0 / (objects_count * 3 * 2)) ** 0.5
        d_1 = (sum_1 / (objects_count * 3 * 2)) ** 0.5

        if d_0 < 1e-5 or d_1 < 1e-5:
            dt_0 = 1e-4
        else:
            dt_0 = d_0 / d_1

        x_1 = x + (dt_0 / 100) * v
        v_1 = v + (dt_0 / 100) * a
        a_1 = acceleration(objects_count, x_1, m, G)

        # Calculate d_2 to measure how much the derivatives have changed.
        sum_2 = np.sum(np.square((v_1 - v) / tolerance_scale_x)) + np.sum(
            np.square((a_1 - a) / tolerance_scale_v)
        )
        d_2 = (sum_2 / (objects_count * 3 * 2)) ** 0.5 / dt_0

        if max(d_1, d_2) <= 1e-15:
            dt_1 = max(1e-6, dt_0 * 1e-3)
        else:
            dt_1 = (0.01 / max(d_1, d_2)) ** (1.0 / (1.0 + power))
        dt = min(100 * dt_0, dt_1)

        return dt * 1e-2

    @staticmethod
    def rk_embedded_butcher_tableaus_numpy(order):
        """
        Butcher tableaus for embedded rk

        :raise ValueError: If order is not in [45, 54, 78, 65]
        :return: power, power_test, coeff, weights, weights_test
        :rtype: numpy.array
        """
        # Select integrator
        # 45) Runge-Kutta-Fehleberg 4(5)
        # 54) Dormand-Prince 5(4)
        # 78) Runge-Kutta-Fehlberg 7(8)
        # 65) Verner's method 6(5), DVERK

        match order:
            # RUNGE-KUTTA-FEHLBERG 4(5)
            case 45:
                # Order
                power = 4
                power_test = 5
                # nodes = np.array([1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 0.5])
                coeff = np.array(
                    [
                        [1.0 / 4.0, 0.0, 0.0, 0.0, 0.0],
                        [3.0 / 32.0, 9.0 / 32.0, 0.0, 0.0, 0.0],
                        [1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0, 0.0, 0.0],
                        [439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0, 0.0],
                        [
                            -8.0 / 27.0,
                            2.0,
                            -3544.0 / 2565.0,
                            1859.0 / 4104.0,
                            -11.0 / 40.0,
                        ],
                    ]
                )

                weights = np.array(
                    [25.0 / 216.0, 0.0, 1408.0 / 2565.0, 2197.0 / 4104.0, -0.2, 0.0]
                )
                weights_test = np.array(
                    [
                        16.0 / 135.0,
                        0.0,
                        6656.0 / 12825.0,
                        28561.0 / 56430.0,
                        -9.0 / 50.0,
                        2.0 / 55.0,
                    ]
                )

            # DORMAND-PRINCE 5(4)
            case 54:
                # order
                power = 5
                power_test = 4
                # nodes = np.array([1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0])
                coeff = np.array(
                    [
                        [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
                        [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
                        [
                            19372.0 / 6561.0,
                            -25360.0 / 2187.0,
                            64448.0 / 6561.0,
                            -212.0 / 729.0,
                            0.0,
                            0.0,
                        ],
                        [
                            9017.0 / 3168.0,
                            -355.0 / 33.0,
                            46732.0 / 5247.0,
                            49.0 / 176.0,
                            -5103.0 / 18656.0,
                            0.0,
                        ],
                        [
                            35.0 / 384.0,
                            0.0,
                            500.0 / 1113.0,
                            125.0 / 192.0,
                            -2187.0 / 6784.0,
                            11.0 / 84.0,
                        ],
                    ]
                )
                weights = np.array(
                    [
                        35.0 / 384.0,
                        0.0,
                        500.0 / 1113.0,
                        125.0 / 192.0,
                        -2187.0 / 6784.0,
                        11.0 / 84.0,
                        0.0,
                    ]
                )
                weights_test = np.array(
                    [
                        5179.0 / 57600.0,
                        0.0,
                        7571.0 / 16695.0,
                        393.0 / 640.0,
                        -92097.0 / 339200.0,
                        187.0 / 2100.0,
                        1.0 / 40.0,
                    ]
                )

            # RUNGE-KUTTA-FEHLBERG 7(8)
            case 78:
                # Order
                power = 7
                power_test = 8
                # nodes = np.array(
                #     [
                #         2.0 / 27.0,
                #         1.0 / 9.0,
                #         1.0 / 6.0,
                #         5.0 / 12.0,
                #         1.0 / 2.0,
                #         5.0 / 6.0,
                #         1.0 / 6.0,
                #         2.0 / 3.0,
                #         1.0 / 3.0,
                #         1.0,
                #         0.0,
                #         1.0,
                #     ]
                # )
                coeff = np.array(
                    [
                        [
                            2.0 / 27.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                        ],
                        [
                            1.0 / 36.0,
                            1.0 / 12.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                        ],
                        [
                            1.0 / 24.0,
                            0.0,
                            1.0 / 8.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                        ],
                        [
                            5.0 / 12.0,
                            0.0,
                            -25.0 / 16.0,
                            25.0 / 16.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                        ],
                        [
                            1.0 / 20.0,
                            0.0,
                            0.0,
                            1.0 / 4.0,
                            1.0 / 5.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                        ],
                        [
                            -25.0 / 108.0,
                            0.0,
                            0.0,
                            125.0 / 108.0,
                            -65.0 / 27.0,
                            125.0 / 54.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                        ],
                        [
                            31.0 / 300.0,
                            0.0,
                            0.0,
                            0.0,
                            61.0 / 225.0,
                            -2.0 / 9.0,
                            13.0 / 900.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                        ],
                        [
                            2.0,
                            0.0,
                            0.0,
                            -53.0 / 6.0,
                            704.0 / 45.0,
                            -107.0 / 9.0,
                            67.0 / 90.0,
                            3.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                        ],
                        [
                            -91.0 / 108.0,
                            0.0,
                            0.0,
                            23.0 / 108.0,
                            -976.0 / 135.0,
                            311.0 / 54.0,
                            -19.0 / 60.0,
                            17.0 / 6.0,
                            -1.0 / 12.0,
                            0.0,
                            0.0,
                            0.0,
                        ],
                        [
                            2383.0 / 4100.0,
                            0.0,
                            0.0,
                            -341.0 / 164.0,
                            4496.0 / 1025.0,
                            -301.0 / 82.0,
                            2133.0 / 4100.0,
                            45.0 / 82.0,
                            45.0 / 164.0,
                            18.0 / 41.0,
                            0.0,
                            0.0,
                        ],
                        [
                            3.0 / 205.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            -6.0 / 41.0,
                            -3.0 / 205.0,
                            -3.0 / 41.0,
                            3.0 / 41.0,
                            6.0 / 41.0,
                            0.0,
                            0.0,
                        ],
                        [
                            -1777.0 / 4100.0,
                            0.0,
                            0.0,
                            -341.0 / 164.0,
                            4496.0 / 1025.0,
                            -289.0 / 82.0,
                            2193.0 / 4100.0,
                            51.0 / 82.0,
                            33.0 / 164.0,
                            19.0 / 41.0,
                            0.0,
                            1.0,
                        ],
                    ]
                )

                weights = np.array(
                    [
                        41.0 / 840.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        34.0 / 105.0,
                        9.0 / 35.0,
                        9.0 / 35.0,
                        9.0 / 280.0,
                        9.0 / 280.0,
                        41.0 / 840.0,
                        0.0,
                        0.0,
                    ]
                )
                weights_test = np.array(
                    [
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        34.0 / 105.0,
                        9.0 / 35.0,
                        9.0 / 35.0,
                        9.0 / 280.0,
                        9.0 / 280.0,
                        0.0,
                        41.0 / 840.0,
                        41.0 / 840.0,
                    ]
                )

            # VERNER 6(5) DVERK
            case 65:
                # Order
                power = 6
                power_test = 7
                # nodes = np.array(
                #     [1.0 / 6.0, 4.0 / 15.0, 2.0 / 3.0, 5.0 / 6.0, 1.0, 1.0 / 15.0, 1.0]
                # )
                coeff = np.array(
                    [
                        [1.0 / 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        [4.0 / 75.0, 16.0 / 75.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        [5.0 / 6.0, -8.0 / 3.0, 5.0 / 2.0, 0.0, 0.0, 0.0, 0.0],
                        [
                            -165.0 / 64.0,
                            55.0 / 6.0,
                            -425.0 / 64.0,
                            85.0 / 96.0,
                            0.0,
                            0.0,
                            0.0,
                        ],
                        [
                            12.0 / 5.0,
                            -8.0,
                            4015.0 / 612.0,
                            -11.0 / 36.0,
                            88.0 / 255.0,
                            0.0,
                            0.0,
                        ],
                        [
                            -8263.0 / 15000.0,
                            124.0 / 75.0,
                            -643.0 / 680.0,
                            -81.0 / 250.0,
                            2484.0 / 10625.0,
                            0.0,
                            0.0,
                        ],
                        [
                            3501.0 / 1720.0,
                            -300.0 / 43.0,
                            297275.0 / 52632.0,
                            -319.0 / 2322.0,
                            24068.0 / 84065.0,
                            0.0,
                            3850.0 / 26703.0,
                        ],
                    ]
                )

                weights = np.array(
                    [
                        3.0 / 40.0,
                        0.0,
                        875.0 / 2244.0,
                        23.0 / 72.0,
                        264.0 / 1955.0,
                        0.0,
                        125.0 / 11592.0,
                        43.0 / 616.0,
                    ]
                )
                weights_test = np.array(
                    [
                        13.0 / 160.0,
                        0.0,
                        2375.0 / 5984.0,
                        5.0 / 16.0,
                        12.0 / 85.0,
                        3.0 / 44.0,
                        0.0,
                        0.0,
                    ]
                )

            case _:
                raise ValueError

        return power, power_test, coeff, weights, weights_test
