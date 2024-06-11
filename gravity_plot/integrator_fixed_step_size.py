"""
Fixed step size integrator:
Euler, Euler-Cromer, Runge-Kutta 4th order, Leapfrog
"""

import ctypes
import math
import threading

import numpy as np

from common import acceleration
from progress_bar import progress_bar_c_lib_fixed_integrator
from progress_bar import Progress_bar_with_data_size


class FIXED_STEP_SIZE_INTEGRATOR:
    def __init__(self, simulator):
        self.store_every_n = simulator.store_every_n
        if simulator.is_c_lib:
            self.c_lib = simulator.c_lib
            (
                self.sol_state,
                self.sol_time,
                self.sol_dt,
                self.m,
                self.G,
                self.objects_count,
            ) = self.simulation_c_lib(
                simulator.system.encode("utf-8"),
                simulator.integrator,
                simulator.dt,
                simulator.tf,
                simulator.x,  #######################
                simulator.v,  #                     #
                simulator.m,  #  For Custom system  #
                simulator.G,  #                     #
                simulator.objects_count,  #######################
            )
        else:
            self.simulation_numpy(
                simulator.integrator,
                simulator.objects_count,
                simulator.x,
                simulator.v,
                simulator.m,
                simulator.G,
                simulator.dt,
                simulator.tf,
            )

    def simulation_c_lib(
        self,
        system_name,
        integrator,
        dt,
        tf,
        custom_sys_x,
        custom_sys_v,
        custom_sys_m,
        custom_sys_G,
        custom_sys_objects_count,
    ):
        class Solutions(ctypes.Structure):
            _fields_ = [
                ("sol_state", ctypes.POINTER(ctypes.c_double)),
                ("sol_time", ctypes.POINTER(ctypes.c_double)),
                ("sol_dt", ctypes.POINTER(ctypes.c_double)),
                ("m", ctypes.POINTER(ctypes.c_double)),
                ("G", ctypes.c_double),
                ("objects_count", ctypes.c_int),
            ]

        self.c_lib.euler.restype = Solutions
        self.c_lib.euler_cromer.restype = Solutions
        self.c_lib.rk4.restype = Solutions
        self.c_lib.leapfrog.restype = Solutions

        npts = int(np.floor((tf / dt))) + 1  # + 1 for t0
        if self.store_every_n != 1:
            store_npts = math.floor((npts - 1) / self.store_every_n) + 1  # + 1 for t0
        else:
            store_npts = npts

        store_count = ctypes.c_int(0)

        progress_bar_thread = threading.Thread(
            target=progress_bar_c_lib_fixed_integrator, args=(store_npts, store_count)
        )
        progress_bar_thread.start()

        # parameters are double no matter what "real" is defined
        match integrator:
            case "euler":
                solutions = self.c_lib.euler(
                    ctypes.c_char_p(system_name),
                    ctypes.c_double(dt),
                    ctypes.c_int64(npts),
                    ctypes.c_int(store_npts),
                    ctypes.c_int(self.store_every_n),
                    ctypes.byref(store_count),
                    custom_sys_x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    custom_sys_v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    custom_sys_m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    ctypes.c_double(custom_sys_G),
                    ctypes.c_int(custom_sys_objects_count),
                )
            case "euler_cromer":
                solutions = self.c_lib.euler_cromer(
                    ctypes.c_char_p(system_name),
                    ctypes.c_double(dt),
                    ctypes.c_int64(npts),
                    ctypes.c_int(store_npts),
                    ctypes.c_int(self.store_every_n),
                    ctypes.byref(store_count),
                    custom_sys_x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    custom_sys_v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    custom_sys_m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    ctypes.c_double(custom_sys_G),
                    ctypes.c_int(custom_sys_objects_count),
                )
            case "rk4":
                solutions = self.c_lib.rk4(
                    ctypes.c_char_p(system_name),
                    ctypes.c_double(dt),
                    ctypes.c_int64(npts),
                    ctypes.c_int(store_npts),
                    ctypes.c_int(self.store_every_n),
                    ctypes.byref(store_count),
                    custom_sys_x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    custom_sys_v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    custom_sys_m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    ctypes.c_double(custom_sys_G),
                    ctypes.c_int(custom_sys_objects_count),
                )
            case "leapfrog":
                solutions = self.c_lib.leapfrog(
                    ctypes.c_char_p(system_name),
                    ctypes.c_double(dt),
                    ctypes.c_int64(npts),
                    ctypes.c_int(store_npts),
                    ctypes.c_int(self.store_every_n),
                    ctypes.byref(store_count),
                    custom_sys_x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    custom_sys_v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    custom_sys_m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    ctypes.c_double(custom_sys_G),
                    ctypes.c_int(custom_sys_objects_count),
                )

        # Close the thread forcefully if not closed
        if store_count.value < (store_npts - 1):
            store_count.value = store_npts - 1
        progress_bar_thread.join()

        # Convert C arrays to numpy arrays
        return_sol_state = (
            np.ctypeslib.as_array(
                solutions.sol_state,
                shape=(store_npts, solutions.objects_count * 6),
            )
            .copy()
            .reshape(store_npts, solutions.objects_count * 6)
        )

        return_sol_time = np.ctypeslib.as_array(
            solutions.sol_time, shape=(store_npts,)
        ).copy()
        return_sol_dt = np.ctypeslib.as_array(
            solutions.sol_dt, shape=(store_npts,)
        ).copy()
        return_m = np.ctypeslib.as_array(
            solutions.m, shape=(solutions.objects_count,)
        ).copy()

        # Free memory
        self.c_lib.free_memory_real(solutions.sol_state)
        self.c_lib.free_memory_real(solutions.sol_time)
        self.c_lib.free_memory_real(solutions.sol_dt)
        self.c_lib.free_memory_real(solutions.m)

        return (
            return_sol_state,
            return_sol_time,
            return_sol_dt,
            return_m,
            solutions.G,
            solutions.objects_count,
        )

    def simulation_numpy(self, integrator, objects_count, x, v, m, G, dt, tf):
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

        progress_bar = Progress_bar_with_data_size()
        store_count = 0
        with progress_bar:
            task = progress_bar.add_task("", total=store_npts, store_count=1)
            match integrator:
                case "euler":
                    for count in range(npts - 1):
                        a = acceleration(objects_count, x, m, G)
                        x, v, x_err_comp_sum, v_err_comp_sum = self.euler(
                            x, v, a, dt, x_err_comp_sum, v_err_comp_sum
                        )
                        if (count + 1) % self.store_every_n == 0:
                            self.sol_state[store_count + 1] = np.concatenate(
                                (
                                    np.reshape(x, objects_count * 3),
                                    np.reshape(v, objects_count * 3),
                                )
                            )
                            store_count += 1

                        if (count + 2) == npts:
                            self.sol_state[-1] = np.concatenate(
                                (
                                    np.reshape(x, objects_count * 3),
                                    np.reshape(v, objects_count * 3),
                                )
                            )

                        progress_bar.update(
                            task, completed=store_count, store_count=store_count + 1
                        )
                    progress_bar.update(
                        task, completed=store_npts, store_count=store_npts
                    )

                case "euler_cromer":
                    for count in range(npts - 1):
                        a = acceleration(objects_count, x, m, G)
                        x, v, x_err_comp_sum, v_err_comp_sum = self.euler_cromer(
                            x, v, a, dt, x_err_comp_sum, v_err_comp_sum
                        )
                        if (count + 1) % self.store_every_n == 0:
                            self.sol_state[store_count + 1] = np.concatenate(
                                (
                                    np.reshape(x, objects_count * 3),
                                    np.reshape(v, objects_count * 3),
                                )
                            )
                            store_count += 1

                        if (count + 2) == npts:
                            self.sol_state[-1] = np.concatenate(
                                (
                                    np.reshape(x, objects_count * 3),
                                    np.reshape(v, objects_count * 3),
                                )
                            )

                        progress_bar.update(
                            task, completed=store_count, store_count=store_count + 1
                        )
                    progress_bar.update(
                        task, completed=store_npts, store_count=store_npts
                    )

                case "rk4":
                    for count in range(npts - 1):
                        x, v, x_err_comp_sum, v_err_comp_sum = self.rk4(
                            objects_count,
                            x,
                            v,
                            m,
                            G,
                            dt,
                            x_err_comp_sum,
                            v_err_comp_sum,
                        )
                        if (count + 1) % self.store_every_n == 0:
                            self.sol_state[store_count + 1] = np.concatenate(
                                (
                                    np.reshape(x, objects_count * 3),
                                    np.reshape(v, objects_count * 3),
                                )
                            )
                            store_count += 1

                        if (count + 2) == npts:
                            self.sol_state[-1] = np.concatenate(
                                (
                                    np.reshape(x, objects_count * 3),
                                    np.reshape(v, objects_count * 3),
                                )
                            )

                        progress_bar.update(
                            task, completed=store_count, store_count=store_count + 1
                        )
                    progress_bar.update(
                        task, completed=store_npts, store_count=store_npts
                    )

                case "leapfrog":
                    a = acceleration(objects_count, x, m, G)
                    for count in range(npts - 1):
                        x, v, a, x_err_comp_sum, v_err_comp_sum = self.leapfrog(
                            objects_count,
                            x,
                            v,
                            a,
                            m,
                            dt,
                            G,
                            x_err_comp_sum,
                            v_err_comp_sum,
                        )
                        if (count + 1) % self.store_every_n == 0:
                            self.sol_state[store_count + 1] = np.concatenate(
                                (
                                    np.reshape(x, objects_count * 3),
                                    np.reshape(v, objects_count * 3),
                                )
                            )
                            store_count += 1

                        if (count + 2) == npts:
                            self.sol_state[-1] = np.concatenate(
                                (
                                    np.reshape(x, objects_count * 3),
                                    np.reshape(v, objects_count * 3),
                                )
                            )

                        progress_bar.update(
                            task, completed=store_count, store_count=store_count + 1
                        )
                    progress_bar.update(
                        task, completed=store_npts, store_count=store_npts
                    )

    @staticmethod
    def euler(x, v, a, dt, x_err_comp_sum, v_err_comp_sum):
        x0 = x.copy()
        v0 = v.copy()

        x_err_comp_sum += v * dt
        v_err_comp_sum += a * dt

        x = x0 + x_err_comp_sum
        v = v0 + v_err_comp_sum

        x_err_comp_sum += x0 - x
        v_err_comp_sum += v0 - v

        return x, v, x_err_comp_sum, v_err_comp_sum

    @staticmethod
    def euler_cromer(x, v, a, dt, x_err_comp_sum, v_err_comp_sum):
        v0 = v.copy()

        v_err_comp_sum += a * dt
        v = v0 + v_err_comp_sum
        v_err_comp_sum += v0 - v

        x0 = x.copy()

        x_err_comp_sum += v * dt
        x = x0 + x_err_comp_sum
        x_err_comp_sum += x0 - x

        return x, v, x_err_comp_sum, v_err_comp_sum

    @staticmethod
    def rk4(objects_count, x, v, m, G, dt, x_err_comp_sum, v_err_comp_sum):
        vk1 = acceleration(objects_count, x, m, G)
        xk1 = v

        vk2 = acceleration(objects_count, x + 0.5 * xk1 * dt, m, G)
        xk2 = v + 0.5 * vk1 * dt

        vk3 = acceleration(objects_count, x + 0.5 * xk2 * dt, m, G)
        xk3 = v + 0.5 * vk2 * dt

        vk4 = acceleration(objects_count, x + xk3 * dt, m, G)
        xk4 = v + vk3 * dt

        v0 = v.copy()
        x0 = x.copy()

        v_err_comp_sum += dt * (vk1 + 2 * vk2 + 2 * vk3 + vk4) / 6.0
        x_err_comp_sum += dt * (xk1 + 2 * xk2 + 2 * xk3 + xk4) / 6.0

        v = v0 + v_err_comp_sum
        x = x0 + x_err_comp_sum

        v_err_comp_sum += v0 - v
        x_err_comp_sum += x0 - x

        return x, v, x_err_comp_sum, v_err_comp_sum

    @staticmethod
    def leapfrog(objects_count, x, v, a, m, dt, G, x_err_comp_sum, v_err_comp_sum):
        v0 = v.copy()
        x0 = x.copy()

        a_0 = a
        x_err_comp_sum += v * dt + a_0 * 0.5 * dt * dt
        x = x0 + x_err_comp_sum
        x_err_comp_sum += x0 - x

        a_1 = acceleration(objects_count, x, m, G)
        v_err_comp_sum += (a_0 + a_1) * 0.5 * dt
        v = v0 + v_err_comp_sum
        v_err_comp_sum += v0 - v

        return x, v, a_1, x_err_comp_sum, v_err_comp_sum
