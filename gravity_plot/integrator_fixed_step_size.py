"""
Fixed step size integrator:
Euler, Euler-Cromer, Runge-Kutta 4th order, Leapfrog
"""

import ctypes
import math

import numpy as np

from common import acceleration
from progress_bar import Progress_bar


class FIXED_STEP_SIZE_INTEGRATOR:
    def __init__(
        self, simulator, integrator, objects_count, x, v, m, G, dt, tf, is_c_lib
    ):
        self.store_every_n = simulator.store_every_n
        if is_c_lib == True:
            self.c_lib = simulator.c_lib
            self.simulation_c_lib(integrator, objects_count, x, v, m, G, dt, tf)
        elif is_c_lib == False:
            self.simulation_numpy(integrator, objects_count, x, v, m, G, dt, tf)

    def simulation_c_lib(self, integrator, objects_count, x, v, m, G, dt, tf):
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

        progress_bar = Progress_bar()
        count = ctypes.c_int64(0)
        store_count = ctypes.c_int(0)
        with progress_bar:
            task = progress_bar.add_task("", total=store_npts)
            match integrator:
                case "euler":
                    while count.value < npts:
                        self.c_lib.euler(
                            ctypes.c_int(objects_count),
                            x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            ctypes.c_double(dt),
                            ctypes.c_int64(npts),
                            m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            ctypes.c_double(G),
                            self.sol_state.ctypes.data_as(
                                ctypes.POINTER(ctypes.c_double)
                            ),
                            ctypes.c_int(self.store_every_n),
                            ctypes.c_int(store_npts),
                            ctypes.byref(count),
                            ctypes.byref(store_count),
                            x_err_comp_sum.ctypes.data_as(
                                ctypes.POINTER(ctypes.c_double)
                            ),
                            v_err_comp_sum.ctypes.data_as(
                                ctypes.POINTER(ctypes.c_double)
                            ),
                        )
                        progress_bar.update(task, completed=store_count.value)
                    progress_bar.update(task, completed=store_npts)
                case "euler_cromer":
                    while count.value < npts:
                        self.c_lib.euler_cromer(
                            ctypes.c_int(objects_count),
                            x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            ctypes.c_double(dt),
                            ctypes.c_int64(npts),
                            m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            ctypes.c_double(G),
                            self.sol_state.ctypes.data_as(
                                ctypes.POINTER(ctypes.c_double)
                            ),
                            ctypes.c_int(self.store_every_n),
                            ctypes.c_int(store_npts),
                            ctypes.byref(count),
                            ctypes.byref(store_count),
                            x_err_comp_sum.ctypes.data_as(
                                ctypes.POINTER(ctypes.c_double)
                            ),
                            v_err_comp_sum.ctypes.data_as(
                                ctypes.POINTER(ctypes.c_double)
                            ),
                        )
                        progress_bar.update(task, completed=store_count.value)
                    progress_bar.update(task, completed=store_npts)
                case "rk4":
                    while count.value < npts:
                        self.c_lib.rk4(
                            ctypes.c_int(objects_count),
                            x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            ctypes.c_double(dt),
                            ctypes.c_int64(npts),
                            m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            ctypes.c_double(G),
                            self.sol_state.ctypes.data_as(
                                ctypes.POINTER(ctypes.c_double)
                            ),
                            ctypes.c_int(self.store_every_n),
                            ctypes.c_int(store_npts),
                            ctypes.byref(count),
                            ctypes.byref(store_count),
                            x_err_comp_sum.ctypes.data_as(
                                ctypes.POINTER(ctypes.c_double)
                            ),
                            v_err_comp_sum.ctypes.data_as(
                                ctypes.POINTER(ctypes.c_double)
                            ),
                        )
                        progress_bar.update(task, completed=store_count.value)
                    progress_bar.update(task, completed=store_npts)
                case "leapfrog":
                    while count.value < npts:
                        self.c_lib.leapfrog(
                            ctypes.c_int(objects_count),
                            x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            ctypes.c_double(dt),
                            ctypes.c_int64(npts),
                            m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            ctypes.c_double(G),
                            self.sol_state.ctypes.data_as(
                                ctypes.POINTER(ctypes.c_double)
                            ),
                            ctypes.c_int(self.store_every_n),
                            ctypes.c_int(store_npts),
                            ctypes.byref(count),
                            ctypes.byref(store_count),
                            x_err_comp_sum.ctypes.data_as(
                                ctypes.POINTER(ctypes.c_double)
                            ),
                            v_err_comp_sum.ctypes.data_as(
                                ctypes.POINTER(ctypes.c_double)
                            ),
                        )
                        progress_bar.update(task, completed=store_count.value)
                    progress_bar.update(task, completed=store_npts)

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

        progress_bar = Progress_bar()
        store_count = 0
        with progress_bar:
            match integrator:
                case "euler":
                    for count in progress_bar.track(range(npts - 1)):
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

                case "euler_cromer":
                    for count in progress_bar.track(range(npts - 1)):
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

                case "rk4":
                    for count in progress_bar.track(range(npts - 1)):
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

                case "leapfrog":
                    a = acceleration(objects_count, x, m, G)
                    for count in progress_bar.track(range(npts - 1)):
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
