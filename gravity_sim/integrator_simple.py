"""
Simple integrators:
    Euler, Euler-Cromer, Runge-Kutta 4th order, Leapfrog
"""

import ctypes

from common import acceleration


class SimpleIntegrator:
    def simulation(self, simulator, integrator, objects_count, m, G, dt, time_speed):
        if simulator.is_c_lib == True:
            match integrator:
                case "euler":
                    if (
                        simulator.is_initialize == True
                        and simulator.is_initialize_integrator == "euler"
                    ):
                        simulator.is_initialize = False

                    simulator.c_lib.euler(
                        ctypes.c_int(objects_count),
                        simulator.x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        simulator.v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        ctypes.c_double(G),
                        ctypes.c_double(dt),
                        ctypes.c_int(time_speed),
                    )

                case "euler_cromer":
                    if (
                        simulator.is_initialize == True
                        and simulator.is_initialize_integrator == "euler_cromer"
                    ):
                        simulator.is_initialize = False

                    simulator.c_lib.euler_cromer(
                        ctypes.c_int(objects_count),
                        simulator.x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        simulator.v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        ctypes.c_double(G),
                        ctypes.c_double(dt),
                        ctypes.c_int(time_speed),
                    )

                case "rk4":
                    if (
                        simulator.is_initialize == True
                        and simulator.is_initialize_integrator == "rk4"
                    ):
                        simulator.is_initialize = False

                    simulator.c_lib.rk4(
                        ctypes.c_int(objects_count),
                        simulator.x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        simulator.v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        ctypes.c_double(G),
                        ctypes.c_double(dt),
                        ctypes.c_int(time_speed),
                    )

                case "leapfrog":
                    if (
                        simulator.is_initialize == True
                        and simulator.is_initialize_integrator == "leapfrog"
                    ):
                        simulator.c_lib.acceleration(
                            ctypes.c_int(objects_count),
                            simulator.x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            simulator.a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            ctypes.c_double(G),
                        )
                        simulator.is_initialize = False

                    simulator.c_lib.leapfrog(
                        ctypes.c_int(objects_count),
                        simulator.x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        simulator.v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        simulator.a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        ctypes.c_double(G),
                        ctypes.c_double(dt),
                        ctypes.c_int(time_speed),
                    )

        elif simulator.is_c_lib == False:
            match integrator:
                case "euler":
                    if (
                        simulator.is_initialize == True
                        and simulator.is_initialize_integrator == "euler"
                    ):
                        simulator.is_initialize = False

                    simulator.x, simulator.v = self._euler(
                        objects_count, simulator.x, simulator.v, m, G, dt, time_speed
                    )

                case "euler_cromer":
                    if (
                        simulator.is_initialize == True
                        and simulator.is_initialize_integrator == "euler_cromer"
                    ):
                        simulator.is_initialize = False

                    simulator.x, simulator.v = self._euler_cromer(
                        objects_count, simulator.x, simulator.v, m, G, dt, time_speed
                    )

                case "rk4":
                    if (
                        simulator.is_initialize == True
                        and simulator.is_initialize_integrator == "rk4"
                    ):
                        simulator.is_initialize = False

                    simulator.x, simulator.v = self._rk4(
                        objects_count,
                        simulator.x,
                        simulator.v,
                        m,
                        G,
                        dt,
                        time_speed,
                    )

                case "leapfrog":
                    if (
                        simulator.is_initialize == True
                        and simulator.is_initialize_integrator == "leapfrog"
                    ):
                        simulator.a = acceleration(objects_count, simulator.x, m, G)
                        simulator.is_initialize = False

                    simulator.x, simulator.v, simulator.a = self._leapfrog(
                        objects_count,
                        simulator.x,
                        simulator.v,
                        simulator.a,
                        m,
                        G,
                        dt,
                        time_speed,
                    )

    @staticmethod
    def _euler(objects_count, x, v, m, G, dt, time_speed):
        for _ in range(time_speed):
            a = acceleration(objects_count, x, m, G)
            x += v * dt
            v += a * dt

        return x, v

    @staticmethod
    def _euler_cromer(objects_count, x, v, m, G, dt, time_speed):
        for _ in range(time_speed):
            a = acceleration(objects_count, x, m, G)
            v += a * dt
            x += v * dt

        return x, v

    @staticmethod
    def _rk4(objects_count, x, v, m, G, dt, time_speed):
        for _ in range(time_speed):
            vk1 = acceleration(objects_count, x, m, G)
            xk1 = v

            vk2 = acceleration(objects_count, x + 0.5 * xk1 * dt, m, G)
            xk2 = v + 0.5 * vk1 * dt

            vk3 = acceleration(objects_count, x + 0.5 * xk2 * dt, m, G)
            xk3 = v + 0.5 * vk2 * dt

            vk4 = acceleration(objects_count, x + xk3 * dt, m, G)
            xk4 = v + vk3 * dt

            v = v + dt * (vk1 + 2 * vk2 + 2 * vk3 + vk4) / 6.0
            x = x + dt * (xk1 + 2 * xk2 + 2 * xk3 + xk4) / 6.0

        return x, v

    @staticmethod
    def _leapfrog(objects_count, x, v, a, m, G, dt, time_speed):
        a_1 = a
        for _ in range(time_speed):
            a_0 = a_1
            x = x + v * dt + a_0 * 0.5 * dt * dt
            a_1 = acceleration(objects_count, x, m, G)
            v = v + (a_0 + a_1) * 0.5 * dt

        return x, v, a_1
