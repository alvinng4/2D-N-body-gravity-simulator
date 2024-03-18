import numba as nb  # Note: nb.njit cannot works on functions inside a class
import numpy as np

from grav_obj import Grav_obj

from integrator_fixed_step_size import FIXED_STEP_SIZE_INTEGRATOR
from integrator_rk_embedded import RK_EMBEDDED
from integrator_ias15 import IAS15

from common import acceleration

class Simulator:
    def __init__(self, grav_sim):
        self.stats = grav_sim.stats
        self.settings = grav_sim.settings

        self.m = []
        self.x = []
        self.v = []
        self.a = []

        self.fixed_step_size_integrator = FIXED_STEP_SIZE_INTEGRATOR()
        self.rk_embdded_integrator = RK_EMBEDDED()
        self.ias15_integrator = IAS15()

        self.is_initialize = True
        self.set_all_integrators_false()
        self.is_rk4 = True  # Default integrator
        self.current_integrator = "rk4"
        self.is_initialize_integrator = "rk4"

    def run_simulation(self, grav_sim):
        if self.is_initialize == True:
            self.initialize_problem(grav_sim)

        # Simple euler is enough when there is no interaction
        if self.stats.objects_count == 1:
            self.fixed_step_size_integrator.simulation(
                self,
                "euler",
                self.stats.objects_count,
                self.m,
                Grav_obj.G,
                self.settings.dt,
                self.settings.time_speed,
            )
            self.stats.simulation_time += (
                self.settings.dt * self.settings.time_speed
            )
        else:
            match self.current_integrator:
                # Fixed step size integrators
                case "euler" | "euler_cromer" | "rk4" | "leapfrog":
                    self.fixed_step_size_integrator.simulation(
                        self,
                        self.current_integrator,
                        self.stats.objects_count,
                        self.m,
                        Grav_obj.G,
                        self.settings.dt,
                        self.settings.time_speed,
                        )
                    
                    self.stats.simulation_time += (
                        self.settings.dt * self.settings.time_speed
                    )
                # Embedded RK methods
                case "rkf45" | "dopri" | "dverk" | "rkf78":
                    self.rk_embdded_integrator.simulation(
                        self,
                        self.stats.objects_count,
                        self.m,
                        Grav_obj.G,
                        self.settings.tolerance,
                        self.settings.tolerance,
                        self.settings.expected_time_scale,
                        self.settings.rk_max_iteration,
                        self.settings.rk_min_iteration,
                    )

                case "ias15":
                    self.ias15_integrator.simulation(
                        self,
                        self.stats.objects_count,
                        self.m,
                        Grav_obj.G,
                        self.settings.tolerance,
                        self.settings.expected_time_scale,
                        self.settings.rk_max_iteration,
                        self.settings.rk_min_iteration,
                    )

        self.stats.total_energy = total_energy(
            self.stats.objects_count, self.x, self.v, self.m, Grav_obj.G
        )

    def initialize_problem(self, grav_sim):
        """
        Initialize x, v and m
        """
        objects_count = grav_sim.stats.objects_count
        self.x = np.zeros((objects_count, 3))
        self.v = np.zeros((objects_count, 3))
        self.m = np.zeros(objects_count)
        for j in range(objects_count):
            self.x[j] = np.array(
                [grav_sim.grav_objs.sprites()[j].params[f"r{i + 1}"] for i in range(3)]
            )
            self.v[j] = np.array(
                [grav_sim.grav_objs.sprites()[j].params[f"v{i + 1}"] for i in range(3)]
            )
            self.m[j] = grav_sim.grav_objs.sprites()[j].params["m"]

    def unload_value(self, grav_sim):
        """
        Unload the position and velocity values back the to main system
        """
        for j in range(self.stats.objects_count):
            grav_sim.grav_objs.sprites()[j].params["r1"] = self.x[j][0]
            grav_sim.grav_objs.sprites()[j].params["r2"] = self.x[j][1]
            grav_sim.grav_objs.sprites()[j].params["r3"] = self.x[j][2]
            grav_sim.grav_objs.sprites()[j].params["v1"] = self.v[j][0]
            grav_sim.grav_objs.sprites()[j].params["v2"] = self.v[j][1]
            grav_sim.grav_objs.sprites()[j].params["v3"] = self.v[j][2]

    def set_all_integrators_false(self):
        self.is_euler = False
        self.is_euler_cromer = False
        self.is_rk4 = False
        self.is_leapfrog = False
        self.is_rkf45 = False
        self.is_dopri = False
        self.is_dverk = False
        self.is_rkf78 = False
        self.is_ias15 = False

    def check_current_integrator(self):
        """
        Check what integrators are currently chosen
        """
        if self.is_euler == True:
            self.current_integrator = "euler"
        elif self.is_euler_cromer == True:
            self.current_integrator = "euler_cromer"
        elif self.is_rk4 == True:
            self.current_integrator = "rk4"
        elif self.is_leapfrog == True:
            self.current_integrator = "leapfrog"
        elif self.is_rkf45 == True:
            self.current_integrator = "rkf45"
        elif self.is_dopri == True:
            self.current_integrator = "dopri"
        elif self.is_dverk == True:
            self.current_integrator = "dverk"
        elif self.is_rkf78 == True:
            self.current_integrator = "rkf78"
        elif self.is_ias15 == True:
            self.current_integrator = "ias15"


@nb.njit
def _rk_embedded_initial_time_step(
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

    Modified: Return dt * 1e-3 since this function gives initial dt thats too large

    Reference: Moving Planets Around: An Introduction to N-Body Simulations Applied to Exoplanetary Systems
    Chapter 6, Page 92 - 94
    """
    tolerance_scale_x = abs_tolerance + rel_tolerance * np.abs(x)
    tolerance_scale_v = abs_tolerance + rel_tolerance * np.abs(v)
    sum_0 = np.sum(np.square(x / tolerance_scale_x)) + np.sum(
        np.square(v / tolerance_scale_v)
    )
    sum_1 = np.sum(np.square(v / tolerance_scale_x)) + np.sum(
        np.square(a / tolerance_scale_v)
    )
    d_0 = np.sqrt(sum_0 / (objects_count * 3 * 2))
    d_1 = np.sqrt(sum_1 / (objects_count * 3 * 2))

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
    d_2 = np.sqrt(sum_2 / (objects_count * 3 * 2)) / dt_0

    if max(d_1, d_2) <= 1e-15:
        dt_1 = max(1e-6, dt_0 * 1e-3)
    else:
        dt_1 = (0.01 / max(d_1, d_2)) ** (1.0 / (1.0 + power))
    dt = min(100 * dt_0, dt_1)

    return dt * 1e-3


@nb.njit
def rk_embedded(
    objects_count: int,
    x,
    v,
    m,
    G,
    expected_time_scale: float,
    simulation_time,
    actual_dt,
    power,
    power_test,
    coeff,
    weights,
    weights_test,
    max_iteration: int,
    min_iteration: int,
    abs_tolerance: float,
    rel_tolerance: float,
):
    # Initializing
    t = simulation_time
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

    for i in range(max_iteration):
        # Calculate xk and vk
        vk[0] = acceleration(objects_count, x, m, G)
        xk[0] = np.copy(v)
        for stage in range(1, stages):
            temp_v = np.zeros((objects_count, 3))
            temp_x = np.zeros((objects_count, 3))
            for j in range(stage):
                temp_v += coeff[stage - 1][j] * vk[j]
                temp_x += coeff[stage - 1][j] * xk[j]
            vk[stage] = acceleration(objects_count, x + actual_dt * temp_x, m, G)
            xk[stage] = v + actual_dt * temp_v

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
        v_1 = v + actual_dt * temp_v
        x_1 = x + actual_dt * temp_x
        error_estimation_delta_v *= actual_dt
        error_estimation_delta_x *= actual_dt

        # Error calculation
        tolerance_scale_v = (
            abs_tolerance + np.maximum(np.abs(v), np.abs(v_1)) * rel_tolerance
        )
        tolerance_scale_x = (
            abs_tolerance + np.maximum(np.abs(x), np.abs(x_1)) * rel_tolerance
        )

        # Sum up all the elements of x/tol and v/tol, square and divide by the total number of elements
        sum = np.sum(np.square(error_estimation_delta_x / tolerance_scale_x)) + np.sum(
            np.square(error_estimation_delta_v / tolerance_scale_v)
        )
        error = np.sqrt(sum / (objects_count * 3 * 2))

        if error <= 1 or actual_dt == expected_time_scale * 1e-12:
            t += actual_dt
            x = x_1
            v = v_1

        if error == 0.0: # Prevent extreme cases where the error is smaller than machine zero
            dt_new = actual_dt
        else:    
            dt_new = actual_dt * safety_fac / error ** (1.0 / (1.0 + min_power))
        # Prevent dt to be too small or too large relative to the last time step
        if dt_new > safety_fac_max * actual_dt:
            actual_dt = safety_fac_max * actual_dt
        elif dt_new < safety_fac_min * actual_dt:
            actual_dt = safety_fac_min * actual_dt
        else:
            actual_dt = dt_new

        if dt_new / expected_time_scale < 1e-12:
            actual_dt = expected_time_scale * 1e-12

        if i >= min_iteration and t > (simulation_time + expected_time_scale * 1e-5):
            return x, v, t, actual_dt

    # Return values once it reaches max iterations
    return x, v, t, actual_dt


@nb.njit
def total_energy(objects_count, x, v, m, G):
    E = 0
    for j in range(0, objects_count):
        E += 0.5 * m[j] * np.linalg.norm(v[j]) ** 2
        for k in range(0, objects_count):
            if j < k:
                R = x[j] - x[k]
                norm = np.linalg.norm(R)
                if norm != 0:
                    E -= G * m[j] * m[k] / norm 
                else:
                    return np.nan
    return E


@nb.njit
def rk_embedded_butcher_tableaus(order):
    """
    Butcher tableaus for embedded rk

    Reference: Moving Planets Around: An Introduction to N-Body Simulations Applied to Exoplanetary Systems
    Chapter 6, Page 100 - 101

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
                    [-8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0],
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
                    [2.0 / 27.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
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
