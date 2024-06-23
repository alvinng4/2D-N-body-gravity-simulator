import csv
import ctypes
import math
from pathlib import Path
import timeit
import sys

import numpy as np

from progress_bar import Progress_bar
from integrator_fixed_step_size import FIXED_STEP_SIZE_INTEGRATOR
from integrator_rk_embedded import RK_EMBEDDED
from integrator_ias15 import IAS15


class Simulator:
    # Conversion factor from km^3 s^-2 to AU^3 d^-2
    CONVERSION_FACTOR = (86400**2) / (149597870.7**3)
    # GM values (km^3 s^-2)
    # ref: https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf
    GM_SI = {
        "Sun": 132712440041.279419,
        "Mercury": 22031.868551,
        "Venus": 324858.592000,
        "Earth": 398600.435507,
        "Mars": 42828.375816,
        "Jupiter": 126712764.100000,
        "Saturn": 37940584.841800,
        "Uranus": 5794556.400000,
        "Neptune": 6836527.100580,
        "Moon": 4902.800118,
        "Pluto": 975.500000,
        "Ceres": 62.62890,
        "Vesta": 17.288245,
    }
    # GM values (AU^3 d^-2)
    GM = {
        "Sun": 132712440041.279419 * CONVERSION_FACTOR,
        "Mercury": 22031.868551 * CONVERSION_FACTOR,
        "Venus": 324858.592000 * CONVERSION_FACTOR,
        "Earth": 398600.435507 * CONVERSION_FACTOR,
        "Mars": 42828.375816 * CONVERSION_FACTOR,
        "Jupiter": 126712764.100000 * CONVERSION_FACTOR,
        "Saturn": 37940584.841800 * CONVERSION_FACTOR,
        "Uranus": 5794556.400000 * CONVERSION_FACTOR,
        "Neptune": 6836527.100580 * CONVERSION_FACTOR,
        "Moon": 4902.800118 * CONVERSION_FACTOR,
        "Pluto": 975.500000 * CONVERSION_FACTOR,
        "Ceres": 62.62890 * CONVERSION_FACTOR,
        "Vesta": 17.288245 * CONVERSION_FACTOR,
    }
    # Solar system masses (M_sun^-1)
    SOLAR_SYSTEM_MASSES = {
        "Sun": 1.0,
        "Mercury": GM_SI["Mercury"] / GM_SI["Sun"],
        "Venus": GM_SI["Venus"] / GM_SI["Sun"],
        "Earth": GM_SI["Earth"] / GM_SI["Sun"],
        "Mars": GM_SI["Mars"] / GM_SI["Sun"],
        "Jupiter": GM_SI["Jupiter"] / GM_SI["Sun"],
        "Saturn": GM_SI["Saturn"] / GM_SI["Sun"],
        "Uranus": GM_SI["Uranus"] / GM_SI["Sun"],
        "Neptune": GM_SI["Neptune"] / GM_SI["Sun"],
        "Moon": GM_SI["Moon"] / GM_SI["Sun"],
        "Pluto": GM_SI["Pluto"] / GM_SI["Sun"],
        "Ceres": GM_SI["Ceres"] / GM_SI["Sun"],
        "Vesta": GM_SI["Vesta"] / GM_SI["Sun"],
    }
    # Gravitational constant (kg^-1 m^3 s^-2):
    CONSTANT_G_SI = 6.67430e-11
    # Gravitational constant (M_sun^-1 AU^3 d^-2):
    CONSTANT_G = GM["Sun"]

    # Solar system position and velocities data
    # Units: AU-D
    # Coordinate center: Solar System Barycenter
    # Data dated on A.D. 2024-Jan-01 00:00:00.0000 TDB
    # Computational data generated by NASA JPL Horizons System https://ssd.jpl.nasa.gov/horizons/
    SOLAR_SYSTEM_POS = {
        "Sun": [-7.967955691533730e-03, -2.906227441573178e-03, 2.103054301547123e-04],
        "Mercury": [
            -2.825983269538632e-01,
            1.974559795958082e-01,
            4.177433558063677e-02,
        ],
        "Venus": [
            -7.232103701666379e-01,
            -7.948302026312400e-02,
            4.042871428174315e-02,
        ],
        "Earth": [-1.738192017257054e-01, 9.663245550235138e-01, 1.553901854897183e-04],
        "Mars": [-3.013262392582653e-01, -1.454029331393295e00, -2.300531433991428e-02],
        "Jupiter": [3.485202469657674e00, 3.552136904413157e00, -9.271035442798399e-02],
        "Saturn": [8.988104223143450e00, -3.719064854634689e00, -2.931937777323593e-01],
        "Uranus": [1.226302417897505e01, 1.529738792480545e01, -1.020549026883563e-01],
        "Neptune": [
            2.983501460984741e01,
            -1.793812957956852e00,
            -6.506401132254588e-01,
        ],
        "Moon": [-1.762788124769829e-01, 9.674377513177153e-01, 3.236901585768862e-04],
        "Pluto": [1.720200478843485e01, -3.034155683573043e01, -1.729127607100611e00],
        "Ceres": [-1.103880510367569e00, -2.533340440444230e00, 1.220283937721780e-01],
        "Vesta": [-8.092549658731499e-02, 2.558381434460076e00, -6.695836142398572e-02],
    }
    SOLAR_SYSTEM_VEL = {
        "Sun": [4.875094764261564e-06, -7.057133213976680e-06, -4.573453713094512e-08],
        "Mercury": [
            -2.232165900189702e-02,
            -2.157207103176252e-02,
            2.855193410495743e-04,
        ],
        "Venus": [
            2.034068201002341e-03,
            -2.020828626592994e-02,
            -3.945639843855159e-04,
        ],
        "Earth": [
            -1.723001232538228e-02,
            -2.967721342618870e-03,
            6.382125383116755e-07,
        ],
        "Mars": [1.424832259345280e-02, -1.579236181580905e-03, -3.823722796161561e-04],
        "Jupiter": [
            -5.470970658852281e-03,
            5.642487338479145e-03,
            9.896190602066252e-05,
        ],
        "Saturn": [
            1.822013845554067e-03,
            5.143470425888054e-03,
            -1.617235904887937e-04,
        ],
        "Uranus": [
            -3.097615358317413e-03,
            2.276781932345769e-03,
            4.860433222241686e-05,
        ],
        "Neptune": [
            1.676536611817232e-04,
            3.152098732861913e-03,
            -6.877501095688201e-05,
        ],
        "Moon": [
            -1.746667306153906e-02,
            -3.473438277358121e-03,
            -3.359028758606074e-05,
        ],
        "Pluto": [2.802810313667557e-03, 8.492056438614633e-04, -9.060790113327894e-04],
        "Ceres": [
            8.978653480111301e-03,
            -4.873256528198994e-03,
            -1.807162046049230e-03,
        ],
        "Vesta": [
            -1.017876585480054e-02,
            -5.452367109338154e-04,
            1.255870551153315e-03,
        ],
    }

    def __init__(self, grav_plot):
        self.is_c_lib = grav_plot.is_c_lib
        if grav_plot.is_c_lib:
            self.c_lib = grav_plot.c_lib

        self.is_exit = grav_plot.is_exit
        self.x = np.zeros(0)
        self.v = np.zeros(0)
        self.m = np.zeros(0)
        self.objects_count = 0
        self.G = 0.0

    def initialize_system_numpy(self, grav_plot):
        self.G = self.CONSTANT_G

        # Read information of the customized system
        if self.system not in grav_plot.default_systems:
            file_path = Path(__file__).parent / "customized_systems.csv"
            try:
                with open(file_path, "r") as file:
                    reader = csv.reader(file)
                    for row in reader:
                        if self.system == row[0]:
                            self.objects_count = int(row[1])
                            self.m = np.zeros(self.objects_count)
                            for i in range(self.objects_count):
                                self.m[i] = row[2 + i]

                            self.x = np.zeros((self.objects_count, 3))
                            self.v = np.zeros((self.objects_count, 3))
                            for i in range(self.objects_count):
                                for j in range(3):
                                    self.x[i][j] = row[
                                        2 + self.objects_count + i * 3 + j
                                    ]
                                    self.v[i][j] = row[
                                        2
                                        + self.objects_count
                                        + self.objects_count * 3
                                        + i * 3
                                        + j
                                    ]

            except FileNotFoundError:
                sys.exit(
                    "Error: customized_systems.csv not found in gravity_plot. Terminating program."
                )

        else:
            # Pre-defined systems
            match self.system:
                case "circular_binary_orbit":
                    R1 = np.array([1.0, 0.0, 0.0])
                    R2 = np.array([-1.0, 0.0, 0.0])
                    V1 = np.array([0.0, 0.5, 0.0])
                    V2 = np.array([0.0, -0.5, 0.0])
                    self.x = np.array([R1, R2])
                    self.v = np.array([V1, V2])
                    self.m = np.array([1.0 / self.G, 1.0 / self.G])
                    self.objects_count = 2

                case "eccentric_binary_orbit":
                    R1 = np.array([1.0, 0.0, 0.0])
                    R2 = np.array([-1.25, 0.0, 0.0])
                    V1 = np.array([0.0, 0.5, 0.0])
                    V2 = np.array([0.0, -0.625, 0.0])
                    self.x = np.array([R1, R2])
                    self.v = np.array([V1, V2])
                    self.m = np.array([1.0 / self.G, 0.8 / self.G])
                    self.objects_count = 2

                case "3d_helix":
                    R1 = np.array([0.0, 0.0, -1.0])
                    R2 = np.array([-math.sqrt(3.0) / 2.0, 0.0, 0.5])
                    R3 = np.array([math.sqrt(3.0) / 2.0, 0.0, 0.5])
                    v0 = math.sqrt(1.0 / math.sqrt(3))
                    V1 = np.array([-v0, 0.5, 0.0])
                    V2 = np.array([0.5 * v0, 0.5, (math.sqrt(3.0) / 2.0) * v0])
                    V3 = np.array([0.5 * v0, 0.5, -(math.sqrt(3.0) / 2.0) * v0])
                    self.x = np.array([R1, R2, R3])
                    self.v = np.array([V1, V2, V3])
                    self.m = np.array([1.0 / self.G, 1.0 / self.G, 1.0 / self.G])
                    self.objects_count = 3

                case "sun_earth_moon":
                    self.m = np.array(
                        [
                            self.SOLAR_SYSTEM_MASSES["Sun"],
                            self.SOLAR_SYSTEM_MASSES["Earth"],
                            self.SOLAR_SYSTEM_MASSES["Moon"],
                        ]
                    )
                    R_CM = (
                        1
                        / np.sum(self.m)
                        * (
                            self.m[0] * np.array(self.SOLAR_SYSTEM_POS["Sun"])
                            + self.m[1] * np.array(self.SOLAR_SYSTEM_POS["Earth"])
                            + self.m[2] * np.array(self.SOLAR_SYSTEM_POS["Moon"])
                        )
                    )
                    V_CM = (
                        1
                        / np.sum(self.m)
                        * (
                            self.m[0] * np.array(self.SOLAR_SYSTEM_VEL["Sun"])
                            + self.m[1] * np.array(self.SOLAR_SYSTEM_VEL["Earth"])
                            + self.m[2] * np.array(self.SOLAR_SYSTEM_VEL["Moon"])
                        )
                    )
                    R1 = np.array(self.SOLAR_SYSTEM_POS["Sun"] - R_CM)
                    R2 = np.array(self.SOLAR_SYSTEM_POS["Earth"] - R_CM)
                    R3 = np.array(self.SOLAR_SYSTEM_POS["Moon"] - R_CM)
                    V1 = np.array(self.SOLAR_SYSTEM_VEL["Sun"] - V_CM)
                    V2 = np.array(self.SOLAR_SYSTEM_VEL["Earth"] - V_CM)
                    V3 = np.array(self.SOLAR_SYSTEM_VEL["Moon"] - V_CM)
                    self.x = np.array([R1, R2, R3])
                    self.v = np.array([V1, V2, V3])
                    self.objects_count = 3

                case "figure-8":
                    R1 = np.array([0.970043, -0.24308753, 0.0])
                    R2 = np.array([-0.970043, 0.24308753, 0.0])
                    R3 = np.array([0.0, 0.0, 0.0])
                    V1 = np.array([0.466203685, 0.43236573, 0.0])
                    V2 = np.array([0.466203685, 0.43236573, 0.0])
                    V3 = np.array([-0.93240737, -0.86473146, 0.0])
                    self.x = np.array([R1, R2, R3])
                    self.v = np.array([V1, V2, V3])
                    self.m = np.array([1.0 / self.G, 1.0 / self.G, 1.0 / self.G])
                    self.objects_count = 3

                case "pyth-3-body":
                    R1 = np.array([1.0, 3.0, 0.0])
                    R2 = np.array([-2.0, -1.0, 0.0])
                    R3 = np.array([1.0, -1.0, 0.0])
                    V1 = np.array([0.0, 0.0, 0.0])
                    V2 = np.array([0.0, 0.0, 0.0])
                    V3 = np.array([0.0, 0.0, 0.0])
                    self.x = np.array([R1, R2, R3])
                    self.v = np.array([V1, V2, V3])
                    self.m = np.array([3.0 / self.G, 4.0 / self.G, 5.0 / self.G])
                    self.objects_count = 3

                case "solar_system":
                    self.m = np.array(
                        [
                            self.SOLAR_SYSTEM_MASSES["Sun"],
                            self.SOLAR_SYSTEM_MASSES["Mercury"],
                            self.SOLAR_SYSTEM_MASSES["Venus"],
                            self.SOLAR_SYSTEM_MASSES["Earth"],
                            self.SOLAR_SYSTEM_MASSES["Mars"],
                            self.SOLAR_SYSTEM_MASSES["Jupiter"],
                            self.SOLAR_SYSTEM_MASSES["Saturn"],
                            self.SOLAR_SYSTEM_MASSES["Uranus"],
                            self.SOLAR_SYSTEM_MASSES["Neptune"],
                        ]
                    )
                    R_CM = (
                        1
                        / np.sum(self.m)
                        * (
                            self.m[0] * np.array(self.SOLAR_SYSTEM_POS["Sun"])
                            + self.m[1] * np.array(self.SOLAR_SYSTEM_POS["Mercury"])
                            + self.m[2] * np.array(self.SOLAR_SYSTEM_POS["Venus"])
                            + self.m[3] * np.array(self.SOLAR_SYSTEM_POS["Earth"])
                            + self.m[4] * np.array(self.SOLAR_SYSTEM_POS["Mars"])
                            + self.m[5] * np.array(self.SOLAR_SYSTEM_POS["Jupiter"])
                            + self.m[6] * np.array(self.SOLAR_SYSTEM_POS["Saturn"])
                            + self.m[7] * np.array(self.SOLAR_SYSTEM_POS["Uranus"])
                            + self.m[8] * np.array(self.SOLAR_SYSTEM_POS["Neptune"])
                        )
                    )
                    V_CM = (
                        1
                        / np.sum(self.m)
                        * (
                            self.m[0] * np.array(self.SOLAR_SYSTEM_VEL["Sun"])
                            + self.m[1] * np.array(self.SOLAR_SYSTEM_VEL["Mercury"])
                            + self.m[2] * np.array(self.SOLAR_SYSTEM_VEL["Venus"])
                            + self.m[3] * np.array(self.SOLAR_SYSTEM_VEL["Earth"])
                            + self.m[4] * np.array(self.SOLAR_SYSTEM_VEL["Mars"])
                            + self.m[5] * np.array(self.SOLAR_SYSTEM_VEL["Jupiter"])
                            + self.m[6] * np.array(self.SOLAR_SYSTEM_VEL["Saturn"])
                            + self.m[7] * np.array(self.SOLAR_SYSTEM_VEL["Uranus"])
                            + self.m[8] * np.array(self.SOLAR_SYSTEM_VEL["Neptune"])
                        )
                    )

                    R1 = np.array(self.SOLAR_SYSTEM_POS["Sun"] - R_CM)
                    R2 = np.array(self.SOLAR_SYSTEM_POS["Mercury"] - R_CM)
                    R3 = np.array(self.SOLAR_SYSTEM_POS["Venus"] - R_CM)
                    R4 = np.array(self.SOLAR_SYSTEM_POS["Earth"] - R_CM)
                    R5 = np.array(self.SOLAR_SYSTEM_POS["Mars"] - R_CM)
                    R6 = np.array(self.SOLAR_SYSTEM_POS["Jupiter"] - R_CM)
                    R7 = np.array(self.SOLAR_SYSTEM_POS["Saturn"] - R_CM)
                    R8 = np.array(self.SOLAR_SYSTEM_POS["Uranus"] - R_CM)
                    R9 = np.array(self.SOLAR_SYSTEM_POS["Neptune"] - R_CM)

                    V1 = np.array(self.SOLAR_SYSTEM_VEL["Sun"] - V_CM)
                    V2 = np.array(self.SOLAR_SYSTEM_VEL["Mercury"] - V_CM)
                    V3 = np.array(self.SOLAR_SYSTEM_VEL["Venus"] - V_CM)
                    V4 = np.array(self.SOLAR_SYSTEM_VEL["Earth"] - V_CM)
                    V5 = np.array(self.SOLAR_SYSTEM_VEL["Mars"] - V_CM)
                    V6 = np.array(self.SOLAR_SYSTEM_VEL["Jupiter"] - V_CM)
                    V7 = np.array(self.SOLAR_SYSTEM_VEL["Saturn"] - V_CM)
                    V8 = np.array(self.SOLAR_SYSTEM_VEL["Uranus"] - V_CM)
                    V9 = np.array(self.SOLAR_SYSTEM_VEL["Neptune"] - V_CM)

                    self.x = np.array([R1, R2, R3, R4, R5, R6, R7, R8, R9])
                    self.v = np.array([V1, V2, V3, V4, V5, V6, V7, V8, V9])
                    self.objects_count = 9

                case "solar_system_plus":
                    self.m = np.array(
                        [
                            self.SOLAR_SYSTEM_MASSES["Sun"],
                            self.SOLAR_SYSTEM_MASSES["Mercury"],
                            self.SOLAR_SYSTEM_MASSES["Venus"],
                            self.SOLAR_SYSTEM_MASSES["Earth"],
                            self.SOLAR_SYSTEM_MASSES["Mars"],
                            self.SOLAR_SYSTEM_MASSES["Jupiter"],
                            self.SOLAR_SYSTEM_MASSES["Saturn"],
                            self.SOLAR_SYSTEM_MASSES["Uranus"],
                            self.SOLAR_SYSTEM_MASSES["Neptune"],
                            self.SOLAR_SYSTEM_MASSES["Pluto"],
                            self.SOLAR_SYSTEM_MASSES["Ceres"],
                            self.SOLAR_SYSTEM_MASSES["Vesta"],
                        ]
                    )

                    R_CM = (
                        1
                        / np.sum(self.m)
                        * (
                            self.m[0] * np.array(self.SOLAR_SYSTEM_POS["Sun"])
                            + self.m[1] * np.array(self.SOLAR_SYSTEM_POS["Mercury"])
                            + self.m[2] * np.array(self.SOLAR_SYSTEM_POS["Venus"])
                            + self.m[3] * np.array(self.SOLAR_SYSTEM_POS["Earth"])
                            + self.m[4] * np.array(self.SOLAR_SYSTEM_POS["Mars"])
                            + self.m[5] * np.array(self.SOLAR_SYSTEM_POS["Jupiter"])
                            + self.m[6] * np.array(self.SOLAR_SYSTEM_POS["Saturn"])
                            + self.m[7] * np.array(self.SOLAR_SYSTEM_POS["Uranus"])
                            + self.m[8] * np.array(self.SOLAR_SYSTEM_POS["Neptune"])
                            + self.m[9] * np.array(self.SOLAR_SYSTEM_POS["Pluto"])
                            + self.m[10] * np.array(self.SOLAR_SYSTEM_POS["Ceres"])
                            + self.m[11] * np.array(self.SOLAR_SYSTEM_POS["Vesta"])
                        )
                    )

                    V_CM = (
                        1
                        / np.sum(self.m)
                        * (
                            self.m[0] * np.array(self.SOLAR_SYSTEM_VEL["Sun"])
                            + self.m[1] * np.array(self.SOLAR_SYSTEM_VEL["Mercury"])
                            + self.m[2] * np.array(self.SOLAR_SYSTEM_VEL["Venus"])
                            + self.m[3] * np.array(self.SOLAR_SYSTEM_VEL["Earth"])
                            + self.m[4] * np.array(self.SOLAR_SYSTEM_VEL["Mars"])
                            + self.m[5] * np.array(self.SOLAR_SYSTEM_VEL["Jupiter"])
                            + self.m[6] * np.array(self.SOLAR_SYSTEM_VEL["Saturn"])
                            + self.m[7] * np.array(self.SOLAR_SYSTEM_VEL["Uranus"])
                            + self.m[8] * np.array(self.SOLAR_SYSTEM_VEL["Neptune"])
                            + self.m[9] * np.array(self.SOLAR_SYSTEM_VEL["Pluto"])
                            + self.m[10] * np.array(self.SOLAR_SYSTEM_VEL["Ceres"])
                            + self.m[11] * np.array(self.SOLAR_SYSTEM_VEL["Vesta"])
                        )
                    )

                    R1 = np.array(self.SOLAR_SYSTEM_POS["Sun"] - R_CM)
                    R2 = np.array(self.SOLAR_SYSTEM_POS["Mercury"] - R_CM)
                    R3 = np.array(self.SOLAR_SYSTEM_POS["Venus"] - R_CM)
                    R4 = np.array(self.SOLAR_SYSTEM_POS["Earth"] - R_CM)
                    R5 = np.array(self.SOLAR_SYSTEM_POS["Mars"] - R_CM)
                    R6 = np.array(self.SOLAR_SYSTEM_POS["Jupiter"] - R_CM)
                    R7 = np.array(self.SOLAR_SYSTEM_POS["Saturn"] - R_CM)
                    R8 = np.array(self.SOLAR_SYSTEM_POS["Uranus"] - R_CM)
                    R9 = np.array(self.SOLAR_SYSTEM_POS["Neptune"] - R_CM)
                    R10 = np.array(self.SOLAR_SYSTEM_POS["Pluto"] - R_CM)
                    R11 = np.array(self.SOLAR_SYSTEM_POS["Ceres"] - R_CM)
                    R12 = np.array(self.SOLAR_SYSTEM_POS["Vesta"] - R_CM)

                    V1 = np.array(self.SOLAR_SYSTEM_VEL["Sun"] - V_CM)
                    V2 = np.array(self.SOLAR_SYSTEM_VEL["Mercury"] - V_CM)
                    V3 = np.array(self.SOLAR_SYSTEM_VEL["Venus"] - V_CM)
                    V4 = np.array(self.SOLAR_SYSTEM_VEL["Earth"] - V_CM)
                    V5 = np.array(self.SOLAR_SYSTEM_VEL["Mars"] - V_CM)
                    V6 = np.array(self.SOLAR_SYSTEM_VEL["Jupiter"] - V_CM)
                    V7 = np.array(self.SOLAR_SYSTEM_VEL["Saturn"] - V_CM)
                    V8 = np.array(self.SOLAR_SYSTEM_VEL["Uranus"] - V_CM)
                    V9 = np.array(self.SOLAR_SYSTEM_VEL["Neptune"] - V_CM)
                    V10 = np.array(self.SOLAR_SYSTEM_VEL["Pluto"] - V_CM)
                    V11 = np.array(self.SOLAR_SYSTEM_VEL["Ceres"] - V_CM)
                    V12 = np.array(self.SOLAR_SYSTEM_VEL["Vesta"] - V_CM)

                    self.x = np.array(
                        [R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12]
                    )
                    self.v = np.array(
                        [V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12]
                    )
                    self.objects_count = 12

    def simulation(self, is_custom_sys):
        print("Simulating the system...")
        start = timeit.default_timer()

        if is_custom_sys:
            self.x.resize((self.objects_count * 3,))
            self.v.resize((self.objects_count * 3,))

        match self.integrator:
            case "euler" | "euler_cromer" | "rk4" | "leapfrog":
                integrator = FIXED_STEP_SIZE_INTEGRATOR(self)

            case "rkf45" | "dopri" | "dverk" | "rkf78":
                integrator = RK_EMBEDDED(self)

            case "ias15":
                integrator = IAS15(self)

        self.sol_state = integrator.sol_state
        self.sol_time = integrator.sol_time
        self.sol_dt = integrator.sol_dt

        if self.is_c_lib:
            self.m = integrator.m
            self.G = integrator.G
            self.objects_count = integrator.objects_count

        stop = timeit.default_timer()
        print(f"Run time: {stop - start:.3f} s")
        print("")

    def compute_energy(self):
        """
        Compute the total energy using the sol_state array
        """
        print("Computing energy...")
        npts = len(self.sol_state)
        self.energy = np.zeros(npts)

        progress_bar = Progress_bar()

        start = timeit.default_timer()
        if self.is_c_lib == True:
            count = ctypes.c_int(0)
            with progress_bar:
                task = progress_bar.add_task("", total=npts)
                while count.value < npts:
                    self.c_lib.compute_energy(
                        ctypes.c_int(self.objects_count),
                        ctypes.c_int(npts),
                        ctypes.byref(count),
                        self.energy.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        self.sol_state.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        self.m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        ctypes.c_double(self.G),
                    )
                    progress_bar.update(task, completed=count.value)

        elif self.is_c_lib == False:
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
