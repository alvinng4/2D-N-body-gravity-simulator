import ctypes

import numpy as np

from common import compute_energy
from grav_obj import Grav_obj
from integrator_simple import SimpleIntegrator
from integrator_rk_embedded import RKEmbedded
from integrator_ias15 import IAS15


class Simulator:
    def __init__(self, grav_sim):
        self.is_c_lib = grav_sim.is_c_lib
        if self.is_c_lib == True:
            self.c_lib = grav_sim.c_lib

        self.stats = grav_sim.stats
        self.settings = grav_sim.settings

        self.m = np.array([])
        self.x = np.array([])
        self.v = np.array([])
        self.a = np.array([])

        self.simple_integrator = SimpleIntegrator()
        self.rk_embedded_integrator = RKEmbedded()
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
            self.simple_integrator.simulation(
                self,
                "euler",
                self.stats.objects_count,
                self.m,
                Grav_obj.G,
                self.settings.dt,
                self.settings.time_speed,
            )
            self.stats.simulation_time += self.settings.dt * self.settings.time_speed
        else:
            match self.current_integrator:
                # Simple integrators
                case "euler" | "euler_cromer" | "rk4" | "leapfrog":
                    self.simple_integrator.simulation(
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
                    self.rk_embedded_integrator.simulation(
                        self,
                        self.stats.objects_count,
                        self.m,
                        Grav_obj.G,
                        self.settings.tolerance,
                        self.settings.tolerance,
                        self.settings.expected_time_scale,
                        self.settings.max_iteration,
                        self.settings.min_iteration,
                    )

                case "ias15":
                    self.ias15_integrator.simulation(
                        self,
                        self.stats.objects_count,
                        self.m,
                        Grav_obj.G,
                        self.settings.tolerance,
                        self.settings.expected_time_scale,
                        self.settings.max_iteration,
                        self.settings.min_iteration,
                    )

        if self.is_c_lib == True:
            try:
                self.stats.total_energy = self.c_lib.compute_energy(
                    ctypes.c_int(self.stats.objects_count),
                    self.x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    self.v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    self.m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    ctypes.c_double(Grav_obj.G),
                )
            except:
                self.stats.total_energy = compute_energy(
                    self.stats.objects_count, self.x, self.v, self.m, Grav_obj.G
                )
        elif self.is_c_lib == False:
            self.stats.total_energy = compute_energy(
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
