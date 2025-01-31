"""
Application user interface (API) for gravity simulator

Usage:
    from gravity_sim import GravitySimulator
    grav_sim = GravitySimulator()

Author:  Ching Yin Ng
"""

import ctypes
import sys
from pathlib import Path
from typing import Optional, Tuple

sys.path.append(str(Path(__file__).parent))

import numpy as np

from . import plotting
from . import utils
from .gravitational_system import GravitationalSystem
from .simulator import Simulator


class GravitySimulatorAPI:
    """Gravity simulator API"""

    def __init__(self, c_lib_path: Optional[str] = None) -> None:
        """
        Initialize gravity simulator API

        Parameters
        ----------
        c_lib_path : str, optional
            Path to C library, by default None
        """
        if c_lib_path is not None:
            self.c_lib = utils.load_c_lib(Path(c_lib_path))
        else:
            self.c_lib = utils.load_c_lib()
        utils.initialize_c_lib(self.c_lib)

        self.simulator = Simulator(self.c_lib)

        self.plot_2d_trajectory = plotting.plot_2d_trajectory
        self.plot_3d_trajectory = plotting.plot_3d_trajectory
        self.animate_2d_traj_gif = plotting.animate_2d_traj_gif
        self.animate_3d_traj_gif = plotting.animate_3d_traj_gif
        self.plot_quantity_against_time = plotting.plot_quantity_against_time
        self.plot_eccentricity_or_inclination = (
            plotting.plot_eccentricity_or_inclination
        )

    def days_to_years(self, days: float | np.ndarray) -> float | np.ndarray:
        return days / self.simulator.DAYS_PER_YEAR

    def years_to_days(self, years: float | np.ndarray) -> float | np.ndarray:
        return years * self.simulator.DAYS_PER_YEAR

    def create_system(self, name: Optional[str] = None) -> GravitationalSystem:
        """Create a gravitational system

        Parameters
        ----------
        name : str, optional
            Name of the system, by default None

        Returns
        -------
        GravitationalSystem object
        """
        return GravitationalSystem(name)

    def launch_simulation(
        self,
        gravitational_system: GravitationalSystem,
        integrator_params: dict,
        acceleration_params: dict,
        storing_params: dict,
        settings: dict,
        tf: float,
    ) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Launch simulation

        Parameters
        ----------
        gravitational_system : GravitationalSystem
        integrator_params : dict
        acceleration_params : dict
        storing_params : dict
        settings : dict
        tf : float
        """
        is_exit_ctypes_bool = ctypes.c_bool(False)
        try:
            self.simulator.launch_simulation(
                gravitational_system,
                integrator_params,
                acceleration_params,
                storing_params,
                settings,
                tf,
                is_exit_ctypes_bool,
            )
        except KeyboardInterrupt:
            is_exit_ctypes_bool.value = True
            raise KeyboardInterrupt

        if storing_params["method"] == "default":
            return (
                self.simulator.sol_state_,
                self.simulator.sol_time_,
                self.simulator.sol_dt_,
            )

        return None

    def compute_energy(
        self,
        gravitational_system: GravitationalSystem,
        sol_state: np.ndarray,
    ) -> np.ndarray:
        """Compute energy of the system

        Parameters
        ----------
        gravitational_system : GravitationalSystem
        sol_state : np.ndarray

        Returns
        -------
        sol_energy : np.ndarray
            Energy of the system at each time step
        """
        objects_count = gravitational_system.objects_count
        m = gravitational_system.m
        G = gravitational_system.G

        is_exit_ctypes_bool = ctypes.c_bool(False)
        try:
            return self.simulator.compute_energy(
                objects_count, m, G, sol_state, is_exit_ctypes_bool
            )
        except KeyboardInterrupt:
            is_exit_ctypes_bool.value = True
            raise KeyboardInterrupt

    def compute_linear_momentum(
        self,
        gravitational_system: GravitationalSystem,
        sol_state: np.ndarray,
    ) -> np.ndarray:
        """Compute linear momentum of the system

        Parameters
        ----------
        gravitational_system : GravitationalSystem
        sol_state : np.ndarray

        Returns
        -------
        sol_linear_momentum : np.ndarray
            Linear momentum of the system at each time step
        """
        objects_count = gravitational_system.objects_count
        m = gravitational_system.m

        is_exit_ctypes_bool = ctypes.c_bool(False)
        try:
            return self.simulator.compute_linear_momentum(
                objects_count, m, sol_state, is_exit_ctypes_bool
            )
        except KeyboardInterrupt:
            is_exit_ctypes_bool.value = True
            raise KeyboardInterrupt

    def compute_angular_momentum(
        self,
        gravitational_system: GravitationalSystem,
        sol_state: np.ndarray,
    ) -> np.ndarray:
        """Compute angular momentum of the system

        Parameters
        ----------
        gravitational_system : GravitationalSystem
        sol_state : np.ndarray

        Returns
        -------
        sol_angular_momentum : np.ndarray
            Angular momentum of the system at each time step
        """
        objects_count = gravitational_system.objects_count
        m = gravitational_system.m

        is_exit_ctypes_bool = ctypes.c_bool(False)
        try:
            return self.simulator.compute_angular_momentum(
                objects_count, m, sol_state, is_exit_ctypes_bool
            )
        except KeyboardInterrupt:
            is_exit_ctypes_bool.value = True
            raise KeyboardInterrupt

    def compute_eccentricity(
        self,
        gravitational_system: GravitationalSystem,
        sol_state: np.ndarray,
    ) -> np.ndarray:
        """Compute the eccentricity using the sol_state array,
        assuming that the first object is the central object

        Parameters
        ----------
        gravitational_system : GravitationalSystem
        sol_state : np.ndarray

        Returns
        -------
        sol_eccentricity : np.ndarray
            Eccentricity of the system at each time step
        """
        objects_count = gravitational_system.objects_count
        m = gravitational_system.m
        G = gravitational_system.G

        return self.simulator.compute_eccentricity(objects_count, m, G, sol_state)

    def compute_inclination(
        self,
        gravitational_system: GravitationalSystem,
        sol_state: np.ndarray,
    ) -> np.ndarray:
        """Compute the inclination using the sol_state array,
        assuming that the first object is the central object

        Parameters
        ----------
        gravitational_system : GravitationalSystem
        sol_state : np.ndarray

        Returns
        -------
        sol_inclination : np.ndarray
            Inclination of the system at each time step
        """
        objects_count = gravitational_system.objects_count

        return self.simulator.compute_inclination(objects_count, sol_state)
