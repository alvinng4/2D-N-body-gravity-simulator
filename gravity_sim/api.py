"""
Application user interface (API) for gravity simulator

Usage:
    ffrom gravity_sim import GravitySimulatorAPI
    grav_sim = GravitySimulatorAPI()

Author:  Ching Yin Ng
"""

import copy
import ctypes
import sys
import warnings
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

        self.BUILT_IN_SYSTEMS = (
            self.simulator.RECOMMENDED_SETTINGS_BUILT_IN_SYSTEMS.keys()
        )
        self.AVAILABLE_INTEGRATORS = self.simulator.AVAILABLE_INTEGRATORS
        self.AVAILABLE_ACCELERATION_METHODS = (
            self.simulator.AVAILABLE_ACCELERATION_METHODS
        )
        self.AVAILABLE_STORING_METHODS = self.simulator.AVAILABLE_STORING_METHODS

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
        tf: float,
        gravitational_system: GravitationalSystem,
        integrator_params: dict,
        acceleration_params: Optional[dict] = None,
        storing_params: Optional[dict] = None,
        settings: Optional[dict] = None,
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
        if acceleration_params is None:
            acceleration_params = {
                "method": "pairwise",
                "softening_length": 0.0,
            }
        if storing_params is None:
            storing_params = {
                "method": "default",
                "storing_freq": 1,
            }
        if settings is None:
            settings = {
                "verbose": 2,
                "disable_progress_bar": False,
                "make_copy_params": True,
                "make_copy_system": True,
            }
        if "make_copy_params" in settings:
            if settings["make_copy_params"]:
                integrator_params = integrator_params.copy()
                acceleration_params = acceleration_params.copy()
                storing_params = storing_params.copy()
                settings = settings.copy()
            if "make_copy_system" in settings:
                if settings["make_copy_system"]:
                    gravitational_system = copy.deepcopy(gravitational_system)

        self._check_simulation_input(
            gravitational_system,
            integrator_params,
            acceleration_params,
            storing_params,
            settings,
            tf,
        )
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

    def _check_simulation_input(
        self,
        gravitational_system: GravitationalSystem,
        integrator_params: dict,
        acceleration_params: dict,
        storing_params: dict,
        settings: dict,
        tf: float,
    ) -> None:
        """Check simulation input

        Parameters
        ----------
        gravitational_system : GravitationalSystem
        integrator_params : dict
        acceleration_params : dict
        storing_params : dict
        settings : dict
        tf : float

        Raises
        ------
        TypeError
            If input types are incorrect
        ValueError
            If input values are invalid
        """
        ### GravitationalSystem
        if not isinstance(gravitational_system, GravitationalSystem):
            raise TypeError(
                f"Expected GravitationalSystem object, but got {type(gravitational_system)}"
            )
        if gravitational_system.objects_count == 0:
            raise ValueError("No objects in the gravitational system")

        ### integrator_params ###
        if not isinstance(integrator_params, dict):
            raise TypeError(f"Expected dict, but got {type(integrator_params)}")
        if "integrator" not in integrator_params:
            raise ValueError('integrator_params must have key "integrator"')
        if integrator_params["integrator"] not in self.simulator.AVAILABLE_INTEGRATORS:
            raise ValueError(
                f'integrator_params["integrator"] must be one of {self.simulator.AVAILABLE_INTEGRATORS}'
            )

        if (
            integrator_params["integrator"]
            in self.simulator.FIXED_STEP_SIZE_INTEGRATORS
        ):
            if "dt" not in integrator_params:
                raise ValueError(
                    'integrator_params must have key "dt" for fixed-step-size integrators'
                )
            if not isinstance(integrator_params["dt"], (int, float)):
                raise TypeError(
                    f"Expected int or float, but got {type(integrator_params['dt'])}"
                )
            if integrator_params["dt"] <= 0.0:
                raise ValueError('integrator_params["dt"] must be positive')
            if "tolerance" in integrator_params:
                warnings.warn(
                    'integrator_params["tolerance"] is not used for fixed-step-size integrators'
                )
            if "initial_dt" in integrator_params:
                warnings.warn(
                    'integrator_params["initial_dt"] is not used for fixed-step-size integrators'
                )
        elif (
            integrator_params["integrator"]
            in self.simulator.ADAPTIVE_STEP_SIZE_INTEGRATORS
        ):
            if "tolerance" not in integrator_params:
                raise ValueError(
                    'integrator_params must have key "tolerance" for adaptive-step-size integrators'
                )
            if not isinstance(integrator_params["tolerance"], (int, float)):
                raise TypeError(
                    f"Expected int or float, but got {type(integrator_params['tolerance'])}"
                )
            if integrator_params["tolerance"] <= 0.0:
                raise ValueError('integrator_params["tolerance"] must be positive')
            if "initial_dt" in integrator_params:
                if not isinstance(integrator_params["initial_dt"], (int, float)):
                    raise TypeError(
                        f"Expected int or float, but got {type(integrator_params['initial_dt'])}"
                    )
                if integrator_params["initial_dt"] <= 0.0:
                    raise ValueError('integrator_params["initial_dt"] must be positive')
            if "dt" in integrator_params:
                warnings.warn(
                    'integrator_params["dt"] is not used for adaptive-step-size integrators'
                )

        if integrator_params["integrator"] != "whfast":
            if "whfast_kepler_tol" in integrator_params:
                warnings.warn(
                    'integrator_params["whfast_kepler_tol"] is only used for WHFast integrator'
                )
            if "whfast_kepler_max_iter" in integrator_params:
                warnings.warn(
                    'integrator_params["whfast_kepler_max_iter"] is only used for WHFast integrator'
                )
            if "whfast_kepler_auto_remove" in integrator_params:
                warnings.warn(
                    'integrator_params["whfast_kepler_auto_remove"] is only used for WHFast integrator'
                )
            if "whfast_auto_remove_tol" in integrator_params:
                warnings.warn(
                    'integrator_params["whfast_auto_remove_tol"] is only used for WHFast integrator'
                )
        else:
            if "whfast_kepler_tol" in integrator_params:
                if not isinstance(integrator_params["whfast_kepler_tol"], (int, float)):
                    raise TypeError(
                        f"Expected int or float, but got {type(integrator_params['whfast_kepler_tol'])}"
                    )
                else:
                    if integrator_params["whfast_kepler_tol"] <= 0.0:
                        raise ValueError(
                            'integrator_params["whfast_kepler_tol"] must be positive'
                        )
            if "whfast_kepler_max_iter" in integrator_params:
                if not isinstance(integrator_params["whfast_kepler_max_iter"], int):
                    raise TypeError(
                        f"Expected int, but got {type(integrator_params['whfast_kepler_max_iter'])}"
                    )
                else:
                    if integrator_params["whfast_kepler_max_iter"] <= 0:
                        raise ValueError(
                            'integrator_params["whfast_kepler_max_iter"] must be positive'
                        )
            if "whfast_kepler_auto_remove" in integrator_params:
                if not isinstance(integrator_params["whfast_kepler_auto_remove"], bool):
                    raise TypeError(
                        f"Expected bool, but got {type(integrator_params['whfast_kepler_auto_remove'])}"
                    )
            if "whfast_auto_remove_tol" in integrator_params:
                if not isinstance(
                    integrator_params["whfast_auto_remove_tol"], (int, float)
                ):
                    raise TypeError(
                        f"Expected int or float, but got {type(integrator_params['whfast_auto_remove_tol'])}"
                    )
                else:
                    if integrator_params["whfast_auto_remove_tol"] <= 0.0:
                        raise ValueError(
                            'integrator_params["whfast_auto_remove_tol"] must be positive'
                        )

        for key in ["dt", "tolerance", "initial_dt"]:
            if key not in integrator_params:
                integrator_params[key] = 0.0

        if "whfast_kepler_tol" not in integrator_params:
            integrator_params["whfast_kepler_tol"] = 1e-12
        if "whfast_kepler_max_iter" not in integrator_params:
            integrator_params["whfast_kepler_max_iter"] = 500
        if "whfast_kepler_auto_remove" not in integrator_params:
            integrator_params["whfast_kepler_auto_remove"] = False
        if "whfast_auto_remove_tol" not in integrator_params:
            integrator_params["whfast_auto_remove_tol"] = 1e-8

        ### acceleration_params ###
        if not isinstance(acceleration_params, dict):
            raise TypeError(f"Expected dict, but got {type(acceleration_params)}")
        if "method" not in acceleration_params:
            raise ValueError('acceleration_params must have key "method"')
        if (
            acceleration_params["method"]
            not in self.simulator.AVAILABLE_ACCELERATION_METHODS
        ):
            raise ValueError(
                f'acceleration_params["method"] must be one of {self.simulator.AVAILABLE_ACCELERATION_METHODS}'
            )
        if "softening_length" in acceleration_params:
            if not isinstance(acceleration_params["softening_length"], (int, float)):
                raise TypeError(
                    f"Expected int or float, but got {type(acceleration_params['softening_length'])}"
                )
            if acceleration_params["softening_length"] < 0.0:
                raise ValueError(
                    'acceleration_params["softening_length"] must be non-negative'
                )
        else:
            acceleration_params["softening_length"] = 0.0

        if "order" in acceleration_params:
            if not isinstance(acceleration_params["order"], int):
                raise TypeError(
                    f"Expected int, but got {type(acceleration_params['order'])}"
                )
            if acceleration_params["order"] < 0:
                raise ValueError('acceleration_params["order"] must be non-negative')
        else:
            acceleration_params["order"] = 0

        if "opening_angle" in acceleration_params:
            if not isinstance(acceleration_params["opening_angle"], (int, float)):
                raise TypeError(
                    f"Expected int or float, but got {type(acceleration_params['opening_angle'])}"
                )
            if acceleration_params["opening_angle"] <= 0.0:
                raise ValueError(
                    'acceleration_params["opening_angle"] must be positive'
                )
        else:
            acceleration_params["opening_angle"] = 0.5

        ### storing_params ###
        if not isinstance(storing_params, dict):
            raise TypeError(f"Expected dict, but got {type(storing_params)}")
        if "method" not in storing_params:
            raise ValueError('storing_params must have key "method"')
        if storing_params["method"] not in self.simulator.AVAILABLE_STORING_METHODS:
            raise ValueError(
                f'storing_params["method"] must be one of {self.simulator.AVAILABLE_STORING_METHODS}'
            )
        if storing_params["method"] in ["default", "disabled"]:
            if "flush_path" in storing_params:
                warnings.warn(
                    'storing_params["flush_path"] is not used for default storing method'
                )
        if storing_params["method"] == "flush":
            if "flush_path" not in storing_params:
                raise ValueError(
                    'storing_params must have key "flush_path" for flush storing method'
                )
            if not (
                isinstance(storing_params["flush_path"], str)
                or isinstance(storing_params["flush_path"], Path)
            ):
                raise TypeError(
                    f"Expected str or Path, but got {type(storing_params['flush_path'])}"
                )

        if storing_params["method"] != "disabled":
            if "storing_freq" not in storing_params:
                raise ValueError('storing_params must have key "storing_freq"')
            else:
                if not isinstance(storing_params["storing_freq"], int):
                    raise TypeError(
                        f"Expected int, but got {type(storing_params['storing_freq'])}"
                    )
                if storing_params["storing_freq"] <= 0:
                    raise ValueError('storing_params["storing_freq"] must be positive')
        else:
            if "storing_freq" in storing_params:
                warnings.warn(
                    'storing_params["storing_freq"] is not used for disabled storing method'
                )

        ### settings ###
        if "verbose" not in settings:
            settings["verbose"] = 2
        else:
            if not isinstance(settings["verbose"], int):
                raise TypeError(f"Expected int, but got {type(settings['verbose'])}")
        if "disable_progress_bar" not in settings:
            settings["disable_progress_bar"] = False
        else:
            if not isinstance(settings["disable_progress_bar"], bool):
                raise TypeError(
                    f"Expected bool, but got {type(settings['disable_progress_bar'])}"
                )

        settings["make_copy_params"] = False
        settings["make_copy_system"] = False

        ### tf ###
        if not isinstance(tf, (int, float)):
            raise TypeError(f"Expected int or float, but got {type(tf)}")

        if tf <= 0.0:
            raise ValueError("tf must be positive")

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
