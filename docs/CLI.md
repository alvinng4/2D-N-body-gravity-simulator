# Command-line interface (CLI) documentation for Gravity-Simulator

* [Quick Start](#quick-start)
    - [Python version](#python-version)
    - [Installation](#installation)
    - [Some notes](#some-notes)
* [Running the program in terminal](#running-the-program-in-terminal)
* [Default systems](#default-systems)
* [Integrators](#integrators)
    - [Simple methods](#simple-methods)
    - [Embedded Runge-Kutta methods](#embdedded-runge-kutta-methods)
    - [IAS15](#IAS15)
    - [WHFast](#whfast)
* [Making changes to saved systems](#making-changes-to-saved-systems)
* [Saving the results](#saving-the-results)

## Quick Start

### Python version
This program requires Python version 3.10 or higher. 

### Installation
1. Download the source files, or clone this repository by running the following command in terminal:
    ```
    git clone https://github.com/alvinng4/Gravity-Simulator
    ```
2. Install the required packages by
    ```
    pip install .
    ```
    If the installation is not successful, install the following packages manually:
    ```
    matplotlib==3.8.3
    numpy==1.26.4
    rich==13.7.1
    ```
3. I have compiled the C library for you, but you may want to recompile the C library if 
    * The program failed to load the C library, and you don't want to use the NumPy option (slow)
    * You want to use CUDA GPU acceleration

    To compile the C library, simply go to the src folder and run
    ```
    make [CC=gcc] [USE_CUDA=1]
    ```
    Then, move the `c_lib.dylib`, `c_lib.dll` or `c_lib.so` file into the gravity_sim folder (One of them will be generated depending on your operation system).
    ```
    mv c_lib.dylib ../gravity_sim
    mv c_lib.dll   ../gravity_sim
    mv c_lib.so    ../gravity_sim
    ```
    If you wish to use CUDA GPU acceleration, you will also need to recompile
    the library with the `USE_CUDA=1` flag, which requires `nvcc` to be installed in your system.

### Some notes
* The default unit for this project is solar masses, AU and days, with $G = 0.00029591220828411956 \text{ M}_\odot^{-1} \text{ AU}^3 \text{ day}^{-2}$.
* Animations, simulation results, etc. will be stored to `gravity_sim/result` by default, unless a file path is specified.
* API would be more suitable for complex simulations

## Running the program in terminal

Once the installation is completed, navigate to the source directory in terminal and run
```
python gravity_sim [-n|--numpy]
```
`-n, --numpy`: run the program with NumPy instead of C library

Simply follows the instructions from the program. Enjoy!

## Default systems
Some systems are available by default.
| System | Description |
|:-------|:------------| 
| circular_binary_orbit | A circular orbit formed by two stars |
| eccentric_binary_orbit | An eccentric orbit formed by two stars |
| 3d_helix | An upward helix consists of three stars |
| sun_earth_moon | The Sun, Earth, and Moon system |
| figure-8 | A "figure-8" orbit involving three stars  |
| pyth-3-body | Three stars arranged in a triangle with length ratios of 3, 4, and 5. It is a highly chaotic orbit with close encounters that can be used to test the difference between fixed and variable step size integrators. |
| solar_system | Solar System with the Sun and the planets |
| solar_system_plus | solar_system with the inclusion of Pluto, Ceres, and Vesta  |

## Integrators 
### Simple methods
Below are four simple fixed step size methods to simulate the system with a given step size $\text{d}t$.
| Simple methods |
|:-----------|
| Euler |
| Euler Cromer |
| Fourth Order Runge-Kutta (RK4) |
| Leapfrog |

### Embedded Runge-Kutta methods
Embedded RK methods are adaptive methods that decides the step size automatically based on the estimated error.
It can resolve close encounters but fail to conserve energy over long time scele.

| Embdedded Runge-Kutta methods | Recommended tolerance* |
|:-----------|:-------------|
| Runge–Kutta–Fehlberg 4(5) | $10^{-8}$ to $10^{-14}$ |
| Dormand–Prince method (DOPRI) 5(4) | $10^{-8}$ to $10^{-14}$ |
| Verner's method (DVERK) 6(5) | $10^{-8}$ to $10^{-14}$ |
| Runge–Kutta–Fehlberg 7(8) | $10^{-4}$ to $10^{-8}$ |

*For reference only

### IAS15
IAS15 (Implicit integrator with Adaptive time Stepping, 15th order) is a highly optimized integrator with extremely high accuracy. It is the default method for this project.

The recommended tolerance* is $10^{-9}$. Since the integrator is 15th order, changing the tolerance
results in little improvement in performance, but a huge penalty in accuracy. Therefore, it is not
recommended to change this tolerance.

*For reference only

### WHFast
WHFast is a second order symplectic method with fixed step size, which conserves energy over long integration period. This integrator cannot resolve close encounter.

## Making changes to saved systems

If you wish to make any changes to saved systems, you can access the file at 
```
gravity_simulator/gravity_sim/customized_systems.csv
``` 
The data follow the format
```
Name, Gravitational constant, Number of objects, m1, ..., x1, y1, z1, ..., vx1, vy1, vz1, ...
```

## Saving the results
If you saved the results, the numerical data will be stored in the following folder:
```
Gravity-Simulator/gravity_sim/results
```
The file will starts with the metadata which starts with `#`.
Missing information will be saved as `None`.
More rows may be added in the future.

Below is an example:
```
# Data saved on (YYYY-MM-DD): 2024-07-27
# System Name: solar_system
# Integrator: IAS15
# Number of objects: 9
# Gravitational constant: 0.00029591220828411956
# Simulation time (days): 73048.4378
# dt (days): None
# Tolerance: 1e-09
# Data size: 64596
# Store every nth point: 1
# Run time (s): 1.4667500000214204
# masses: 1.0 1.6601208254589484e-07 2.447838287796944e-06 3.0034896154649684e-06 3.2271560829322774e-07 0.0009547919099414248 0.00028588567002459455 4.36624961322212e-05 5.151383772654274e-05
```
Then, the actual data will be saved in the default unit (solar masses, AU and days), and follow this format:
```
time, dt, total energy, x1, y1, z1, ... vx1, vy1, vz1, ...
```
The saved data file can be read by the program.
Even if the metadata is corrupted or missing, the program can still read the data file, although some information could be missing.
