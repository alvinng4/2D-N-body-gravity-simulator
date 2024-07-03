import re

import numpy as np


def get_bool(msg: str) -> bool:
    """
    Get boolean input from user

    Args: msg (str): Message to display to user
    Print: "{msg} (y/n): "
    Returns: bool
    """
    while True:
        if matches := re.search(
            r"^\s*(yes|no|y|n)\s*$", input(f"{msg} (y/n): "), re.IGNORECASE
        ):
            if matches.group(1).lower() in ["y", "yes"]:
                return True
            elif matches.group(1).lower() in ["n", "no"]:
                return False

        print("Invalid input. Please try again.\n")


def get_int(
    msg: str,
    larger_than: int = None,
    smaller_than: int = None,
    allow_cancel: bool = False,
) -> int:
    """
    Get integer input from user

    Args: msg (str): Message to display to user
          larger_than (int): Value must be larger than this
          smaller_than (int): Value must be smaller than this
          allow_cancel (bool): Allow user to enter "cancel" to cancel input
    Print: "{msg}"
    Returns: int, None if allow_cancel is True and user enters "cancel"
    """
    while True:
        try:
            value = input(f"{msg}")
            if allow_cancel and value.strip().lower() == "cancel":
                return None
            else:
                value = int(value)

            if (larger_than is None or value > larger_than) and (
                smaller_than is None or value < smaller_than
            ):
                return value

            if larger_than is not None and value <= larger_than:
                print("Value too small! Please try again.")
                print()

            elif smaller_than is not None and value >= smaller_than:
                print("Value too big! Please try again.")
                print()

        except ValueError:
            print("Invalid input. Please try again.")
            print()


def get_float(
    msg: str,
    larger_than: float = None,
    smaller_than: float = None,
    allow_cancel: bool = False,
) -> float:
    """
    Get integer input from user

    Args: msg (str): Message to display to user
          larger_than (float): Value must be larger than this
          smaller_than (float): Value must be smaller than this
          allow_cancel (bool): Allow user to enter "cancel" to cancel input
    Print: "{msg}"
    Returns: float, None if allow_cancel is True and user enters "cancel"
    """
    while True:
        try:
            value = input(f"{msg}")
            if allow_cancel and value.strip().lower() == "cancel":
                return None
            else:
                value = float(value)

            if (larger_than is None or value > larger_than) and (
                smaller_than is None or value < smaller_than
            ):
                return value

            if larger_than is not None and value <= larger_than:
                print("Value too small! Please try again.")
                print()

            elif smaller_than is not None and value >= smaller_than:
                print("Value too big! Please try again.")
                print()

        except ValueError:
            print("Invalid input. Please try again.")
            print()


def acceleration(objects_count, x, m, G):
    """
    Calculate acceleration by a = - GM/r^3 vec{r}
    """
    # Allocate memory
    a = np.zeros((objects_count, 3))

    # Calculations
    for j in range(objects_count):
        for k in range(j + 1, objects_count):
            R = x[j] - x[k]
            temp_value = G * R / np.linalg.norm(R) ** 3
            a[j] += -temp_value * m[k]
            a[k] += temp_value * m[j]

    return a
