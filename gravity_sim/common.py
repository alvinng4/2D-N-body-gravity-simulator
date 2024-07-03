import numpy as np


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


def compute_energy(objects_count, x, v, m, G):
    E = 0
    for j in range(objects_count):
        E += 0.5 * m[j] * np.linalg.norm(v[j]) ** 2
        for k in range(j + 1, objects_count):
            R = x[j] - x[k]
            norm = np.linalg.norm(R)
            if norm != 0:
                E -= G * m[j] * m[k] / norm
            else:
                return np.nan
    return E
