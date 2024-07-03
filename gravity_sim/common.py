import numpy as np


def acceleration(objects_count, x, m, G):
    """
    Calculate acceleration by a = - GM/r^3 vec{r} using a matrix approach

    Old version of code using nested for loops:
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
    """
    # Compute the displacement vector
    r_ij = x[np.newaxis, :, :] - x[:, np.newaxis, :]

    # Compute the distance squared
    r_squared = np.sum(r_ij ** 2, axis=2)
    softening = 1e-12   # To prevent r = 0
    r_squared[r_squared < softening] = softening

    # Compute 1/r^3
    inv_r_cubed = r_squared ** (-1.5)

    # Set diagonal elements to 0 to avoid self-interaction
    np.fill_diagonal(inv_r_cubed, 0.0)

    # Compute the acceleration
    a = G * np.sum(r_ij * inv_r_cubed[:,:,np.newaxis] * m[np.newaxis, :, np.newaxis], axis=1)

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
