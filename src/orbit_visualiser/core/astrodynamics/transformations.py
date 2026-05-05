import numpy as np
from numpy.typing import NDArray

def raan_rot_mat(raan: float) -> NDArray[np.float64]:
    c_raan = np.cos(raan)
    s_raan = np.sin(raan)
    return np.array([[c_raan,  s_raan,  0],
                     [-s_raan, c_raan,  0],
                     [0,       0,       1]])

def inclination_rot_mat(i: float) -> NDArray[np.float64]:
    c_i = np.cos(i)
    s_i = np.sin(i)
    return np.array([[1,    0,   0],
                     [0,  c_i, s_i],
                     [0, -s_i, c_i]])

def arg_of_periapsis_rot_mat(omega: float) -> NDArray[np.float64]:
    c_omega = np.cos(omega)
    s_omega = np.sin(omega)
    return np.array([[c_omega,  s_omega,  0],
                     [-s_omega, c_omega,  0],
                     [0,        0,        1]])

def eci_to_perifocal_trans_mat(raan: float, i: float, omega: float) -> NDArray[np.float64]:
    """
    Returns the transformation matrix from the Earth Centred Inertial frame to the Perifocal frame.

    Parameters
    ----------
    raan : float
        Right ascension of the ascending node (rad)
    i : float
        Orbital inclination (rad)
    omega : float
        Argument of periapsis (rad)

    Returns
    -------
    NDArray[np.float64]
        ECI to Perifocal transformation matrix
    """
    return np.matmul(np.matmul(arg_of_periapsis_rot_mat(omega), inclination_rot_mat(i)), raan_rot_mat(raan))

def perifocal_to_eci_trans_mat(raan: float, i: float, omega: float) -> NDArray[np.float64]:
    """
    Returns the transformation matrix from the Perifocal frame to the Earth Centred Inertial frame.

    Parameters
    ----------
    raan : float
        Right ascension of the ascending node (rad)
    i : float
        Orbital inclination (rad)
    omega : float
        Argument of periapsis (rad)

    Returns
    -------
    NDArray[np.float64]
        Perifocal to ECI transformation matrix
    """
    return np.transpose(eci_to_perifocal_trans_mat(raan, i, omega))