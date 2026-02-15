from math import floor
import numpy as np

def floor_to_decimals(num: float, decimals: int) -> float:
    """
    Floors a float to the given number of decimal points. For example floor_to_decimals(1.236, 2)
    returns 1.23.

    Parameters
    ----------
    num : float
        Float to be floored
    decimals : int
        Number of decimals to floor up to

    Returns
    -------
    float
        Floored float
    """
    factor = 10**decimals
    return floor(num*factor)/factor

def infinity(nu: float) -> np.float64:
    """
    Returns a signed infinity according to the sign of input.

    Parameters
    ----------
    nu : float
        Input float determining the sign of the infinity

    Returns
    -------
    np.float64
        The signed infinity
    """
    return np.sign(nu)*np.inf