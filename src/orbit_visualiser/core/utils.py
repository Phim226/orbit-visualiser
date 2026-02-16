import numpy as np
from numpy.typing import NDArray

def promote_to_3d(vector: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    Promotes a 2D numpy array to 3D, and returns vectors of other shapes unchanged.

    Parameters
    ----------
    vector : NDArray[np.float64]
        Numpy array

    Returns
    -------
    NDArray[np.float64]
        Promoted array
    """
    if vector.shape[-1] == 2:
        return np.append(vector, 0.0)

    return vector