from enum import Enum

class OrbitType(Enum):
    """
    Enum representing the conic type of an analytical Keplerian orbit.
    """
    CIRCULAR = 1
    ELLIPTICAL = 2
    PARABOLIC = 3
    HYPERBOLIC = 4