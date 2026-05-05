from enum import Enum

class OrbitType(Enum):
    """
    Enum representing the conic type of an analytical Keplerian orbit.
    """
    CIRCULAR = 1
    ELLIPTICAL = 2
    PARABOLIC = 3
    HYPERBOLIC = 4

class OrbitMotion(Enum):
    """
    Enum representing whether an analytical Keplerian orbit is prograde or retrograde.
    """
    PROGRADE = 1
    RETROGRADE = 2