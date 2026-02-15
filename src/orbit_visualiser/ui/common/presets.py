from dataclasses import dataclass

@dataclass
class Preset():
    eccentricity: float
    radius_of_periapsis: float
    gravitational_parameter: float
    true_anomaly: float
    radius: float

initial_config: Preset = Preset(0.0, 50_000.0, 398600.0, 0.0, 6378.0)