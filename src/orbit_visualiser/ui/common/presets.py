from dataclasses import dataclass

@dataclass
class Preset():
    eccentricity: float
    radius_of_periapsis: float
    true_anomaly: float
    right_ascension_of_the_ascending_node: float
    inclination: float
    argument_of_periapsis: float
    gravitational_parameter: float
    radius: float

initial_config: Preset = Preset(0.0, 50_000.0, 0.0, 0.0, 0.0, 0.0, 398600.0, 6378.0)