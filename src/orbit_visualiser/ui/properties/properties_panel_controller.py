import numpy as np
from orbit_visualiser.ui.properties.properties_panel_builder import PropertiesBuilder
from orbit_visualiser.ui.data_access import OrbitDataAccess
from orbit_visualiser.core import OrbitType

# TODO: Properly display very large values in the properties panel without them being cut off.

class PropertiesController():


    def __init__(self, builder: PropertiesBuilder, oda: OrbitDataAccess):
        self._builder = builder
        self._oda = oda

    def update_display(self) -> None:
        for property, spec in list(self._builder.property_specs.items()):
            new_value = spec.getter(self._oda.satellite)
            unit = spec.units
            getattr(self._builder, f"{property}_str").set(self.format_display_value(new_value, unit))

    @staticmethod
    def format_display_value(value: float | OrbitType, unit: str | None) -> str:
        formatted_value: str = "Value not formatted"

        if unit is None:
            formatted_value = value.name.lower().capitalize()
            return formatted_value

        elif np.isneginf(value):
            formatted_value = f"-∞ {unit}"
            return formatted_value

        elif np.isinf(value):
            formatted_value = f"∞ {unit}"
            return formatted_value

        elif np.isnan(value):
            formatted_value = "n/a"
            return formatted_value


        if np.isclose(value, 0):
            value = 0.00

        format_mapping: dict[str, list[str]] = {
            f"{value:6.0f} {unit}": ["km", "km²/s", "s"],
            f"{value:6.2f} {unit}": ["°", "km/s", "km²/s²"],
            f"{value:6.6f} {unit}": ["°/s"]
        }

        for formatting, units in format_mapping.items():
            if unit in units:
                formatted_value = formatting
                break

        return formatted_value

    def draw_periapsis() -> None:
        pass