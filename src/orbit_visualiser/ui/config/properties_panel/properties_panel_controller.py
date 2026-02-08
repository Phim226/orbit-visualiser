import numpy as np
from orbit_visualiser.ui.config.properties_panel.properties_panel_builder import PropertiesBuilder

# TODO: Allow for temporary increase in slider scale when inputting manual values.
# TODO: Allow for fractional manual inputs.
# TODO: Remove any leading 0s from manual inputs.
# TODO: Refactor to a lazy/cached recalculation model. Currently everything is recalculated on every variable change.
# TODO: Properly display very large values in the properties panel without them being cut off.

class PropertiesController():


    def __init__(self, builder: PropertiesBuilder):
        self._builder = builder

    def update_display(self) -> None:
        for property, spec in list(self._builder.property_specs.items()):
            new_value = spec.getter(spec.obj)
            unit = spec.units
            getattr(self._builder, f"{property}_str").set(self.format_display_value(new_value, unit))

    def format_display_value(self, value: float | str, units: str | None) -> str:
        if units is None:
            return value.capitalize()

        if np.isclose(value, 0):
            value = 0.00

        if np.isneginf(value):
            return f"-∞ {units}"

        elif np.isinf(value):
            return f"∞ {units}"

        elif np.isnan(value):
            return "n/a"

        elif units in ["km", "km²/s", "s"]:
            return f"{value:6.0f} {units}"

        elif units in ["°", "km/s", "km²/s²"]:
            return f"{value:6.2f} {units}"

        elif units in ["°/s"]:
            return f"{value:6.6f} {units}"

    def draw_periapsis() -> None:
        pass