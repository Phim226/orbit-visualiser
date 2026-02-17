from __future__ import annotations
from typing import Callable, TypeVar, Generic, TYPE_CHECKING
from dataclasses import dataclass
if TYPE_CHECKING:
    from orbit_visualiser.core import Orbit, Satellite, CentralBody

T = TypeVar("T")

@dataclass(frozen = True)
class PropertySpec(Generic[T]):
    label: str
    units: str | None
    getter: Callable[[T], float | bool]

@dataclass(frozen = True)
class VariableSpec(PropertySpec):
    init_value: float
    slider_lims: tuple[int]
    decimal_places: int
    entry_pos: tuple[int]