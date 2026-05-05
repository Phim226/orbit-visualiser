from __future__ import annotations
from typing import Callable, TypeVar, Generic, Literal
from dataclasses import dataclass

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
    init_state: Literal["normal", "disabled"]