"""Typed configuration for time evolution.

The dataclasses deliberately serialize to the dictionaries consumed by the
legacy numerical kernels.  This keeps the kernels stable while providing a
validated, discoverable public API.
"""

from dataclasses import asdict, dataclass, fields
from typing import Any, Mapping

__all__ = ["EvolutionConfig"]


@dataclass
class EvolutionConfig:
    lambda_value: float = 0
    num_threads: int = 1
    gridlength: float = 10
    resol: int = 128
    step_factor: float = 1.0
    t0: float = 0
    tmax: float = 1
    rmax: float = 7.6
    Plim: float = 5.6
    cmass: float = 0
    plott0: bool = False
    Boverlap: bool = True
    methodEnerg: int = 1
    info: bool = False

    def __post_init__(self) -> None:
        if isinstance(self.resol, bool) or not isinstance(self.resol, int) or self.resol <= 0:
            raise ValueError("resol must be a positive integer")
        if isinstance(self.num_threads, bool) or not isinstance(self.num_threads, int) or self.num_threads <= 0:
            raise ValueError("num_threads must be a positive integer")
        for name in ("gridlength", "step_factor", "tmax"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or value <= 0:
                raise ValueError(f"{name} must be a positive number")
        for name in ("lambda_value", "t0", "rmax", "Plim", "cmass"):
            if isinstance(getattr(self, name), bool) or not isinstance(getattr(self, name), (int, float)):
                raise TypeError(f"{name} must be numeric")
        if isinstance(self.methodEnerg, bool) or not isinstance(self.methodEnerg, int):
            raise TypeError("methodEnerg must be an integer")
        for name in ("plott0", "Boverlap", "info"):
            if not isinstance(getattr(self, name), bool):
                raise TypeError(f"{name} must be a boolean")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_legacy(cls, values: Mapping[str, Any] | None = None) -> "EvolutionConfig":
        values = {} if values is None else dict(values)
        allowed = {field.name for field in fields(cls)}
        unknown = set(values) - allowed
        if unknown:
            raise ValueError(f"Unknown evolution configuration keys: {sorted(unknown)}")
        return cls(**values)
