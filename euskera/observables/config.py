"""Typed configuration for conserved quantities and diagnostics."""

from dataclasses import asdict, dataclass, fields
from typing import Any, Mapping

__all__ = ["DiagnosticsConfig"]


@dataclass
class DiagnosticsConfig:
    Numb_Part: bool = True
    Energ: bool = True
    Pi: bool = False
    Ji: bool = False
    Frequency: bool = False

    def __post_init__(self) -> None:
        for name in ("Numb_Part", "Energ", "Pi", "Ji", "Frequency"):
            if not isinstance(getattr(self, name), bool):
                raise TypeError(f"{name} must be a boolean")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_legacy(cls, values: Mapping[str, Any] | None = None) -> "DiagnosticsConfig":
        values = {} if values is None else dict(values)
        allowed = {field.name for field in fields(cls)}
        unknown = set(values) - allowed
        if unknown:
            raise ValueError(f"Unknown diagnostics configuration keys: {sorted(unknown)}")
        return cls(**values)
