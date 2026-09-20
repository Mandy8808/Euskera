"""Typed configuration for simulation output."""

from dataclasses import asdict, dataclass, field, fields
from typing import Any, Mapping

__all__ = ["OutputConfig"]


_DEFAULT_DATA_SAVE = {
    "grid": True, "save_rho": False, "save_psi": False, "save_phi": False,
    "save_plane_rho": False, "save_plane_psi": False, "save_plane_phi": False,
    "save_line_rho": True, "save_line_psi": True, "save_line_phi": True,
    "save_energies": True,
}


@dataclass
class OutputConfig:
    format: str = "npz"
    address: str = "Data"
    save_number: int = 10
    data_save: dict[str, bool] = field(default_factory=lambda: dict(_DEFAULT_DATA_SAVE))

    def __post_init__(self) -> None:
        if not isinstance(self.format, str) or not self.format:
            raise ValueError("format must be a non-empty string")
        if not isinstance(self.address, str) or not self.address:
            raise ValueError("address must be a non-empty string")
        if isinstance(self.save_number, bool) or not isinstance(self.save_number, int) or self.save_number <= 0:
            raise ValueError("save_number must be a positive integer")
        if not isinstance(self.data_save, Mapping):
            raise TypeError("data_save must be a mapping")
        if any(not isinstance(key, str) or not isinstance(value, bool) for key, value in self.data_save.items()):
            raise TypeError("data_save keys must be strings and values booleans")
        merged_data_save = dict(_DEFAULT_DATA_SAVE)
        merged_data_save.update(self.data_save)
        self.data_save = merged_data_save

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_legacy(cls, values: Mapping[str, Any] | None = None) -> "OutputConfig":
        values = {} if values is None else dict(values)
        allowed = {field.name for field in fields(cls)}
        unknown = set(values) - allowed
        if unknown:
            raise ValueError(f"Unknown output configuration keys: {sorted(unknown)}")
        return cls(**values)
