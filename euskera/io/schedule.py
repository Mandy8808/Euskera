"""Validated temporal selections and per-output sampling rules."""

from dataclasses import dataclass, field
import math

__all__ = ["All", "Last", "TimeRange", "Final", "SaveRule"]


@dataclass(frozen=True)
class All:
    """Select every eligible sample, including zero."""
    kind: str = field(default="all", init=False)


@dataclass(frozen=True)
class Last:
    """Select the last count eligible samples, including the final sample."""
    count: int
    kind: str = field(default="last", init=False)

    def __post_init__(self):
        _positive_integer(self.count, "count")


@dataclass(frozen=True)
class TimeRange:
    """Select samples in an inclusive physical-time interval."""
    start: float
    stop: float
    kind: str = field(default="range", init=False)

    def __post_init__(self):
        for value in (self.start, self.stop):
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError("TimeRange bounds must be finite numbers")
        if not 0 <= self.start <= self.stop:
            raise ValueError("TimeRange requires 0 <= start <= stop")


@dataclass(frozen=True)
class Final:
    """Select only the final sample."""
    kind: str = field(default="final", init=False)


def _positive_integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")


@dataclass
class SaveRule:
    """One geometry, variables and temporal selection.

    Diagnostics use variables from DiagnosticsConfig and no orientations.
    ``every`` subsamples the common calendar before temporal selection; the
    final sample is always eligible. ``save_initial`` adds sample zero.
    """
    geometry: str
    variables: list[str]
    orientations: list[str] = field(default_factory=list)
    selection: object = field(default_factory=All)
    every: int = 1
    save_initial: bool = False

    def __post_init__(self):
        if self.geometry not in ("volume", "plane", "line", "diagnostics"):
            raise ValueError("Unknown output geometry")
        allowed = ("Numb_Part", "Energ", "Pi", "Frequency") if self.geometry == "diagnostics" else ("rho", "psi", "phi")
        if not isinstance(self.variables, (list, tuple)) or not self.variables or any(v not in allowed for v in self.variables):
            raise ValueError(f"variables must be a nonempty sequence drawn from {allowed}")
        axes = {"plane": ("xy", "xz", "yz"), "line": ("x", "y", "z")}.get(self.geometry, ())
        if not isinstance(self.orientations, (list, tuple)) or any(o not in axes for o in self.orientations):
            raise ValueError("Invalid orientations for geometry")
        if axes and not self.orientations:
            raise ValueError("Planes and lines require explicit orientations")
        if len(set(self.variables)) != len(self.variables) or len(set(self.orientations)) != len(self.orientations):
            raise ValueError("Duplicate variables or orientations")
        if isinstance(self.selection, dict):
            values = dict(self.selection)
            kind = values.pop("kind", None)
            types = {"all": All, "last": Last, "range": TimeRange, "final": Final}
            if kind not in types:
                raise ValueError("Unknown temporal selection")
            self.selection = types[kind](**values)
        if not isinstance(self.selection, (All, Last, TimeRange, Final)):
            raise TypeError("selection must be All, Last, TimeRange or Final")
        _positive_integer(self.every, "every")
        if not isinstance(self.save_initial, bool):
            raise TypeError("save_initial must be boolean")

    def matches(self, index, total, time, duration):
        """Evaluate a selection without allocating a full sampling calendar."""
        if index == 0 and self.save_initial:
            return True
        if index % self.every and index != total:
            return False
        selection = self.selection
        if isinstance(selection, All):
            return True
        if isinstance(selection, Final):
            return index == total
        if isinstance(selection, Last):
            size = total // self.every + 1 + bool(total % self.every)
            rank = size - 1 if index == total else index // self.every
            return rank >= max(0, size - selection.count)
        tolerance = 8 * math.ulp(max(duration, 1.0))
        return selection.start - tolerance <= time <= selection.stop + tolerance
