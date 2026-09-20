from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from scripts.generate_package_map import public_symbols


ROOT = Path(__file__).parents[1]
MAP = ROOT / "docs" / "package-map.md"


def test_package_map_is_current():
    generated = subprocess.run(
        [sys.executable, "scripts/generate_package_map.py"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert generated.returncode == 0, generated.stdout + generated.stderr
    result = subprocess.run(
        [sys.executable, "scripts/generate_package_map.py", "--check"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_public_symbols_exclude_imports_and_private_names():
    source = """
from elsewhere import Imported
__all__ = ["Visible", "CONSTANT", "Imported"]
CONSTANT = 1
_private = 2
class Visible: pass
def hidden(): pass
"""
    symbols = public_symbols(source)
    assert [(symbol.name, symbol.line) for symbol in symbols] == [
        ("CONSTANT", 4),
        ("Visible", 6),
    ]


def test_package_map_contains_representative_modules_and_symbols():
    content = MAP.read_text(encoding="utf-8")
    assert "[`euskera/evolution/evolve.py`](../euskera/evolution/evolve.py)" in content
    assert "[`evolve`](../euskera/evolution/evolve.py#L" in content
    assert "[`Models`](../euskera/models/models.py#L" in content
    assert "[`__init__`](../euskera/evolution/__init__.py#L" not in content
