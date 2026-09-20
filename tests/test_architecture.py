"""Architecture and compatibility checks for the canonical package layout."""

import ast
import importlib
from pathlib import Path

import euskera.backgrounds as canonical


def test_background_canonical_modules_export_public_symbols():
    for module_name in (
        "systems",
        "solvers",
        "multifrequency",
        "asymptotics",
        "profiles",
        "fitting",
    ):
        module = importlib.import_module(f"euskera.backgrounds.{module_name}")
        assert module.__all__
        assert all(hasattr(module, name) for name in module.__all__)


def test_background_import_graph_has_no_cycles():
    root = Path(__file__).parents[1] / "euskera" / "backgrounds"
    graph = {}
    for path in root.glob("*.py"):
        module = f"euskera.backgrounds.{path.stem}"
        tree = ast.parse(path.read_text())
        deps = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.level:
                target = node.module or ""
                if target in {"systems", "solvers", "multifrequency", "asymptotics", "profiles", "fitting"}:
                    deps.add(f"euskera.backgrounds.{target}")
        graph[module] = deps

    visiting, visited = set(), set()

    def visit(module):
        if module in visiting:
            raise AssertionError(f"cycle detected at {module}")
        if module in visited:
            return
        visiting.add(module)
        for dependency in graph.get(module, ()):
            visit(dependency)
        visiting.remove(module)
        visited.add(module)

    for module in graph:
        visit(module)
