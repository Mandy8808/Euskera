from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from euskera.visualization import (
    plotPerf,
    plotUsingDiscSol,
    plotUsingPerf,
    plotUsingSol,
)


def test_plot_using_perf_and_discrete_solution(monkeypatch):
    monkeypatch.setattr(plt, "show", lambda: None)
    r = np.linspace(0.1, 1.0, 4)
    perf = [r, r**2]
    fig, axes = plt.subplots()
    axes = plotUsingPerf(axes, perf)
    assert len(axes) == 1
    fig2, axes2 = plotUsingDiscSol(perf)
    assert fig2 is not None and len(axes2) == 1
    plt.close(fig)
    plt.close(fig2)


def test_plot_using_sol(monkeypatch):
    monkeypatch.setattr(plt, "show", lambda: None)
    solution = SimpleNamespace(
        t=np.linspace(0.0, 1.0, 4),
        y=np.array([[1.0, 0.8, 0.5, 0.2], [0.0, 0.0, 0.0, 0.0],
                    [0.1, 0.1, 0.1, 0.1], [0.0, 0.0, 0.0, 0.0]]),
    )
    fig, axes = plotUsingSol(solution)
    assert fig is not None and len(axes.lines) == 1
    plt.close(fig)


def test_plot_perf_small_case(monkeypatch):
    monkeypatch.setattr(plt, "show", lambda: None)
    fig, axes = plotPerf([1.0, 0.0, 0.1, 0.0], 0.1, [1, 0.0, 0.0])
    assert fig is not None
    assert axes.lines
    plt.close(fig)


def test_background_imports_are_acyclic():
    import euskera.backgrounds
    import euskera.visualization.background_plots
