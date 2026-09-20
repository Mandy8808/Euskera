from pathlib import Path
import subprocess
import sys


def test_markdown_links():
    root = Path(__file__).parents[1]
    result = subprocess.run(
        [sys.executable, "scripts/check_markdown_links.py"],
        cwd=root,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
