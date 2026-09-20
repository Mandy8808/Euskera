# Development and validation

Install the editable package with the test extra:

```bash
python -m pip install -e ".[test]"
```

Run the focused and scientific tests from the repository root:

```bash
python -m pytest
python -m pytest -m scientific
python -m compileall euskera
```

The lightweight Markdown link check is:

```bash
python scripts/check_markdown_links.py
```

It scans Markdown files, ignores external URLs and anchors, and verifies that
relative file targets exist. Documentation changes should preserve the
dictionary-based model API, link to notebooks without moving them, and avoid
committing generated simulation output. The CI workflow is documented in
[`.github/workflows/ci.yml`](../.github/workflows/ci.yml).
