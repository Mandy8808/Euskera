"""Check relative links in Markdown files without third-party dependencies."""

from __future__ import annotations

import re
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
LINK = re.compile(r'(?<!!)\[[^\]]+\]\(([^)\s]+)(?:\s+[\'"][^)]*)?\)')


def main() -> int:
    errors: list[str] = []
    for markdown in sorted(ROOT.rglob("*.md")):
        if any(part in {".git", ".venv"} for part in markdown.parts):
            continue
        for target in LINK.findall(markdown.read_text(encoding="utf-8")):
            parsed = urlsplit(target)
            if parsed.scheme or parsed.netloc or target.startswith("#"):
                continue
            path = (markdown.parent / unquote(parsed.path)).resolve()
            if not path.exists():
                errors.append(f"{markdown.relative_to(ROOT)}: {target}")
    if errors:
        print("Broken Markdown links:")
        print("\n".join(errors))
        return 1
    print("Markdown links: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
