#!/usr/bin/env python3
"""Check local HTML links and fragments in a generated documentation site."""
from __future__ import annotations

from html.parser import HTMLParser
from pathlib import Path
import sys
from urllib.parse import unquote, urljoin, urlparse


class Links(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []
        self.ids: set[str] = set()

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = dict(attrs)
        if values.get("href"):
            self.links.append(values["href"] or "")
        if values.get("id"):
            self.ids.add(values["id"] or "")
        if values.get("name"):
            self.ids.add(values["name"] or "")


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: check_site_links.py SITE_DIR")
    root = Path(sys.argv[1]).resolve()
    pages = list(root.rglob("*.html"))
    parsed = {}
    for page in pages:
        parser = Links()
        parser.feed(page.read_text(encoding="utf-8", errors="replace"))
        parsed[page] = parser
    errors: list[str] = []
    for page, parser in parsed.items():
        base = page.as_uri()
        for link in parser.links:
            target = urlparse(urljoin(base, link))
            if target.scheme not in ("file",):
                continue
            path = Path(unquote(target.path))
            if not path.exists():
                errors.append(f"{page.relative_to(root)}: missing {link}")
            elif target.fragment and path.suffix == ".html":
                dest = path.resolve()
                if dest in parsed and unquote(target.fragment) not in parsed[dest].ids:
                    errors.append(f"{page.relative_to(root)}: missing anchor {link}")
    if errors:
        raise SystemExit("\n".join(errors))
    print(f"Checked local links and anchors across {len(pages)} HTML pages.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
