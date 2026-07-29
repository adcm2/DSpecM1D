#!/usr/bin/env python3

from pathlib import Path
import re


REPORT_DIRECTORY = Path(__file__).resolve().parent
INVENTORY = REPORT_DIRECTORY.parent / "TESTS.md"
OUTPUT = REPORT_DIRECTORY / "test_inventory.tex"


def escape(text):
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(character, character)
                   for character in text)


def inline(text):
    pieces = re.split(r"(`[^`]+`)", text)
    output = []
    for piece in pieces:
        if piece.startswith("`") and piece.endswith("`"):
            output.append(r"\texttt{" + escape(piece[1:-1]) + "}")
        else:
            output.append(escape(piece))
    return "".join(output)


def main():
    sections = []
    current = None
    normal = 0
    paper = 0
    pattern = re.compile(
        r"^- \*\*(normal|paper-validation)\*\* `([^`]+)` "
        r"— `([^`]+)` — (.+)$"
    )
    for line in INVENTORY.read_text(encoding="utf-8").splitlines():
        if line.startswith("## "):
            current = [line[3:], []]
            sections.append(current)
            continue
        match = pattern.match(line)
        if not match:
            continue
        if current is None:
            raise RuntimeError("Inventory entry appears before a section")
        classification, name, source, description = match.groups()
        current[1].append((classification, name, source, description))
        if classification == "normal":
            normal += 1
        else:
            paper += 1

    if (normal, paper) != (58, 11):
        raise RuntimeError(
            f"Expected 58 normal and 11 paper tests, got {normal} and {paper}"
        )

    lines = [
        "% Generated from ../TESTS.md by generate_test_appendix.py.",
        r"\small",
    ]
    for heading, entries in sections:
        if not entries:
            continue
        lines.extend([
            r"\subsection*{" + escape(heading) + "}",
            r"\begin{itemize}[leftmargin=1.4em,itemsep=0.35em]",
        ])
        for classification, name, source, description in entries:
            label = (
                "normal"
                if classification == "normal"
                else "paper validation"
            )
            lines.extend([
                r"\item {\footnotesize\texttt{" + escape(name) + r"}}\\",
                r"\emph{" + label + r"; \texttt{" + escape(source)
                + r"}.} " + inline(description),
            ])
        lines.append(r"\end{itemize}")
    OUTPUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {normal + paper} entries to {OUTPUT}")


if __name__ == "__main__":
    main()
