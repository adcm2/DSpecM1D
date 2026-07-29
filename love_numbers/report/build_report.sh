#!/usr/bin/env sh

set -eu

report_directory=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$report_directory"

python3 generate_test_appendix.py

if command -v latexmk >/dev/null 2>&1; then
    latexmk -pdf -interaction=nonstopmode -halt-on-error \
        love_numbers_report.tex
else
    pdflatex -interaction=nonstopmode -halt-on-error \
        love_numbers_report.tex
    if command -v bibtex >/dev/null 2>&1; then
        bibtex love_numbers_report
    fi
    pdflatex -interaction=nonstopmode -halt-on-error \
        love_numbers_report.tex
    pdflatex -interaction=nonstopmode -halt-on-error \
        love_numbers_report.tex
fi
