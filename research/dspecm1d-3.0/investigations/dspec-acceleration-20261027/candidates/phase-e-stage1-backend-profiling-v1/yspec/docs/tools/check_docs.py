#!/usr/bin/env python3
"""Check local HTML links, anchors, assets and authored code block wrappers."""
import argparse
import ast
import copy
import hashlib
import json
import re
import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

PLOT_ASSETS = (
    'docs/site/assets/first-seismogram.svg',
    'docs/site/assets/traces-small.png',
    'docs/site/assets/kernels-radial-small.png',
    'docs/site/assets/kernels-section-small.png',
)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def plot_provenance_errors(manifest, root, site):
    """Check current source hashes and the built copies of tracked plot assets."""
    errors = []
    figures = manifest.get('figures', [])
    ids = [figure.get('id') for figure in figures]
    if len(figures) != 21 or len(set(ids)) != 21:
        errors.append('plot provenance must contain 21 unique figure IDs')
    incomplete = [figure.get('id', '<unnamed>') for figure in figures
                  if figure.get('render_status') != 'rendered']
    if incomplete:
        errors.append('plot provenance has incomplete renders: ' +
                      ', '.join(incomplete))

    def verify_hash(label, relative, expected):
        path = (root / relative).resolve()
        try:
            path.relative_to(root.resolve())
        except ValueError:
            errors.append(f'{label}: source path escapes repository: {relative}')
            return
        try:
            actual = sha256(path)
        except OSError as error:
            errors.append(f'{label}: cannot read {relative}: {error}')
            return
        if actual != expected:
            errors.append(f'{label}: hash mismatch for {relative}: expected '
                          f'{expected}, found {actual}')

    for figure in figures:
        ident = figure.get('id', '<unnamed>')
        for relative, expected in figure.get('input_hashes', {}).items():
            verify_hash(f'{ident} input', relative, expected)
        model = figure.get('model')
        if model and figure.get('model_sha256'):
            verify_hash(f'{ident} model', model, figure['model_sha256'])
        source = figure.get('plotter_source')
        if source and figure.get('plotter_sha256'):
            verify_hash(f'{ident} plotter', source, figure['plotter_sha256'])
        for relative, expected in figure.get('generator_hashes', {}).items():
            verify_hash(f'{ident} generator', relative, expected)

    for label in ('shared_plot_helpers', 'shared_plotter_modules'):
        for relative, expected in manifest.get(label, {}).items():
            verify_hash(label, relative, expected)

    tracked = {(figure.get('render_output') or {}).get('tracked_path'): figure
               for figure in figures
               if (figure.get('render_output') or {}).get('tracked_path')}
    if set(tracked) != set(PLOT_ASSETS):
        errors.append('plot provenance must list exactly the four tracked site assets')
    for relative in PLOT_ASSETS:
        figure = tracked.get(relative)
        if not figure:
            continue
        output_hash = figure.get('render_output', {}).get('sha256')
        source = (root / relative).resolve()
        try:
            actual_source = sha256(source)
        except OSError as error:
            errors.append(f'plot asset: cannot read {relative}: {error}')
            continue
        if actual_source != output_hash:
            errors.append(f'plot asset hash mismatch for {relative}: expected '
                          f'{output_hash}, found {actual_source}')
        built_relative = Path(relative).relative_to('docs/site')
        built = site / built_relative
        try:
            actual_built = sha256(built)
        except OSError as error:
            errors.append(f'built plot asset: cannot read {built_relative}: {error}')
            continue
        if actual_built != output_hash:
            errors.append(f'built plot asset hash mismatch for {built_relative}: '
                          f'expected {output_hash}, found {actual_built}')
    return errors



def markdown_output_hash_errors(path, expected, label):
    """Compare a current Markdown output table with its JSON hash source."""
    try:
        lines = path.read_text().splitlines()
    except OSError as error:
        return [f'{label}: cannot read Markdown table: {error}']
    found = {}
    errors = []
    for line in lines:
        if not line.lstrip().startswith('|'):
            continue
        cells = line.strip().strip('|').split('|')
        if not cells:
            continue
        match = re.search(r'`([^`]+)`', cells[0])
        if not match:
            continue
        name = match.group(1)
        if name not in expected:
            continue
        hashes = re.findall(r'(?<![0-9a-f])[0-9a-f]{64}(?![0-9a-f])', line)
        if len(hashes) != 1:
            errors.append(f'{label}: expected one SHA-256 in row {name}')
        elif name in found:
            errors.append(f'{label}: duplicate output row {name}')
        else:
            found[name] = hashes[0]
    for name, digest in expected.items():
        if name not in found:
            errors.append(f'{label}: missing output row {name}')
        elif found[name] != digest:
            errors.append(f'{label}: output hash mismatch for {name}: expected '
                          f'{digest}, found {found[name]}')
    return errors


def check_plot_provenance(site, root):
    path = site / 'examples/all-plots-provenance.json'
    try:
        manifest = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        return [f'plot provenance is invalid or missing: {error}']
    errors = plot_provenance_errors(manifest, root, site)

    markdown_outputs = {
        figure.get('id'): (figure.get('render_output') or {}).get('sha256')
        for figure in manifest.get('figures', [])
    }
    errors.extend(markdown_output_hash_errors(
        site / 'examples/all-plots-provenance.md', markdown_outputs,
        'all-plots Markdown table'))
    phase2c_assets = {
        (figure.get('render_output') or {}).get('tracked_path'):
        (figure.get('render_output') or {}).get('sha256')
        for figure in manifest.get('figures', [])
        if (figure.get('render_output') or {}).get('tracked_path')
    }
    errors.extend(markdown_output_hash_errors(
        site / 'examples/first-seismogram-provenance.md', phase2c_assets,
        'Phase 2C asset Markdown table'))

    # Prove that the guard catches a stale source hash without altering the
    # built manifest or any repository file.
    probe = copy.deepcopy(manifest)
    first = probe['figures'][0]
    input_path = next(iter(first.get('input_hashes', {})), None)
    if input_path:
        first['input_hashes'][input_path] = '0' * 64
        probe_errors = plot_provenance_errors(probe, root, site)
        if not any(' input: hash mismatch' in error for error in probe_errors):
            errors.append('plot provenance stale-hash negative probe was not rejected')
    else:
        errors.append('plot provenance stale-hash negative probe has no input to alter')
    return errors

class Page(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.ids, self.refs, self.pre_depth, self.pre_wrapped = set(), [], 0, []
        self.current_count = 0
    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if 'id' in attrs: self.ids.add(attrs['id'])
        if tag == 'a' and 'name' in attrs: self.ids.add(attrs['name'])
        if tag == 'a' and 'aria-current' in attrs and attrs['aria-current'] == 'page': self.current_count += 1
        if tag in ('a', 'link', 'script', 'img'):
            key = 'href' if tag in ('a', 'link') else 'src'
            if attrs.get(key): self.refs.append((tag, attrs[key]))
        if tag == 'pre':
            self.pre_depth += 1
            self.pre_wrapped.append('code-block' in attrs.get('class', '').split())
    def handle_endtag(self, tag):
        if tag == 'pre' and self.pre_depth:
            self.pre_depth -= 1

def inspect(html):
    page = Page(); page.feed(html.read_text(errors='replace')); return page

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('site', type=Path, help='built docs/site directory')
    args = parser.parse_args()
    site = args.site.resolve()
    if not (site / 'index.html').is_file() or not (site / 'api/index.html').is_file():
        parser.error('expected built site/index.html and site/api/index.html')
    html_files = sorted(site.rglob('*.html'))
    parsed = {p: inspect(p) for p in html_files}
    errors = []
    required_pages = {
        'index.html', 'installation.html', 'first-seismogram.html',
        'python.html', 'plotting.html', 'kernels.html',
        'library-example.html', 'examples.html', 'input-files.html',
        'model-and-output-formats.html', 'methods-and-validation.html',
        'release-notes.html', 'maintainers.html',
    }
    missing_pages = sorted(required_pages - {p.relative_to(site).as_posix() for p in html_files})
    if missing_pages:
        errors.append('missing authored pages: ' + ', '.join(missing_pages))
    required_api = {
        'api/Input_8hpp.html', 'api/Instrument_8hpp.html',
        'api/Run_8hpp.html', 'api/Synthetics_8hpp.html',
        'api/Kernels_8hpp.html', 'api/classYSpec_1_1Model.html',
        'api/classYSpec_1_1RadialSolver.html',
    }
    authored_refs = {
        ref for path, page in parsed.items()
        if 'api' not in path.relative_to(site).parts
        for tag, ref in page.refs if tag == 'a'
    }
    missing_api = sorted(required_api - {urlsplit(ref).path for ref in authored_refs})
    if missing_api:
        errors.append('authored guides must link generated API pages: ' + ', '.join(missing_api))
    required_api_fragments = {
        'api/Input_8hpp.html#a159b2668249241093472863aae07800b',
        'api/Instrument_8hpp.html#ae4d3fd526efb1cdb3c82d9965a3f22d4',
        'api/Run_8hpp.html#a005fd5bf5085a7ce92b77071a426df19',
        'api/Synthetics_8hpp.html#a37f382b0fc276cb40ecc2f2a948e101c',
        'api/Kernels_8hpp.html#a87ca499eeb4cdae1c1060626cd687760',
        'api/Input_8hpp.html#a9f0b0ad165807071cecf1808d796bca3',
        'api/Run_8hpp.html#a94dff6e2a13c21656f02548b0db61ba2',
        'api/Kernels_8hpp.html#a2700b3e83da66c395e2b58f31a290f40',
        'api/classYSpec_1_1Model.html#a11f4503d9302dc5a1f79020f2cc6b72b',
    }
    actual_api_fragments = {
        f'{urlsplit(ref).path}#{urlsplit(ref).fragment}'
        for ref in authored_refs if urlsplit(ref).fragment
    }
    missing_fragments = sorted(required_api_fragments - actual_api_fragments)
    if missing_fragments:
        errors.append('authored guides must link exact generated API symbols: ' + ', '.join(missing_fragments))
    for page_path, page in parsed.items():
        for tag, ref in page.refs:
            url = urlsplit(ref)
            if url.scheme or url.netloc or ref.startswith('//'):
                continue
            target = (page_path.parent / unquote(url.path)).resolve() if url.path else page_path
            try: target.relative_to(site)
            except ValueError:
                errors.append(f'{page_path.relative_to(site)}: local URL escapes site: {ref}')
                continue
            if not target.is_file():
                errors.append(f'{page_path.relative_to(site)}: missing {tag} target {ref}')
                continue
            fragment = unquote(url.fragment)
            if fragment and target.suffix.lower() == '.html':
                target_page = parsed.get(target) or inspect(target)
                if fragment not in target_page.ids:
                    errors.append(f'{page_path.relative_to(site)}: missing anchor {ref}')
        if 'api' not in page_path.relative_to(site).parts:
            rel = page_path.relative_to(site)
            source = page_path.read_text(errors='replace')
            if '{{' in source or '}}' in source:
                errors.append(f'{rel}: unresolved build placeholder remains')
            if page.current_count != 1:
                errors.append(f'{rel}: expected exactly one aria-current="page" navigation link')
            if page.pre_wrapped and not all(page.pre_wrapped):
                errors.append(f'{rel}: each authored pre block must use the copy-code wrapper')
            if not (site / 'assets/site.js').is_file():
                errors.append(f'{rel}: shared code-copy script is missing')
    js = (site / 'assets/site.js').read_text()
    for behavior in ('navigator.clipboard.writeText', 'document.execCommand', 'aria-live'):
        if behavior not in js: errors.append(f'copy script lacks {behavior}')
    try:
        json.loads((site / 'examples/first-seismogram.json').read_text())
    except (OSError, json.JSONDecodeError) as error:
        errors.append(f'tutorial JSON fixture is invalid or missing: {error}')
    if not (site / 'assets/first-seismogram.svg').is_file():
        errors.append('tutorial plot asset is missing')
    try:
        json.loads((site / 'examples/snapshot-small.json').read_text())
    except (OSError, json.JSONDecodeError) as error:
        errors.append(f'snapshot JSON fixture is invalid or missing: {error}')
    try:
        json.loads((site / 'examples/workflows.json').read_text())
    except (OSError, json.JSONDecodeError) as error:
        errors.append(f'workflow JSON fixture is invalid or missing: {error}')
    for fixture in ('kernels-small.json', 'kernels-volume-small.json'):
        try:
            json.loads((site / 'examples' / fixture).read_text())
        except (OSError, json.JSONDecodeError) as error:
            errors.append(f'{fixture} is invalid or missing: {error}')
    python_example = site / 'examples/python-small-calculation.py'
    try:
        ast.parse(python_example.read_text(), filename=str(python_example))
    except (OSError, SyntaxError) as error:
        errors.append(f'Python calculation example is invalid or missing: {error}')
    for required_asset in ('examples/prem.200', 'examples/station.pz',
                           'examples/workflows.json',
                           'examples/library-workflows.cpp',
                           'examples/kernels-small.json',
                           'examples/kernels-volume-small.json',
                           'examples/python-small-calculation.py',
                           'assets/traces-small.png',
                           'assets/kernels-radial-small.png',
                           'assets/kernels-section-small.png',
                           'assets/adjoint-kernels.pdf'):
        if not (site / required_asset).is_file():
            errors.append(f'required downloadable asset is missing: {required_asset}')
    try:
        fixture_hash = hashlib.sha256((site / 'examples/prem.200').read_bytes()).digest()
        canonical_hash = hashlib.sha256((Path(__file__).resolve().parents[2] / 'examples/prem.200').read_bytes()).digest()
        if fixture_hash != canonical_hash:
            errors.append('published PREM model copy differs from examples/prem.200')
    except OSError as error:
        errors.append(f'could not verify the published PREM model copy: {error}')
    root = Path(__file__).resolve().parents[2]
    errors.extend(check_plot_provenance(site, root))
    if errors:
        print('\n'.join(errors), file=sys.stderr); return 1
    print(f'Checked {len(html_files)} HTML pages, local links, anchors, assets, '
          'code-copy wrappers, 21-figure provenance and Markdown hash tables '
          '(including stale-hash probes).')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
