"""Build a static, server-free copy of the Epoch Lab that shows the recorded runs.

The output folder can be served by any static host (GitHub Pages, Netlify, a USB
stick opened through a local web server). It contains the same front end, the
preset target list and every file of every completed recorded run; live
analysis is disabled because it needs the Python server and the NASA archives.

    python build_static_demo.py            # writes ./site
    python build_static_demo.py --out DIR
"""
import argparse
import json
import shutil
from pathlib import Path

from targets import targets

ROOT = Path(__file__).resolve().parent
STATIC = ROOT / 'static'
DEMOS = ROOT / 'artifacts' / 'epoch_lab' / 'demos'
SERVED_SUFFIXES = {'.json', '.png', '.fits', '.csv'}  # same rule as the API's file endpoint


def build(out):
    out = Path(out)
    if out.exists():
        shutil.rmtree(out)
    (out / 'data' / 'demos').mkdir(parents=True)
    page = (STATIC / 'index.html').read_text(encoding='utf-8')
    replacements = (('<html lang="en">', '<html lang="en" data-static="true">'),
                    ('href="/static/app.css"', 'href="app.css"'), ('src="/static/app.js"', 'src="app.js"'))
    for old, new in replacements:
        if page.count(old) != 1:
            raise ValueError(f'static/index.html no longer contains {old!r} exactly once')
        page = page.replace(old, new)
    (out / 'index.html').write_text(page, encoding='utf-8')
    for name in ('app.css', 'app.js'):
        shutil.copy2(STATIC / name, out / name)
    (out / 'data' / 'targets.json').write_text(json.dumps(targets()), encoding='utf-8')

    manifests, copied = [], 0
    for path in sorted(DEMOS.glob('*/manifest.json')):
        manifest = json.loads(path.read_text(encoding='utf-8'))
        manifests.append(manifest)
        if manifest.get('status') != 'COMPLETED':
            continue
        target = out / 'data' / 'demos' / manifest['name']
        target.mkdir()
        for name in manifest.get('files', []) + ['viewer.json']:
            source = path.parent / name
            if source.suffix in SERVED_SUFFIXES and source.is_file() and source.parent == path.parent:
                shutil.copy2(source, target / name)
                copied += 1
    (out / 'data' / 'demos.json').write_text(json.dumps(manifests), encoding='utf-8')
    (out / '.nojekyll').write_text('', encoding='utf-8')  # GitHub Pages: serve files verbatim
    size = sum(f.stat().st_size for f in out.rglob('*') if f.is_file())
    return {'out': str(out), 'runs': len(manifests), 'files': copied, 'bytes': size}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', default=str(ROOT / 'site'), help='output folder (replaced if it exists)')
    summary = build(parser.parse_args().out)
    print(f"{summary['runs']} recorded runs, {summary['files']} files, {summary['bytes'] / 1e6:.1f} MB written to {summary['out']}")


if __name__ == '__main__':
    main()
