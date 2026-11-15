import json
import re
from pathlib import Path

from scripts import build_static_demo
def test_static_copy_contains_every_recorded_run(tmp_path):
    summary = build_static_demo.build(tmp_path / 'site')
    site = Path(summary['out'])
    page = (site / 'index.html').read_text(encoding='utf-8')
    assert 'data-static="true"' in page and 'href="app.css"' in page and 'src="app.js"' in page
    assert '/static/' not in page, 'absolute paths break hosting under a sub-path such as GitHub Pages'
    demos = json.loads((site / 'data' / 'demos.json').read_text(encoding='utf-8'))
    assert {d['name'] for d in demos} == {p.parent.name for p in build_static_demo.DEMOS.glob('*/manifest.json')}
    for demo in demos:
        if demo['status'] == 'COMPLETED':
            for name in demo['files']:
                assert (site / 'data' / 'demos' / demo['name'] / name).is_file(), (demo['name'], name)
    assert json.loads((site / 'data' / 'targets.json').read_text(encoding='utf-8'))


def test_static_paths_in_front_end_match_the_build_layout():
    js = Path('static/app.js').read_text(encoding='utf-8')
    assert "'data/targets.json'" in js and "'data/demos.json'" in js
    assert re.search(r'`data/demos/\$\{m\[1\]\}/\$\{m\[2\]\}`', js)
