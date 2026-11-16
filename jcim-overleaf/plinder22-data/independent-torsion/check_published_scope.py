"""Check the published torsion measurements and whether CASF data changed."""
from pathlib import Path
import argparse
import hashlib
import json
import sqlite3

import pandas as pd

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('benchmark', type=Path)
parser.add_argument('--baseline', default='v2-1p7b-coordtoken-v1')
parser.add_argument('--release', default='plinder-torsion-minimized-v1')
args = parser.parse_args()
results = args.benchmark.resolve() / 'data/results'
baseline, release = [results / 'staging' / name for name in [args.baseline, args.release]]


def digest(path):
    with path.open('rb') as stream:
        value = hashlib.sha256()
        for block in iter(lambda: stream.read(1 << 20), b''):
            value.update(block)
        return value.hexdigest()


checks = {}
for line in (results / 'generation/plinder_independent_torsion/casf_snapshot.sha256').read_text().splitlines():
    expected, path = line.split(maxsplit=1)
    if path.startswith('sqlite:'):
        continue
    assert digest(Path(path)) == expected, path
    checks[path] = expected

casf_hashes = [digest(path / 'casf_analysis_dashboard.sqlite') for path in [baseline, release]]
assert casf_hashes[0] == casf_hashes[1]
connections = [sqlite3.connect(f'file:{path}/extended_casf_analysis.sqlite?mode=ro', uri=True)
               for path in [baseline, release]]
tables = [set(row[0] for row in con.execute("SELECT name FROM sqlite_master WHERE type='table'"))
          for con in connections]
assert tables[0] == tables[1]
unchanged = {}
for name in sorted(tables[0]):
    if 'plinder' in name or 'druglike' in name:
        continue
    # Sort complete row representations so physical SQLite row order is immaterial.
    rows = [sorted(repr(row) for row in con.execute(f'SELECT * FROM "{name}"'))
            for con in connections]
    assert rows[0] == rows[1], name
    unchanged[name] = {'rows': len(rows[0]),
                       'sha256': hashlib.sha256('\n'.join(rows[0]).encode()).hexdigest()}

published = {}
for source, public in [('torsion_independent_raw_druglike', 'torsion_raw_druglike'),
                       ('torsion_independent_minimized_druglike', 'torsion_minimized_22_druglike')]:
    published[source] = {}
    for table, filename, keys in [('plinder_summary', 'summary.csv', []),
                                  ('plinder_per_molecule', 'per_molecule.csv', ['smiles']),
                                  ('plinder_recovery_curve', 'recovery.csv', ['threshold']),
                                  ('plinder_k_efficiency', 'k_efficiency.csv', ['k'])]:
        expected = pd.read_csv(results / 'plinder_panel' / source / filename)
        actual = pd.read_sql_query(f'SELECT * FROM "{table}" WHERE label = ?', connections[1], params=[public])
        columns = [c for c in expected.columns if c not in ['label', 'display_label', 'generator', 'recipe']]
        expected, actual = [frame[columns].where(frame[columns].notna(), float("nan")) for frame in [expected, actual]]
        if keys:
            expected, actual = [frame.sort_values(keys) for frame in [expected, actual]]
        pd.testing.assert_frame_equal(expected.reset_index(drop=True), actual.reset_index(drop=True),
                                      check_dtype=False, check_exact=False, atol=1e-12, rtol=0)
        published[source][table] = len(actual)
for con in connections:
    con.close()

report = dict(baseline=str(baseline), release=str(release),
              database_sha256={str(path / name): digest(path / name)
                               for path in [baseline, release]
                               for name in ['casf_analysis_dashboard.sqlite', 'extended_casf_analysis.sqlite']},
              unchanged_casf_files=checks, unchanged_casf_tables=unchanged,
              published_rows_verified=published, structural_recalculation=False)
(HERE / 'published-scope.json').write_text(json.dumps(report, indent=2) + '\n')
print(f'CASF database and {len(unchanged)} extended CASF tables unchanged; both published PLINDER torsion runs verified.')
