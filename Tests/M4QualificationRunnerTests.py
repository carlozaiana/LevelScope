#!/usr/bin/env python3
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from compare_m4_reports import compare_m4, LIMITS


class M4ComparatorTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.roots = [Path(self.temp.name) / p for p in ('left', 'right')]
        for root, arch in zip(self.roots, ('arm64', 'x86_64')):
            root.mkdir()
            row = dict(method='m4', checks=[dict(passed=True)], build='build', engine='engine',
                       passed=True, architecture=arch, performance=dict(LIMITS, generated_fixture_seconds=10800))
            (root / 'diagnostics.json').write_text(json.dumps(row))
            (root / 'environment.json').write_text(json.dumps(dict(execution_arch=arch)))

    def compare(self):
        with patch('compare_m4_reports.compare_m3', return_value=dict(errors=[], passed=True)):
            return compare_m4(*self.roots)['passed']

    def change(self, fn):
        p = self.roots[1] / 'diagnostics.json'; row = json.loads(p.read_text()); fn(row); p.write_text(json.dumps(row))

    def test_matching(self): self.assertTrue(self.compare())
    def test_wrong_arch(self):
        self.change(lambda r: r.update(architecture='arm64')); self.assertFalse(self.compare())
    def test_failed_check(self):
        self.change(lambda r: r.update(checks=[dict(passed=False)])); self.assertFalse(self.compare())
    def test_missing_budget(self):
        self.change(lambda r: r['performance'].pop('worker_including_identity_ms')); self.assertFalse(self.compare())
    def test_slow_viewport(self):
        self.change(lambda r: r['performance'].update(viewport_prepare_p95_ms=50.0001)); self.assertFalse(self.compare())
    def test_nonfinite(self):
        self.change(lambda r: r['performance'].update(analysis_ms=float('nan'))); self.assertFalse(self.compare())
    def test_shorter_fixture(self):
        self.change(lambda r: r['performance'].update(generated_fixture_seconds=7200)); self.assertFalse(self.compare())
    def test_preserves_regression_failure(self):
        with patch('compare_m4_reports.compare_m3', return_value=dict(errors=['M3 failed'], passed=False)):
            self.assertFalse(compare_m4(*self.roots)['passed'])


if __name__ == '__main__': unittest.main()
