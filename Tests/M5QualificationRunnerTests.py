#!/usr/bin/env python3
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from compare_m5_reports import compare_m5, LIMITS


class M5ComparatorTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.roots = [Path(self.temp.name) / p for p in ('left', 'right')]
        for root, arch in zip(self.roots, ('arm64', 'x86_64')):
            root.mkdir()
            row = dict(method='m5', checks=[dict(passed=True)], build='build', engine='engine',
                       passed=True, architecture=arch, performance=dict(LIMITS, generated_fixture_seconds=10800, section_count=1000))
            row['measured_section_values'] = [dict(sample_rate=48000, signal_kind=0, begin=0, end=48000,
                **{family: dict(contained=7, crossing=3, padding=0, included=7, missing=0,
                                minimum=-23., median=-23., maximum=-23., share=.5)
                   for family in ('m', 's', 'lra')}) for _ in range(24)]
            (root / 'sections.json').write_text(json.dumps(row))
            (root / 'environment.json').write_text(json.dumps(dict(execution_arch=arch)))

    def compare(self):
        with patch('compare_m5_reports.compare_m4', return_value=dict(errors=[], passed=True)):
            return compare_m5(*self.roots)['passed']

    def change(self, fn):
        p = self.roots[1] / 'sections.json'; row = json.loads(p.read_text()); fn(row); p.write_text(json.dumps(row))

    def test_matching(self): self.assertTrue(self.compare())
    def test_wrong_arch(self):
        self.change(lambda r: r.update(architecture='arm64')); self.assertFalse(self.compare())
    def test_failed_check(self):
        self.change(lambda r: r.update(checks=[dict(passed=False)])); self.assertFalse(self.compare())
    def test_missing_budget(self):
        self.change(lambda r: r['performance'].pop('aggregation_ms')); self.assertFalse(self.compare())
    def test_slow_viewport(self):
        self.change(lambda r: r['performance'].update(interaction_p95_ms=100.0001)); self.assertFalse(self.compare())
    def test_nonfinite(self):
        self.change(lambda r: r['performance'].update(aggregation_ms=float('nan'))); self.assertFalse(self.compare())
    def test_shorter_fixture(self):
        self.change(lambda r: r['performance'].update(generated_fixture_seconds=7200)); self.assertFalse(self.compare())
    def test_negative_timing(self):
        self.change(lambda r: r['performance'].update(aggregation_ms=-1)); self.assertFalse(self.compare())
    def test_missing_sections(self):
        self.change(lambda r: r['performance'].update(section_count=999)); self.assertFalse(self.compare())
    def test_missing_values(self):
        self.change(lambda r: r.pop('measured_section_values')); self.assertFalse(self.compare())
    def test_value_mismatch(self):
        self.change(lambda r: r['measured_section_values'][0]['m'].update(median=-22.999999)); self.assertFalse(self.compare())
    def test_value_tolerance(self):
        self.change(lambda r: r['measured_section_values'][0]['m'].update(median=-23.+1e-10)); self.assertTrue(self.compare())
    def test_membership_mismatch(self):
        self.change(lambda r: r['measured_section_values'][0]['m'].update(crossing=4)); self.assertFalse(self.compare())
    def test_unavailable_mismatch(self):
        self.change(lambda r: r['measured_section_values'][0]['m'].update(median=None)); self.assertFalse(self.compare())
    def test_nonfinite_values(self):
        self.change(lambda r: r['measured_section_values'][0]['m'].update(median=float('nan'))); self.assertFalse(self.compare())
    def test_preserves_regression_failure(self):
        with patch('compare_m5_reports.compare_m4', return_value=dict(errors=['M3 failed'], passed=False)):
            self.assertFalse(compare_m5(*self.roots)['passed'])


if __name__ == '__main__': unittest.main()
