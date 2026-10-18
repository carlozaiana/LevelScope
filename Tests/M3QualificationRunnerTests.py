#!/usr/bin/env python3
"""Exercise the additional M3 comparator guards; M2 comparison remains separate."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from compare_m3_reports import compare_m3


class ComparatorTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.left, self.right = (Path(self.temp.name) / n for n in ('left', 'right'))
        self.base = dict(method='test-method', checks=[dict(passed=True)], profile_digest='profile',
                         build='build', engine='engine', scalar_outcomes=[dict(total_gain=3)], passed=True,
                         architecture='arm64')
        for root, arch in ((self.left, 'arm64'), (self.right, 'x86_64')):
            root.mkdir()
            row = dict(self.base, architecture=arch)
            (root / 'proposal.json').write_text(json.dumps(row))
            (root / 'environment.json').write_text(json.dumps(dict(execution_arch=arch)))

    def run_comparison(self):
        with patch('compare_m3_reports.compare', return_value=dict(errors=[], passed=True)):
            return compare_m3(self.left, self.right)

    def alter(self, **fields):
        path = self.right / 'proposal.json'
        data = json.loads(path.read_text()); data.update(fields); path.write_text(json.dumps(data))

    def test_matching(self):
        self.assertTrue(self.run_comparison()['passed'])

    def test_plan_difference(self):
        self.alter(scalar_outcomes=[dict(total_gain=3.0000000000000004)])
        self.assertFalse(self.run_comparison()['passed'])

    def test_method_difference(self):
        self.alter(method='other-version')
        self.assertFalse(self.run_comparison()['passed'])

    def test_failure(self):
        self.alter(passed=False)
        self.assertFalse(self.run_comparison()['passed'])

    def test_architecture(self):
        self.alter(architecture='arm64')
        self.assertFalse(self.run_comparison()['passed'])

    def test_missing_outcomes(self):
        path = self.right / 'proposal.json'
        data = json.loads(path.read_text()); del data['scalar_outcomes']; path.write_text(json.dumps(data))
        self.assertFalse(self.run_comparison()['passed'])

    def test_preserves_existing_comparison_failure(self):
        with patch('compare_m3_reports.compare', return_value=dict(errors=['measurement failure'], passed=False)):
            self.assertFalse(compare_m3(self.left, self.right)['passed'])


if __name__ == '__main__':
    unittest.main()
