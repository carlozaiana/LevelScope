"""Check that qualification rejects bad/missing evidence and asymmetric errors."""
import copy
import unittest
from run_ebu_measurements import evaluate, compare_results


class AcceptanceTests(unittest.TestCase):
    def rule(self, kind='scalar', **extra):
        return dict(kind=kind, metric='TP', expected=0, lower_error=-.4, upper_error=.2, **extra)

    def test_asymmetric_limits(self):
        for value, accepted in [(-.4, True), (.2, True), (-.400001, False), (.200001, False)]:
            self.assertEqual(evaluate(self.rule(), {'TP': value})['passed'], accepted)

    def test_nonfinite_and_missing_readings(self):
        for value in [None, float('nan'), float('inf'), True]:
            self.assertFalse(evaluate(self.rule(), {'TP': value})['passed'])

    def test_curve_requires_every_reading(self):
        rule = self.rule('curve', start_frame=10, end_frame=30, step_frames=10)
        points = [{'end_frame': n, 'value': 0} for n in [10, 20, 30]]
        self.assertTrue(evaluate(rule, {'TP_curve': points})['passed'])
        points[1]['value'] = .3
        self.assertFalse(evaluate(rule, {'TP_curve': points})['passed'])

    def test_missing_duplicate_and_unordered_curve(self):
        rule = self.rule('curve', start_frame=10, end_frame=30, step_frames=10)
        for frames in [[], [10, 30], [10, 20, 20, 30], [30, 20, 10]]:
            with self.assertRaises(ValueError):
                evaluate(rule, {'TP_curve': [{'end_frame': n, 'value': 0} for n in frames]})

    def test_checkpoint_requires_each_expected_level(self):
        rule = self.rule('checkpoints', period_frames=10)
        rule['expected'] = [-2, -1, 0]
        points = [{'end_frame': n, 'TP': v} for n, v in [(10, -2), (20, -1), (30, 0)]]
        self.assertTrue(evaluate(rule, {'checkpoints': points})['passed'])
        points[1]['TP'] = 0
        self.assertFalse(evaluate(rule, {'checkpoints': points})['passed'])
        with self.assertRaises(ValueError):
            evaluate(rule, {'checkpoints': points[:-1]})

    def test_unknown_rule_fails(self):
        with self.assertRaises(ValueError):
            evaluate(self.rule('unknown'), {})

    def test_comparison_checks_interior_curve_values(self):
        metrics = dict(I=0, M=0, S=0, TP=0, SP=0, LRA=0,
                       M_curve=[dict(end_frame=1, value=-2), dict(end_frame=2, value=-1), dict(end_frame=3, value=0)])
        current = dict(file='fixture', sha256='hash', passed=True, measured=metrics, checks=[])
        previous = dict(results=[copy.deepcopy(current)])
        self.assertEqual(compare_results([current], previous), [])
        current['measured']['M_curve'][1]['value'] = -.9
        self.assertTrue(compare_results([current], previous))

    def test_comparison_rejects_failed_or_changed_fixture(self):
        current = dict(file='fixture', sha256='hash', passed=True)
        self.assertTrue(compare_results([current], dict(results=[])))
        self.assertTrue(compare_results([current], dict(results=[dict(file='fixture', sha256='other', passed=True)])))


if __name__ == '__main__':
    unittest.main()
