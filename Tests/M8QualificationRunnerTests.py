#!/usr/bin/env python3
"""Qualification must fail closed on missing/wrong-platform or partial evidence."""
import copy
import unittest

from run_m8_qualification import validate_results


class RunnerGuards(unittest.TestCase):
    def setUp(self):
        self.identity = dict(source_digest='candidate', measurement_build='build')
        self.env = dict(execution_arch='arm64', platform_guard_passed=True,
            product_version=dict(output='15.7.9'), artifact_integrity_errors=[], source_digest='candidate')
        self.core = dict(passed=True, build='build', architecture='arm64',
            performance=dict(capture_ms=1, ui_adoption_ms=1, save_ms=100, parse_ms=100,
                active_cancel_ms=1, rss_growth_bytes=1000, session_bytes=1000,
                ui_callback_p95_ms=10, ui_callback_count=20),
            long_media=dict(frames=518400000, rate=48000, channels=8, active_cancel_ms=1, cancel_progress=.01))
        self.playback = dict(passed=True, build='build', architecture='arm64', reference_cases=5,
            physical_devices_opened=False, checks=[dict(passed=True)], sample_cases=[
                dict(rate=r, channels=c, gain_db=g) for r in (44100,48000,96000,192000)
                for c in (2,6,8) for g in (0.,-6.020599913279624,6.020599913279624,24.)])

    def verify(self):
        validate_results(self.core, self.playback, self.identity, self.env)

    def test_complete(self):
        self.verify()

    def test_incomplete_or_mismatched(self):
        for target, key, value in [
            ('env','platform_guard_passed',False), ('env','source_digest','other'),
            ('env','artifact_integrity_errors',['bad']), ('core','passed',False),
            ('core','build','other'), ('core','architecture','x86_64'),
            ('playback','reference_cases',4), ('playback','physical_devices_opened',True),
            ('playback','sample_cases',[]), ('playback','checks',[])]:
            with self.subTest(target=target, key=key):
                saved = copy.deepcopy(getattr(self,target))
                getattr(self,target)[key] = value
                with self.assertRaises(ValueError): self.verify()
                setattr(self,target,saved)

    def test_intel_on_sequoia_is_not_high_sierra(self):
        self.env['execution_arch']='x86_64'
        self.core['architecture']=self.playback['architecture']='x86_64'
        with self.assertRaises(ValueError): self.verify()

    def test_resource_budgets(self):
        for name in self.core['performance']:
            with self.subTest(name=name):
                previous=self.core['performance'][name]
                self.core['performance'][name]=0 if name=='ui_callback_count' else 10**12
                with self.assertRaises(ValueError): self.verify()
                self.core['performance'][name]=previous

    def test_duplicate_capture_case(self):
        self.playback['sample_cases'][-1]=self.playback['sample_cases'][0]
        with self.assertRaises(ValueError): self.verify()

    def test_missing_long_verification(self):
        self.core['long_media']['frames']=1000
        with self.assertRaises(ValueError): self.verify()


if __name__ == '__main__':
    unittest.main()
