#!/usr/bin/env python3
"""Acceptance policy regressions; synthetic counters, no devices or media."""
import copy
import unittest
from m7_video_acceptance import evaluate


def report():
    r = dict(passed=True, error='', all_routes_muted=True, selected_audio_channels=8,
        layout_capability='pcm-7.1-bs2051-system-i-v1',pre_monitor_nonzero_channel_mask=255,
        pre_monitor_channel_count_valid=True,blank_paints=0,bound_failures=0,
        paint_deadline_failures=0,underrun_blocks=0,attachment_within_budget=True,
        memory_passed=True,hal_output_timestamp=True,skipped_frame_limit=0,
        source_seconds=18.01,requested_seconds=18,maximum_paint_interval_ms=17,
        maximum_frame_pts_skew_seconds=0,maximum_cpu_paint_ms=.2,mean_skew_seconds=0,
        first_minute_mean_skew_seconds=0,last_minute_mean_skew_seconds=0,
        actual_paint_skipped_frames=0,unique_painted_frames=400,paints=960,
        pre_monitor_frames_observed=864000,audition_identity='current:identity',skip_events=[],
        measurement_build='expected',memory_method='prores-memory-v2-user-approved-2026-10-20',
        qualification_eligible=True,allocator_relief_requested=False,memory_query_failures=0,
        first_minute_memory_sample_count=17,last_minute_memory_sample_count=17)
    r['memory_samples']=[dict(seconds=i,rss_bytes=1000,reusable_bytes=200,
                             non_reusable_resident_bytes=800,physical_footprint_bytes=900) for i in range(2,19)]
    for field,value in [('rss',1000),('non_reusable_resident',800),('physical_footprint',900)]:
        r['peak_'+('sampled_playback_rss' if field=='rss' else field)+'_bytes']=value
        r['first_minute_mean_'+field+'_bytes']=value
        r['last_minute_mean_'+field+'_bytes']=value
    return r


class Acceptance(unittest.TestCase):
    def check(self,r,rc=0): return evaluate(r,{'measurement_build':'expected'},rc,25)
    def test_strict_pass_and_input_unchanged(self):
        r=report(); before=copy.deepcopy(r)
        self.assertTrue(self.check(r)['accepted']);self.assertEqual(before,r)
    def test_skip_is_separate_accepted_limitation(self):
        r=report();r.update(passed=False,actual_paint_skipped_frames=1,unique_painted_frames=399,
            maximum_paint_interval_ms=51,skip_events=[dict(frame=12,previous_frame=10,audio_seconds=3.0,paint_interval_ms=51)])
        before=copy.deepcopy(r);v=self.check(r,1)
        self.assertTrue(v['accepted']);self.assertFalse(v['strict_harness_passed']);self.assertEqual(before,r)
    def test_protected_failures(self):
        for key,value in [('underrun_blocks',1),('bound_failures',1),('blank_paints',1),
            ('paint_deadline_failures',1),('maximum_paint_interval_ms',251),
            ('pre_monitor_nonzero_channel_mask',127),('selected_audio_channels',2),
            ('measurement_build','other'),('hal_output_timestamp',False),('all_routes_muted',False),
            ('source_seconds',17),('audition_identity',''),('last_minute_mean_skew_seconds',.1),
            ('memory_passed',False),('qualification_eligible',False),('error','failure')]:
            with self.subTest(key=key):
                r=report();r[key]=value;self.assertFalse(self.check(r)['accepted'])
    def test_missing_and_nonfinite(self):
        for key in ('source_seconds','maximum_frame_pts_skew_seconds','maximum_paint_interval_ms','last_minute_mean_skew_seconds'):
            for value in (None,float('nan'),float('inf')):
                r=report();r[key]=value;self.assertFalse(self.check(r)['accepted'])
    def test_raw_failure_without_skips_is_not_relabelled(self):
        r=report();r['passed']=False;self.assertFalse(self.check(r,1)['accepted'])
        self.assertFalse(self.check(report(),1)['accepted'])
    def test_skip_accounting_and_coverage(self):
        r=report();r.update(passed=False,actual_paint_skipped_frames=2)
        self.assertFalse(self.check(r,1)['accepted'])
        r=report();r['unique_painted_frames']=10;self.assertFalse(self.check(r)['accepted'])
    def test_memory_recomputed(self):
        for key,value in [('rss_bytes',900*1024**2),('non_reusable_resident_bytes',2),('physical_footprint_bytes',float('nan'))]:
            r=report();r['memory_samples'][0][key]=value
            self.assertFalse(self.check(r)['accepted'])
    def test_endurance_requires_growth_evidence(self):
        r=report();r.update(source_seconds=1800,requested_seconds=1800,unique_painted_frames=44950)
        self.assertFalse(self.check(r)['accepted'])


if __name__=='__main__': unittest.main()
