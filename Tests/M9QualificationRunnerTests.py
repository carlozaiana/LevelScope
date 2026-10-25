#!/usr/bin/env python3
"""Fail-closed M9 completeness, architecture, candidate and resource checks."""
import copy
import unittest
from run_m9_qualification import validate_results
from compare_m9_extended_oracles import validate as validate_extended, compare as compare_extended


class Guards(unittest.TestCase):
    def setUp(self):
        self.identity=dict(source_digest='source', measurement_build='build')
        self.env=dict(execution_arch='arm64',platform_guard_passed=True,product_version=dict(output='15.7.9'),artifact_integrity_errors=[],source_digest='source')
        self.report=dict(passed=True,build='build',architecture='arm64',physical_devices_opened=False,checks=[dict(passed=True)],reference_cases=10,cases=[])
        self.report['cases']=[dict(rate=r,channels=c,gain_db=g,kind='M9 mono group Current') for r in (44100,48000,96000,192000) for c in (2,6) for g in (0.,6.020599913279624,-6.020599913279624,-3.305074118944,-60.,24.)]
        self.report['cases'] += [dict(rate=r,channels=c,gain_db=0.,kind='M9 mono group Source') for r in (44100,48000,96000,192000) for c in (2,6)]
        self.report['m8_migrations']=[dict(channels=c,original_unchanged=True) for c in (2,6,8)]
        self.long=dict(passed=True,build='build',architecture='arm64',physical_devices_opened=False,checks=[dict(passed=True)],long_programme=dict(frames=518400000,rate=48000,channels=6,save_ms=1,parse_ms=1,cancel_ms=1,preparation_cancel_ms=1,max_seek_prefill_ms=1,ui_callback_p95_ms=1,rss_growth_bytes=100,session_bytes=100,ui_callback_count=10,owned_payload_bytes=1000,original_payload_bytes=1000))
    def verify(self):validate_results(self.report,self.long,self.identity,self.env)
    def test_complete(self):self.verify()
    def test_missing_or_mismatched(self):
        for name,key,value in [('env','platform_guard_passed',False),('env','artifact_integrity_errors',['bad']),('env','source_digest','other'),('report','passed',False),('report','build','other'),('report','architecture','x86_64'),('report','checks',[]),('report','checks',[dict(passed=False)]),('report','reference_cases',9),('report','cases',[]),('report','physical_devices_opened',True),('long','passed',False)]:
            with self.subTest(name=name,key=key):
                saved=copy.deepcopy(getattr(self,name));getattr(self,name)[key]=value
                with self.assertRaises(ValueError):self.verify()
                setattr(self,name,saved)
    def test_intel_on_sequoia_rejected(self):
        self.env['execution_arch']='x86_64';self.report['architecture']=self.long['architecture']='x86_64'
        with self.assertRaises(ValueError):self.verify()
    def test_duplicate_case(self):
        self.report['cases'][-1]=self.report['cases'][0]
        with self.assertRaises(ValueError):self.verify()
    def test_incomplete_migration(self):
        self.report['m8_migrations'].pop()
        with self.assertRaises(ValueError):self.verify()
    def test_resource_budgets(self):
        media=self.long['long_programme']
        for key in media:
            with self.subTest(key=key):
                previous=media[key];media[key]=0 if key=='ui_callback_count' else 10**12
                with self.assertRaises(ValueError):self.verify()
                media[key]=previous


class ExtendedGuards(unittest.TestCase):
    def setUp(self):
        # Synthetic reports test rejection logic only; never qualification evidence.
        rows=[]
        for channels in (2,6):
            for rate in (44100,48000,96000,192000):
                for silence in (False,True):
                    for domain in ('Source','Current State'):
                        frames=rate*(4 if silence else 64)+17;step=rate//10
                        windows=[[family,n,0. if silence or i%2 else .01] for family,start,end in
                                 ((0,4*step,frames),(1,30*step,frames),(2,30*step,frames+15*step))
                                 for i,n in enumerate(range(start,end+1,step))]
                        rows.append(dict(channels=channels,rate=rate,frames=frames,silence=silence,domain=domain,
                            gain_db=0 if domain=='Source' else -3.305074118944,lra_stable=not silence,
                            integrated_included=0 if silence else 1,lra_included=0 if silence else 1,
                            sample_peak_frames=[0]*channels,true_peak_eighth_frames=[0]*channels,
                            sample_peaks=[0. if silence else .1]*channels,true_peaks=[0. if silence else .1]*channels,
                            windows=windows,integrated_gate=None if silence else -40.,lra_gate=None if silence else -50.,
                            p10=None if silence else -30.,p95=None if silence else -20.,same_platform_full_result_exact=True,
                            metrics={name:dict(status=1 if silence else 0,value=None if silence else 10. if name=='lra' else -20.)
                                     for name in ('integrated','lra','momentary','short_term','sample_peak','true_peak')}))
        self.left=dict(architecture='arm64',os='macOS 15.7.9',passed=True,physical_devices_opened=False,
                       checks=[dict(passed=True)],extended_oracles=rows,build='synthetic-guard-fixture')
        self.right=copy.deepcopy(self.left);self.right.update(architecture='x86_64',os='Mac OS X 10.13.6')
    def test_complete_synthetic_report(self):self.assertTrue(compare_extended(self.left,self.right)['passed'])
    def test_missing_grid_or_availability(self):
        mutations=[lambda r:r['extended_oracles'].pop(),
                   lambda r:r['extended_oracles'][0]['windows'].pop(),
                   lambda r:r['extended_oracles'][0].update(lra_stable=False),
                   lambda r:r['extended_oracles'][0]['metrics']['integrated'].update(value=None),
                   lambda r:r['extended_oracles'][2]['metrics']['lra'].update(value=0),
                   lambda r:r['extended_oracles'][0]['windows'][0].__setitem__(2,float('nan')),
                   lambda r:r.update(architecture='x86_64')]
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                report=copy.deepcopy(self.left);mutation(report)
                with self.assertRaises(ValueError):validate_extended(report)
    def test_measured_difference_rejected(self):
        mutations=[lambda r:r.update(build='different'),
                   lambda r:r['extended_oracles'][0]['metrics']['integrated'].update(value=-20.011),
                   lambda r:r['extended_oracles'][0]['windows'][0].__setitem__(2,.02),
                   lambda r:r['extended_oracles'][0]['true_peaks'].__setitem__(0,.2),
                   lambda r:r['extended_oracles'][0]['sample_peak_frames'].__setitem__(0,1)]
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                report=copy.deepcopy(self.right);mutation(report)
                with self.assertRaises(ValueError):compare_extended(self.left,report)


if __name__=='__main__':unittest.main()
