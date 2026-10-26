#!/usr/bin/env python3
import copy, unittest
from run_m10_qualification import validate_core,validate_long
from compare_m10_reports import compare_peaks
from verify_m10_decoder import decode_spans
class Runner(unittest.TestCase):
    def setUp(self):
        self.identity={'measurement_build':'pinned'}
        self.core=dict(passed=True,build='pinned',architecture='arm64',physical_devices_opened=False,video_cases=10,checks=[{'passed':True}],cases=[dict(domain='Exported artifact',verification_gain_db=0,dither_applications=1,evaluation={'verdict':'PASS'},sample_rate=r,per_channel_peaks=[{}]*c,input_kind=k) for r in (44100,48000,96000,192000) for c in (2,6) for k in ('interleaved','mono-group')])
        self.long=dict(passed=True,build='pinned',architecture='arm64',physical_devices_opened=False,checks=[{'passed':True}],long_export=dict(cancel_ms=10,ui_p95_ms=12,inspection_p95_ms=12,rss_growth_bytes=1024,report_bytes=10000,frames='518400000',artifact_bytes=9331200104,annotations=1000,ui_callbacks=100,rf64_audition_cases=6,rf64_reimport_audition_ms=1000,prerequisite_ms=1000,render_verify_ms=1000,publish_ms=1000))
    def test_valid(self):validate_core(self.core,self.identity,'arm64');validate_long(self.long,self.identity,'arm64')
    def test_fail_closed(self):
        for key,value in [('passed',False),('build','other'),('architecture','x86_64'),('physical_devices_opened',True),('checks',[]),('cases',[]),('video_cases',9)]:
            r=copy.deepcopy(self.core);r[key]=value
            with self.subTest(key=key),self.assertRaises((ValueError,KeyError)):validate_core(r,self.identity,'arm64')
        for key,value in [('verification_gain_db',1),('dither_applications',2),('domain','Current State'),('input_kind','unknown'),('sample_rate',32000),('evaluation',{'verdict':'FAIL'})]:
            r=copy.deepcopy(self.core);r['cases'][0][key]=value
            with self.subTest(key=key),self.assertRaises((ValueError,KeyError)):validate_core(r,self.identity,'arm64')
    def test_long_guards(self):
        for key,value in [('frames','48000'),('artifact_bytes',1000),('annotations',999),('ui_callbacks',0),('rf64_audition_cases',5),('cancel_ms',1001),('ui_p95_ms',101),('inspection_p95_ms',101),('rss_growth_bytes',256*1024*1024+1),('report_bytes',16*1024*1024+1),('render_verify_ms',10800001)]:
            r=copy.deepcopy(self.long);r['long_export'][key]=value
            with self.subTest(key=key),self.assertRaises((ValueError,KeyError)):validate_long(r,self.identity,'arm64')
    def test_peak_comparison_separates_positions_from_metric_rounding(self):
        peak=dict(role='Left',sample_peak=dict(status='available',unit='dBFS',value=-20.),
                  true_peak=dict(status='available',unit='dBTP',value=-20.),
                  sample_peak_frame='127',true_peak_eighth_frame='1017')
        other=copy.deepcopy(peak);other['true_peak']['value']=-20.000000000000007
        self.assertGreater(compare_peaks([peak],[other],True),0)
        for field,value in [('role','Right'),('sample_peak_frame','128'),('true_peak_eighth_frame','1018')]:
            wrong=copy.deepcopy(other);wrong[field]=value
            with self.subTest(field=field),self.assertRaises(AssertionError):compare_peaks([peak],[wrong],True)
        for value in (-20.02,float('nan'),float('inf'),None):
            wrong=copy.deepcopy(other);wrong['true_peak']['value']=value
            with self.subTest(value=value),self.assertRaises((AssertionError,TypeError)):compare_peaks([peak],[wrong],True)
        wrong=copy.deepcopy(other);wrong['true_peak']['status']='unavailable'
        with self.assertRaises(AssertionError):compare_peaks([peak],[wrong],True)
    def test_independent_decoder_positions_stay_inside_medium_and_long_files(self):
        for frames in (1,192000,999999,1000000,48000*1810,518400000):
            spans=decode_spans(frames)
            self.assertEqual(spans[0][0],0)
            self.assertEqual(spans[-1][0]+spans[-1][1],frames)
            for first,count in spans:self.assertTrue(0<=first<frames and 0<count<=frames-first)
        self.assertEqual(len(decode_spans(518400000)),4)
if __name__=='__main__':unittest.main()
