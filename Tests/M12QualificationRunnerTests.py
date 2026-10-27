#!/usr/bin/env python3
import copy, hashlib, json, struct, unittest
from package_m12_qualification import audio_only_session
from run_m12_qualification import REQUIRED, UI_CASES, MANDATORY, validate_runs, validate_ui, validate_report, reuse_compatible

class Inventory(unittest.TestCase):
    def test_audio_only_fixture_preserves_canonical_bytes(self):
        core=b'{"label": "unchanged", "video": {"locator": "/old/movie.mov", "sha256": "abc123", "stream": 1, "offset": "-1001/30000"}, "cache": 3}'
        original=b'HEADER01'+struct.pack('<Q',len(core))+hashlib.sha256(core).hexdigest().encode()+core+b'original-cache-bytes'
        derived=audio_only_session(original);n=struct.unpack_from('<Q',derived,8)[0];expected=core.replace(b'"/old/movie.mov"',b'""').replace(b'"abc123"',b'""').replace(b'"stream": 1',b'"stream": 0')
        self.assertEqual(derived[80:80+n],expected);self.assertEqual(derived[80+n:],b'original-cache-bytes')
        self.assertEqual(derived[16:80],hashlib.sha256(expected).hexdigest().encode())
        self.assertEqual(json.loads(expected)['video']['offset'],'-1001/30000')
        self.assertEqual(audio_only_session(derived),derived)
        with self.assertRaises(ValueError):audio_only_session(original.replace(b'"video"',b'"other"'))
    def test_reuse_requires_exact_production_and_binary(self):
        old=dict(source_files_sha256={'CMakeLists.txt':'a','Source/core.h':'b','Tests/test.cpp':'c'},measurement_build='d',artifacts_sha256={'bin/Test':'e'})
        new=copy.deepcopy(old);new['source_files_sha256']['Tests/test.cpp']='updated test';self.assertTrue(reuse_compatible(old,new,'Test'))
        for kind in ('source','build','binary'):
            bad=copy.deepcopy(new)
            if kind=='source':bad['source_files_sha256']['Source/core.h']='changed'
            if kind=='build':bad['measurement_build']='changed'
            if kind=='binary':bad['artifacts_sha256']['bin/Test']='changed'
            self.assertFalse(reuse_compatible(old,bad,'Test'))
        self.assertFalse(reuse_compatible(old,new,'Missing'))
    def test_required_runs(self):
        for phase,names in REQUIRED.items():
            rows=[dict(name=x,exit_code=0) for x in names];validate_runs(rows,phase)
            for bad in (rows[:-1],rows+[rows[0]],[dict(x,exit_code=1) for x in rows]):
                with self.assertRaises(ValueError):validate_runs(bad,phase)
            with self.assertRaises(ValueError):validate_runs(rows+[dict(name='video',exit_code=0)],phase)
    def test_candidate_and_matrix(self):
        identity=dict(measurement_build='a');r=dict(passed=True,build='a',architecture='arm64',phase='sequoia-audio',physical_devices_opened=False,video_loaded=False,
          checks=[dict(name=x,passed=True) for x in MANDATORY],cases=[dict(rate=rate,channels=ch,kind='M12 matrix input-'+str(mode)) for rate in (44100,48000,96000,192000) for ch in (2,6) for mode in range(4)])
        validate_report(r,identity,'sequoia-audio')
        for key,value in (('passed',False),('build','wrong'),('architecture','x86_64'),('video_loaded',True),('physical_devices_opened',True),('cases',r['cases'][:-1]),('checks',r['checks'][:-1])):
            with self.assertRaises(ValueError):validate_report(dict(r,**{key:value}),identity,'sequoia-audio')
    def test_ui(self):
        identity=dict(source_digest='a',artifacts_sha256={'LevelScope Standalone.app/Contents/MacOS/LevelScope Standalone':'b'});r=dict(passed=True,source_digest='a',app_sha256='b',architecture='arm64',device_closed=True,routes_muted=True,video_loaded=False,cases=[dict(id=x,passed=True,evidence='receipt') for x in UI_CASES]);validate_ui(r,identity)
        for key,value in (('cases',r['cases'][:-1]),('video_loaded',True),('source_digest','wrong'),('app_sha256','wrong'),('device_closed',False)):
            with self.assertRaises(ValueError):validate_ui(dict(r,**{key:value}),identity)
    def test_long_budgets(self):
        v=dict(frames='518400000',rate=48000,channels=6,intervals=1000,capture_cases=8,total_ms=100,cancel_ms=10,ui_p95_ms=15,additional_peak_rss_bytes=1,session_bytes=1,report_bytes=1,stages_ms=dict(save_ms=1,parse_ms=1))
        r=dict(passed=True,build='a',architecture='arm64',phase='sequoia-long',physical_devices_opened=False,video_loaded=False,checks=[dict(name='long',passed=True)],long_workflow=v);identity=dict(measurement_build='a');validate_report(r,identity,'sequoia-long')
        for key,value in (('intervals',999),('ui_p95_ms',101),('cancel_ms',float('nan')),('total_ms',10800001)):
            bad=copy.deepcopy(r);bad['long_workflow'][key]=value
            with self.assertRaises(ValueError):validate_report(bad,identity,'sequoia-long')
if __name__=='__main__':unittest.main()
