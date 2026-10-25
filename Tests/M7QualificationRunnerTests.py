"""Fail closed on incomplete/parity-violating M7 capture evidence."""
import pathlib
import struct
import tempfile
import unittest
import json
import hashlib
from compare_m7_reports import compare_capture, compare, validate_long_seeks


class CaptureTests(unittest.TestCase):
    def test_long_seek_matrix(self):
        rows=[dict(frame=f,gain_db=g,buffer_ms=5.,all_eight_exact=True)
              for f in (0,123456789,518352000,518399999) for g in (0.,-3.305074118944)]
        validate_long_seeks(rows)
        with self.assertRaises(ValueError): validate_long_seeks(rows[:-1])
        for field,value in [('buffer_ms',2000.0001),('buffer_ms',float('nan')),('all_eight_exact',False),('frame',-1)]:
            changed=[dict(r) for r in rows];changed[-1][field]=value
            with self.assertRaises(ValueError): validate_long_seeks(changed)

    def test_identity(self):
        data = struct.pack('<4f', 0, -.3, 1.25, -2)
        self.assertEqual(compare_capture(data, data, True), 0)

    def test_one_ulp_nonzero_only(self):
        a, b = struct.pack('<I', 0x3f800000), struct.pack('<I', 0x3f800001)
        self.assertEqual(compare_capture(a, b, False), 1)
        with self.assertRaises(ValueError): compare_capture(a, b, True)
        with self.assertRaises(ValueError): compare_capture(a, struct.pack('<I', 0x3f800002), False)

    def test_negative_and_nonfinite(self):
        self.assertEqual(compare_capture(struct.pack('<I', 0xbf800000), struct.pack('<I', 0xbf800001), False), 1)
        with self.assertRaises(ValueError): compare_capture(struct.pack('<f', float('nan')), struct.pack('<f', 0), False)

    def test_missing_and_truncated(self):
        with self.assertRaises(ValueError): compare_capture(b'\0', b'\0', False)
        with tempfile.TemporaryDirectory() as root:
            with self.assertRaises(OSError): compare(pathlib.Path(root), pathlib.Path(root))



class CompleteEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.roots=[pathlib.Path(self.temp.name)/x for x in ('arm','intel')]
        rows=[]
        for owner in (3,4):rows.append(dict(kind='muted highest peak '+str(owner),rate=48000,gain_db=3.0))
        for rate in (44100,48000,96000,192000):
            for fmt in range(9):
                for gain in (0.,-6.020599913279624,6.020599913279624,-3.305074118944,-60.,24.):
                    rows.append(dict(kind='M7 format '+str(fmt),rate=rate,gain_db=gain))
        raw=struct.pack('<8f',*[.01]*8)
        for i,row in enumerate(rows):row.update(channels=8,frames=1,file=str(i)+'.f32',sha256=hashlib.sha256(raw).hexdigest(),lra_stable=False,metrics={'integrated':dict(status=0,value=-23.)})
        for root,arch,version in zip(self.roots,('arm64','x86_64'),('15.7','10.13.6')):
            (root/'m7-samples').mkdir(parents=True)
            for row in rows:(root/'m7-samples'/row['file']).write_bytes(raw)
            report=dict(passed=True,architecture=arch,build='unit-test-only',profile_digest='ef6f5cc071a8c9d5a347da0dd61fcc46e9f24da0aed787bc9f3ea48cb3db7e66',layout_capability='pcm-7.1-bs2051-system-i-v1',mapping_version='ordered-pcm-roles-m7.1',sample_cases=rows)
            (root/'multichannel.json').write_text(json.dumps(report))
            (root/'environment.json').write_text(json.dumps(dict(product_version=dict(output=version),platform_guard_passed=True,execution_arch=arch,source_digest='unit-test-manifest')))
    def tearDown(self):self.temp.cleanup()
    def mutate(self,field,value):
        path=self.roots[1]/'multichannel.json';r=json.loads(path.read_text());r[field]=value;path.write_text(json.dumps(r))
    def test_complete_guard(self):self.assertTrue(compare(*self.roots)['passed'])
    def test_profile_and_capability(self):
        for field,value in [('profile_digest','other'),('layout_capability','count-only-8'),('mapping_version','old')]:
            path=self.roots[1]/'multichannel.json';original=path.read_text();self.mutate(field,value)
            with self.assertRaises(ValueError):compare(*self.roots)
            path.write_text(original)
    def test_missing_rate_and_channel(self):
        path=self.roots[1]/'multichannel.json';r=json.loads(path.read_text());r['sample_cases'].pop();path.write_text(json.dumps(r))
        with self.assertRaises(ValueError):compare(*self.roots)
    def test_tampered_capture(self):
        (self.roots[1]/'m7-samples/0.f32').write_bytes(b'bad')
        with self.assertRaises(ValueError):compare(*self.roots)
    def test_actual_os_guard(self):
        path=self.roots[1]/'environment.json';r=json.loads(path.read_text());r['product_version']['output']='15.7';path.write_text(json.dumps(r))
        with self.assertRaises(ValueError):compare(*self.roots)

if __name__ == '__main__':
    unittest.main()
