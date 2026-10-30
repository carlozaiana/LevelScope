#!/usr/bin/env python3
import copy,unittest
from run_m16_qualification import REQUIRED,GROUPS,LEGACY_BUILD,validate
class Evidence(unittest.TestCase):
    def rows(self):return [dict(name=n,exit_code=1 if n=='m15-reader-rejects-schema6' else 0,expected_exit=1 if n=='m15-reader-rejects-schema6' else 0,passed=n!='m15-reader-rejects-schema6',reason='foreign session has exact applied plan and provenance' if n=='m15-reader-rejects-schema6' else '',build=LEGACY_BUILD if n=='m15-reader-rejects-schema6' else 'frozen',groups=list(GROUPS),legacy_input_schema=6,legacy_input_sha256='a'*64) for n in REQUIRED]
    def test_complete(self):validate(self.rows(),'frozen')
    def test_missing(self):
        for n in REQUIRED:
            with self.assertRaises(ValueError):validate([r for r in self.rows() if r['name']!=n],'frozen')
    def test_duplicate(self):
        a=self.rows();a.append(a[0])
        with self.assertRaises(ValueError):validate(a,'frozen')
    def test_groups(self):
        for group in GROUPS:
            a=self.rows();next(r for r in a if r['name']=='preferences')['groups'].remove(group)
            with self.assertRaises(ValueError):validate(a,'frozen')
    def test_bad_result(self):
        for key,value in [('build','other'),('passed',False),('exit_code',1)]:
            a=self.rows();next(r for r in a if r['name']=='preferences')[key]=value
            with self.assertRaises(ValueError):validate(a,'frozen')
    def test_legacy(self):
        for key,value in [('reason','missing file'),('build','other'),('legacy_input_schema',5),('legacy_input_sha256','')]:
            a=self.rows();next(r for r in a if r['name']=='m15-reader-rejects-schema6')[key]=value
            with self.assertRaises(ValueError):validate(a,'frozen')
if __name__=='__main__':unittest.main()
