#!/usr/bin/env python3
import copy,unittest
from run_m17_qualification import REQUIRED,GROUPS,POLICY_GROUPS,LEGACY_BUILD,validate
class Evidence(unittest.TestCase):
    def rows(self):return [dict(name=n,exit_code=1 if n=='m16-reader-rejects-schema7' else 0,expected_exit=1 if n=='m16-reader-rejects-schema7' else 0,passed=n!='m16-reader-rejects-schema7',reason='foreign session has exact applied plan and provenance' if n=='m16-reader-rejects-schema7' else '',build=LEGACY_BUILD if n=='m16-reader-rejects-schema7' else 'frozen',groups=list(POLICY_GROUPS if n=='preferences' else GROUPS),legacy_input_schema=7,legacy_input_sha256='a'*64) for n in REQUIRED]
    def test_complete(self):validate(self.rows(),'frozen')
    def test_missing(self):
        for n in REQUIRED:
            with self.assertRaises(ValueError):validate([r for r in self.rows() if r['name']!=n],'frozen')
    def test_duplicate(self):
        a=self.rows();a.append(a[0])
        with self.assertRaises(ValueError):validate(a,'frozen')
    def test_groups(self):
        for group in GROUPS:
            a=self.rows();next(r for r in a if r['name']=='requirements')['groups'].remove(group)
            with self.assertRaises(ValueError):validate(a,'frozen')
    def test_bad_result(self):
        for key,value in [('build','other'),('passed',False),('exit_code',1)]:
            a=self.rows();next(r for r in a if r['name']=='requirements')[key]=value
            with self.assertRaises(ValueError):validate(a,'frozen')
    def test_expected_exit_cannot_override_policy(self):
        a=self.rows();r=next(r for r in a if r['name']=='requirements');r['exit_code']=r['expected_exit']=1
        with self.assertRaises(ValueError):validate(a,'frozen')
    def test_policy_groups(self):
        for group in POLICY_GROUPS:
            a=self.rows();next(r for r in a if r['name']=='preferences')['groups'].remove(group)
            with self.assertRaises(ValueError):validate(a,'frozen')
    def test_legacy(self):
        for key,value in [('reason','missing file'),('build','other'),('legacy_input_schema',5),('legacy_input_sha256','')]:
            a=self.rows();next(r for r in a if r['name']=='m16-reader-rejects-schema7')[key]=value
            with self.assertRaises(ValueError):validate(a,'frozen')
if __name__=='__main__':unittest.main()
