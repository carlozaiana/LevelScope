#!/usr/bin/env python3
import unittest
from run_m15_qualification import validate,REQUIRED,validate_marker_groups,REQUIRED_MARKER_GROUPS,validate_reuse
class RequiredEvidence(unittest.TestCase):
    def rows(self):
        return [dict(name=n,exit_code=1 if n=='m14-reader-rejects-schema5' else 0,expected_exit=1 if n=='m14-reader-rejects-schema5' else 0,passed=n!='m14-reader-rejects-schema5',build='tested',reason='unsupported schema' if n=='m14-reader-rejects-schema5' else '') for n in REQUIRED]
    def test_marker_groups(self):
        validate_marker_groups(list(REQUIRED_MARKER_GROUPS))
        for group in REQUIRED_MARKER_GROUPS:
            with self.assertRaises(ValueError):validate_marker_groups([g for g in REQUIRED_MARKER_GROUPS if g!=group])
    def test_reuse_identity(self):
        a={"source_files_sha256":{"Source/a.cpp":"same","CMakeLists.txt":"same","Tests/MarkerTests.cpp":"new"},"measurement_build":"same","artifacts_sha256":{"bin/target":"same"}}
        import copy
        b=copy.deepcopy(a);b["source_files_sha256"]["Tests/MarkerTests.cpp"]="old";validate_reuse(a,b,"target")
        for key in ("source_files_sha256","artifacts_sha256"):
            bad=copy.deepcopy(b);bad[key]["Source/a.cpp" if key=="source_files_sha256" else "bin/target"]="different"
            with self.assertRaises(ValueError):validate_reuse(a,bad,"target")
    def test_complete(self):validate(self.rows(),'tested')
    def test_missing(self):
        for name in REQUIRED:
            with self.assertRaises(ValueError):validate([r for r in self.rows() if r['name']!=name],'tested')
    def test_duplicate(self):
        rows=self.rows();rows.append(rows[0])
        with self.assertRaises(ValueError):validate(rows,'tested')
    def test_wrong_build(self):
        rows=self.rows();next(r for r in rows if r['name']=='markers')['build']='stale'
        with self.assertRaises(ValueError):validate(rows,'tested')
    def test_failed_report(self):
        rows=self.rows();next(r for r in rows if r['name']=='markers')['passed']=False
        with self.assertRaises(ValueError):validate(rows,'tested')
    def test_legacy_missing_reason(self):
        rows=self.rows();next(r for r in rows if r['name']=='m14-reader-rejects-schema5')['reason']=''
        with self.assertRaises(ValueError):validate(rows,'tested')
if __name__=='__main__':unittest.main()
