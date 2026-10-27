#!/usr/bin/env python3
"""Close the approved phase only after base receipts and schema-2 coverage pass."""
import argparse, hashlib, json, pathlib, sys

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser()
    for name in ('package','short','long','ui','migration','output'):p.add_argument('--'+name,required=True,type=pathlib.Path)
    a=p.parse_args();sys.path.insert(0,str(a.package/'Tests'))
    from complete_m12_sequoia import complete
    result=complete(a.package,a.short,a.long,a.ui)
    identity=json.loads((a.package/'identity.json').read_text());summary=a.migration/'summary.json';m=json.loads(summary.read_text())
    assert m['passed'] is True and m['source_digest']==identity['source_digest'] and m['measurement_build']==identity['measurement_build']
    assert m['video_loaded'] is False and m['physical_devices_opened'] is False
    assert m['binary_sha256']==identity['artifacts_sha256']['bin/LevelScopeProgrammeTests']
    assert len(m['cases'])==2 and {r['channels'] for r in m['cases']}=={2,6}
    evidence={str(summary.resolve()):sha(summary)}
    for r in m['cases']:
        assert r['passed'] is True and r['originals_unchanged'] is True and r['exact_settings_roles_alignment_and_gain_preserved'] is True
        for suffix,key in (('.json','native_report_sha256'),('.log','log_sha256'),('.levelscope','output_sha256')):
            f=a.migration/(str(r['channels'])+suffix);assert sha(f)==r[key];evidence[str(f.resolve())]=sha(f)
        native=json.loads((a.migration/(str(r['channels'])+'.json')).read_text())
        assert native['passed'] and native['build']==identity['measurement_build'] and native['architecture']=='arm64' and native['reference_cases']==0 and native['physical_devices_opened'] is False
        for name,digest in r['originals_sha256'].items():assert sha(pathlib.Path(name))==digest
    result['supplemental_required_schema2_migration']=evidence
    result['completion_script_sha256']=sha(pathlib.Path(__file__))
    result['no_required_sequoia_checks_outstanding']=True
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
if __name__=='__main__':main()
