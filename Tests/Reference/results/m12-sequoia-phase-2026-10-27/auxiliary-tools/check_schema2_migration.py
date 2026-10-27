#!/usr/bin/env python3
"""Bounded supplemental schema-2 migration using the frozen native reader."""
import argparse, hashlib, json, pathlib, struct, subprocess

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def core(p):
    b=p.read_bytes();n=struct.unpack_from('<Q',b,8)[0]
    assert hashlib.sha256(b[80:80+n]).hexdigest().encode()==b[16:80]
    return json.loads(b[80:80+n])
def without_locators(v):
    # Identical-content relinking refreshes filesystem locators/timestamps. The
    # Tucan archive has second-resolution mtimes; neither field is audio identity.
    if isinstance(v,dict):return {k:without_locators(x) for k,x in v.items() if k not in ('locator','modified')}
    if isinstance(v,list):return [without_locators(x) for x in v]
    return v

def main():
    p=argparse.ArgumentParser()
    for name in ('package','inputs','output'):p.add_argument('--'+name,required=True,type=pathlib.Path)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    identity=json.loads((a.package/'identity.json').read_text());binary=a.package/'bin/LevelScopeProgrammeTests'
    assert sha(binary)==identity['artifacts_sha256']['bin/LevelScopeProgrammeTests']
    rows=[]
    for channels in (2,6):
        original=a.inputs/str(channels)/'programme.levelscope';old=core(original)
        assert old['schema']==2 and old['programme']['mode']=='mono_files' and not old['video']['locator']
        members=[original.parent/pathlib.Path(x['locator']).name for x in old['programme']['files']]
        before={str(x):sha(x) for x in [original]+members}
        dest=a.output/(str(channels)+'.levelscope');report=a.output/(str(channels)+'.json');log=a.output/(str(channels)+'.log')
        argv=['/usr/bin/arch','-arm64',str(binary),'--only','--prepare-session',str(original),'--destination',str(dest),'--output',str(report)]
        with log.open('w') as f:subprocess.run(argv,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=60)
        result=json.loads(report.read_text());assert result['passed'] and result['build']==identity['measurement_build'] and result['architecture']=='arm64' and result['reference_cases']==0 and result['physical_devices_opened'] is False
        new=core(dest);assert new['schema']==4 and new['programme']['schema']==4
        for key in ('documentId','profile','sections','video','audition'):assert new[key]==old[key],key
        for key in ('initialized','gain','planSchema','engine'):assert new['current'][key]==old['current'][key],key
        assert new['current']['sectionGains'] is None
        for key in old['programme']:
            if key!='schema':assert without_locators(new['programme'][key])==without_locators(old['programme'][key]),key
        assert before=={str(x):sha(x) for x in [original]+members}
        rows.append(dict(channels=channels,passed=True,argv=argv,originals_sha256=before,output_sha256=sha(dest),native_report_sha256=sha(report),log_sha256=sha(log),exact_settings_roles_alignment_and_gain_preserved=True,originals_unchanged=True,filesystem_metadata_policy='Verified identical-content relinking may refresh locator and mtime; all audio identities, confirmation, settings and hashes must match.'))
    summary=dict(passed=True,source_digest=identity['source_digest'],measurement_build=identity['measurement_build'],binary_sha256=sha(binary),cases=rows,video_loaded=False,physical_devices_opened=False,scope='Supplemental actual Sequoia schema-2 stereo/5.1 migration; same frozen candidate; no rebuild or long rerun')
    (a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
if __name__=='__main__':main()
