#!/usr/bin/env python3
"""Acceptance is true only after every separate M10 gate has recorded success."""
import pathlib,json,hashlib,datetime
r=pathlib.Path(__file__).resolve().parent;i=json.loads((r/'candidate-03/identity.json').read_text());gates=[]
def gate(name,files,field='passed',value=True):
 evidence=[];errors=[]
 for name_ in files:
  p=r/name_
  if not p.is_file():errors.append('Missing '+name_);continue
  d=json.loads(p.read_text());ok=d.get(field)==value
  if not ok:errors.append('Not passed: '+name_)
  evidence.append({'path':name_,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'required_field':field,'required_value':value})
 gates.append({'name':name,'passed':not errors,'errors':errors,'evidence':evidence})
gate('Both actual-platform core and regressions',[p+'-candidate-03/summary.json' for p in ('sequoia','high-sierra')],'automated_passed')
gate('Both real three-hour RF64 and independent decoder',[p+'-candidate-03-long/summary.json' for p in ('sequoia','high-sierra')],'automated_passed')
gate('Cross-platform measurements and captured samples',['m10-cross-platform-candidate-03.json','m10-long-cross-platform-candidate-03.json','m1-m5-cross-platform-candidate-03.json','m6-cross-platform-candidate-03.json','m7-m9-core-captures-candidate-03.json','m8-cross-platform-candidate-03.json','m9-extended-cross-platform-candidate-03.json'])
gate('Fresh bidirectional export exchange and independent recipient decoders',['exchange-high-sierra-to-sequoia-candidate-03-acceptance.json','exchange-sequoia-to-high-sierra-candidate-03-acceptance.json'])
gate('Silent actual-platform UI',['ui-'+p+'-candidate-03/observations.json' for p in ('sequoia','high-sierra')],'status','passed')
gate('Coordinated physical video under accepted skip policy',[p+'-physical-'+phase+'-candidate-03/summary.json' for p in ('sequoia','high-sierra') for phase in ('short','endurance')],'accepted')
gate('Independent aggregate of original physical case receipts',['physical-video-summary-candidate-03.json'])
gate('Long-file per-channel values and exact locations',['m10-long-detail-comparison-candidate-03.json'])
gate('Frozen packages and current source integrity',['local-integrity-candidate-03.json','high-sierra-integrity-candidate-03.json','m9-package-preservation-candidate-03.json'])
preserve=json.loads((r/'repository-preservation-final-candidate-03.json').read_text());ok=all(preserve[k] is True for k in ('source_manifest_matches','existing_tags_unchanged','worktrees_unchanged','head_unchanged','unrelated_unchanged'));gates.append({'name':'Existing tags, branches, worktrees and unrelated work','passed':ok,'evidence':[{'path':'repository-preservation-final-candidate-03.json','sha256':hashlib.sha256((r/'repository-preservation-final-candidate-03.json').read_bytes()).hexdigest()}]})
errors=[]
for platform,arch,prefix in [('sequoia','arm64','15.'),('high-sierra','x86_64','10.13.')]:
 for suffix in ('','-long'):
  p=r/(platform+'-candidate-03'+suffix)/'environment.json'
  if not p.exists():errors.append('Missing environment: '+str(p));continue
  e=json.loads(p.read_text())
  if e['source_digest']!=i['source_digest'] or e['build']!=i['measurement_build'] or e['architecture']!=arch or not e['os'].startswith(prefix):errors.append('Wrong actual-platform identity: '+str(p))
gates.append({'name':'Exact actual-platform source/build identity','passed':not errors,'errors':errors})
qualified=all(g['passed'] for g in gates);d={'milestone':'M10 verified PCM export and reporting','qualified':qualified,'candidate':'candidate-03','source_digest':i['source_digest'],'base_commit':i['base_commit'],'measurement_build':i['measurement_build'],'identity_sha256':hashlib.sha256((r/'candidate-03/identity.json').read_bytes()).hexdigest(),'app_sha256':i['artifacts_sha256']['LevelScope Standalone.app/Contents/MacOS/LevelScope Standalone'],'gates':gates,'assessed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'publication_and_merge':'Not performed; awaiting later user readiness approval','accepted_limitations':['Occasional picture skips accepted only with continuous audio and authoritative-clock synchronization; counts/timing retained.','Full physical multichannel listening remains separately pending.','Limited PCM24 stereo/conventional5.1 RIFF/RF64; no complete broadcaster/BWF/BW64/metadata certification.','No resampling, automatic correction, batch/network delivery, video/lossy export, or VST3 parity.']}
(r/'acceptance-candidate-03.json').write_text(json.dumps(d,indent=2)+'\n')
print('QUALIFIED' if qualified else 'INCOMPLETE');print('\n'.join(g['name'] for g in gates if not g['passed']))
