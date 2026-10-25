"""Post-run text-evidence audit; not a compiled or frozen candidate input."""
from pathlib import Path
import hashlib,json,sys
R=Path(__file__).resolve().parents[4]
B=R/'build/m7-completion-2026-10-25'
E=R/'Tests/Reference/results/m7-completion-2026-10-25'
C=E
sys.path.insert(0,str(R/'Tests'))
from m7_video_acceptance import evaluate,POLICY
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
identity=json.loads((C/'candidate-identity.json').read_text())
fixture_manifest=json.loads((R/'Tests/Video/results/prores-2026-10-20/fixtures/manifest.json').read_text())
assert not [p for p,h in identity['source_files_sha256'].items() if sha(R/p)!=h]
rows=[]
for platform,arch,device,rate in [('sequoia','arm64','Mac Studio Speakers',48000),('high-sierra','x86_64','Built-in Output',44100)]:
 for phase,count in [('short',20),('endurance',3)]:
  root=E/(platform+'-video-'+phase)
  summary=json.loads((root/'summary.json').read_text())
  assert summary['accepted'] and summary['audio_configuration_unchanged'] and summary['completed_cases']==summary['expected_cases']==count
  env=json.loads((root/'identity.json').read_text())
  assert env['source_digest']==identity['source_digest'] and env['identity_sha256']==sha(C/'candidate-identity.json')
  assert env['binary_sha256']==identity['artifacts_sha256']['LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification']
  assert env['architecture']==arch and env['device']==device and env['sample_rate']==rate and env['all_routes_muted']
  runs=json.loads((root/'runs.json').read_text()); assert len(runs)==count and len(set(v['name'] for v in runs))==count
  expected_names={('prores-current-'+Path(name).stem) for name,m in fixture_manifest.items() if m.get('required')}|{'prores-current-reference','h264-cfr-source','h264-vfr-source','h264-cfr-current','h264-vfr-current'} if phase=='short' else {'prores-current-'+p+'-1800' for p in ('proxy','lt','422')}
  assert {v['name'] for v in runs}==expected_names
  assert sum(v['name'].startswith('prores-current-') for v in runs)==(16 if phase=='short' else 3)
  assert sum(v['name'].startswith('h264-') for v in runs)==(4 if phase=='short' else 0)
  for run in runs:
   path=root/(run['name']+'.json'); report=json.loads(path.read_text());assert sha(path)==run['report_sha256']
   expected=run['expected']
   assert expected['architecture']==arch and expected['device']==device and expected['sample_rate']==rate
   assert expected['measurement_build']==identity['measurement_build']
   domain='Source' if run['name'].endswith('-source') else 'Current State'
   assert expected['audition_state']==domain
   assert expected['processing_gain_db']==(0 if domain=='Source' else -3.305074118944)
   assert expected['audio_sha256']==identity['artifacts_sha256']['m7-fixtures/clock-71-'+str(rate)+'.flac']
   duration=1800 if phase=='endurance' else 100 if run['name']=='prores-current-reference' else 18
   assert expected['requested_seconds']==duration
   result=evaluate(report,expected,run['strict_exit_code'],run['fps']);assert result['accepted'] and result=={k:run[k] for k in result}
   rows.append(dict(platform=platform,phase=phase,name=run['name'],report=str(path.relative_to(R)),sha256=sha(path),
     source_seconds=report['source_seconds'],strict_harness_passed=report['passed'],accepted=True,
     skipped_pictures=report['actual_paint_skipped_frames'],skip_events=report['skip_events'],
     underrun_blocks=report['underrun_blocks'],bound_failures=report['bound_failures'],
     maximum_paint_interval_ms=report['maximum_paint_interval_ms'],
     first_last_mean_skew_difference_seconds=abs(report['last_minute_mean_skew_seconds']-report['first_minute_mean_skew_seconds']),
     peak_sampled_playback_rss_bytes=report['peak_sampled_playback_rss_bytes'],
     memory_passed=report['memory_passed'],memory_growth_assessed=report['memory_growth_assessed']))
assert len(rows)==46 and sum(v['phase']=='endurance' for v in rows)==6
out=dict(policy=POLICY,accepted=True,source_digest=identity['source_digest'],identity_sha256=sha(C/'candidate-identity.json'),
 fresh_cases=len(rows),endurance_cases=6,endurance_requested_seconds_total=10800,
 total_picture_skips=sum(v['skipped_pictures'] for v in rows),physical_71_listening='Pending; excluded from offline/automated qualification',
 all_audio_configurations_unchanged=True,rows=rows)
(E/'fresh-physical-audit.json').write_text(json.dumps(out,indent=2)+'\n')
print('Fresh physical audit accepted:',len(rows),'cases;',out['total_picture_skips'],'picture skips')
