#!/usr/bin/env python3
"""Six fixed local muted reproductions with user-connected Screen Sharing and idle SSH; no remote commands, build or source edits."""
import datetime, hashlib, json, os, pathlib, platform, subprocess, sys, time
ROOT=pathlib.Path('/Users/fuzz/Developer/LevelScope')
OUT=ROOT/'build/m7-qualification-2026-10-24/screen-sharing-idle-ssh-local-six-01'
M7=ROOT/'build/m7-qualification-2026-10-24/candidate-02'
M6=ROOT/'build/m6-qualification-2026-10-23/candidate-01'
SUFFIX='LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification'
BINS={'m6':M6/SUFFIX,'m7':M7/SUFFIX}
EXPECTED={'m6':'c8efc1a6d0b760206c8f99105e553542a4d32fdb9722dc3e19b2172f60bb36e9','m7':'7b0113714b90fa582dd6ae725652630e7bd4f601fecc608fcb4dd15acd0d512e'}
AUDIO=pathlib.Path('/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/build/prores-20261020/fixtures/selected-clock-48000.wav')
VIDEO=ROOT/'Tests/51_SyncTest.mov'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def command(args):
 r=subprocess.run(args,cwd=ROOT,text=True,capture_output=True)
 if r.returncode:raise RuntimeError('Preflight read failed: '+repr(args)+' '+r.stderr)
 return r.stdout.strip()
def write(name,value):
 with (OUT/name).open('x') as f:json.dump(value,f,indent=2);f.write('\n')
def processes():
 records=[]
 for line in command(['/bin/ps','-axo','pid=,comm=']).splitlines():
  pid,comm=line.strip().split(None,1)
  if pathlib.Path(comm).name in ('ssh','Screen Sharing','LevelScope Video Qualification','LevelScope Standalone'):
   records.append({'pid':int(pid),'executable':comm})
 return records
def connection_state():
 records=processes()
 viewers=[r for r in records if pathlib.Path(r['executable']).name=='Screen Sharing']
 clients=[r for r in records if pathlib.Path(r['executable']).name=='ssh']
 forbidden=[r for r in records if pathlib.Path(r['executable']).name not in ('Screen Sharing','ssh')]
 counts=[];viewer_hosts=set();ssh_hosts=set()
 for kind,group in [('screen_sharing',viewers),('ssh',clients)]:
  for process in group:
   result=subprocess.run(['/usr/sbin/lsof','-nP','-a','-p',str(process['pid']),'-iTCP','-sTCP:ESTABLISHED','-Fn'],text=True,capture_output=True)
   connections=[line[1:].strip() for line in result.stdout.splitlines() if line.startswith('n') and '->' in line]
   hosts={line.split('->',1)[1].rsplit(':',1)[0] for line in connections}
   (viewer_hosts if kind=='screen_sharing' else ssh_hosts).update(hosts)
   # Compare destination hosts in memory; retain only counts/booleans, never addresses.
   counts.append({'kind':kind,'pid':process['pid'],'established_tcp_connections':len(connections),'lsof_exit_code':result.returncode})
 shared=bool(viewer_hosts & ssh_hosts)
 matching=bool(ssh_hosts) and ssh_hosts.issubset(viewer_hosts)
 return {'recorded_utc':now(),'viewer_processes':viewers,'ssh_client_processes':clients,'forbidden_candidate_processes':forbidden,'connection_counts':counts,
         'shared_destination_host':shared,'all_ssh_destinations_match_viewer_hosts':matching,
         'both_connected_to_same_host':bool(viewers) and bool(clients) and not forbidden and shared and matching,
         'idle_ssh_basis':'User confirmed idle remote shell; no remote command sent by runner. Boundary socket observations do not measure remote activity.'}
OUT.mkdir(exist_ok=False)
manifest=json.loads((M7/'identity.json').read_text())
source_errors=[p for p,h in manifest['source_files_sha256'].items() if sha(ROOT/p)!=h]
artifact_errors=[p for p,h in manifest['artifacts_sha256'].items() if sha(M7/p)!=h]
pre={'recorded_utc':now(),'purpose':'Six fixed local muted tests with Screen Sharing and idle SSH connected; three M6 and three M7 stereo alternating in the same order as the disconnected baseline; no extra retries.',
 'user_report':'Same post-reboot Mac Studio session. User confirmed Ready after connecting idle SSH to Golem at its shell prompt, with Screen Sharing connected to idle desktop, and agreed to leave both Macs untouched.',
 'boot_time':command(['/usr/sbin/sysctl','kern.boottime']),
 'os':command(['/usr/bin/sw_vers']), 'architecture':platform.machine(),
 'cpu':command(['/usr/sbin/sysctl','-n','machdep.cpu.brand_string']),
 'power_settings_read_only':command(['/usr/bin/pmset','-g','custom']),
 'display_configuration_read_only':json.loads(command(['/usr/sbin/system_profiler','SPDisplaysDataType','-json'])),
 'relevant_processes_before':processes(),'connections_before':connection_state(),
 'head':command(['git','rev-parse','HEAD']),
 'tag_objects_and_targets':command(['git','for-each-ref','--format=%(refname:short) %(objectname) %(*objectname)','refs/tags']),
 'source_digest':manifest['source_digest'],'source_inputs_checked':len(manifest['source_files_sha256']),'source_errors':source_errors,
 'm7_artifacts_checked':len(manifest['artifacts_sha256']),'artifact_errors':artifact_errors,
 'binaries_sha256':{k:sha(p) for k,p in BINS.items()},'audio_sha256':sha(AUDIO),'video_sha256':sha(VIDEO),
 'runner_sha256':sha(pathlib.Path(__file__)),
 'scope':'Internal audio-clock/UI-paint reproduction only; no full M7 qualification or physical 7.1 listening claim.',
 'sequence':[f'{model}-stereo-{i}' for i in (1,2,3) for model in ('m6','m7')]}
write('preflight.json',pre)
assert not source_errors and not artifact_errors
assert pre['binaries_sha256']==EXPECTED
assert pre['audio_sha256']=='4709852fdf45b43636171c07868ad6d0c7199cec2f9a5ae18c2ce33094d1a528'
assert pre['video_sha256']=='0870d5ab93d1492f23524ed1fff9319af1a6a9cf71882543f2190929afeb8e24'
assert pre['architecture']=='arm64'
assert pre['connections_before']['both_connected_to_same_host'],pre['connections_before']
print('PREFLIGHT PASS: binaries, 155 source inputs, 339 artifacts; Screen Sharing and SSH have established TCP connections to the same host; no candidate app process. Starting fixed six-run sequence.',flush=True)
runs=[];fatal=[]
for name in pre['sequence']:
 model=name.split('-')[0]
 boundary_before=connection_state()
 if not boundary_before['both_connected_to_same_host']:
  fatal.append(name+': Screen Sharing/SSH connection condition changed or candidate process present before run; sequence stopped')
  write(name+'-connection-abort.json',boundary_before)
  break
 report=OUT/(name+'.json')
 argv=[str(BINS[model]),'--allow-audio','--muted','--device','Mac Studio Speakers','--seconds','100','--audio',str(AUDIO),'--video',str(VIDEO),'--report',str(report),'--audition-state','source']
 write(name+'-invocation.json',{'name':name,'started_utc':now(),'binary_sha256':EXPECTED[model],'argv':argv,'connection_before':boundary_before})
 print('START',name,now(),flush=True)
 start=time.monotonic()
 with (OUT/(name+'.log')).open('x') as log:
  proc=subprocess.Popen(argv,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
  try:rc=proc.wait(timeout=210)
  except subprocess.TimeoutExpired:
   proc.terminate()
   try:proc.wait(timeout=5)
   except subprocess.TimeoutExpired:proc.kill();proc.wait()
   rc=proc.returncode;fatal.append(name+': timeout; process terminated; sequence stopped')
 row={'name':name,'ended_utc':now(),'exit_code':rc,'wall_seconds':time.monotonic()-start,'binary_sha256':sha(BINS[model])}
 if report.exists():
  data=json.loads(report.read_text())
  row.update({k:data.get(k) for k in ('passed','error','source_seconds','requested_seconds','sample_rate','all_routes_muted','actual_paint_skipped_frames','maximum_paint_interval_ms','bound_failures','underrun_blocks','blank_paints','paint_deadline_failures','maximum_frame_pts_skew_seconds','memory_passed','attachment_within_budget','skip_events')})
  row['report_sha256']=sha(report)
  if not (data.get('all_routes_muted') is True and data.get('sample_rate')==48000 and data.get('device')=='Mac Studio Speakers' and data.get('audio_sha256')==pre['audio_sha256'] and data.get('video_sha256')==pre['video_sha256'] and data.get('source_seconds',0)>=100 and not data.get('error')):
   fatal.append(name+': incomplete run or configuration mismatch; sequence stopped')
 else:fatal.append(name+': no report; sequence stopped')
 row['connection_after']=connection_state()
 if not row['connection_after']['both_connected_to_same_host']:fatal.append(name+': Screen Sharing/SSH connection condition changed or candidate process present after run; sequence stopped')
 if row['binary_sha256']!=EXPECTED[model]:fatal.append(name+': binary changed; sequence stopped')
 runs.append(row);write(name+'-outcome.json',row)
 print('FINISH',name,json.dumps({k:row.get(k) for k in ('exit_code','wall_seconds','passed','actual_paint_skipped_frames','bound_failures','underrun_blocks')}),flush=True)
 if fatal:break
final={'recorded_utc':now(),'source_errors':[p for p,h in manifest['source_files_sha256'].items() if sha(ROOT/p)!=h],
 'artifact_errors':[p for p,h in manifest['artifacts_sha256'].items() if sha(M7/p)!=h],
 'binaries_sha256':{k:sha(p) for k,p in BINS.items()},'relevant_processes_after':processes(),'connections_after':connection_state(),
 'head_unchanged':command(['git','rev-parse','HEAD'])==pre['head'],
 'tags_unchanged':command(['git','for-each-ref','--format=%(refname:short) %(objectname) %(*objectname)','refs/tags'])==pre['tag_objects_and_targets'],
 'device_closure_basis':'Each completed harness disconnects before writing its report; every child process was waited for.'}
final['passed']=not final['source_errors'] and not final['artifact_errors'] and final['binaries_sha256']==EXPECTED and final['connections_after']['both_connected_to_same_host'] and final['head_unchanged'] and final['tags_unchanged']
write('final-integrity.json',final)
summary={'recorded_utc':now(),'scope':pre['scope'],'planned_runs':6,'completed_runs':len(runs),'sequence_complete':len(runs)==6 and not fatal,
 'all_six_passed':len(runs)==6 and not fatal and all(r.get('passed') and r['exit_code']==0 for r in runs),
 'fatal_errors':fatal,'m7_qualified':False,'m7_code_changed':False,'baseline':'screen-sharing-local-six-01','connection_observation':'Established Screen Sharing and SSH TCP connections to matching destination hosts, with no candidate apps, checked before/after each completed run; no continuous network trace. User confirms idle Golem shell and desktop; runner sends no remote commands.','runs':runs,'integrity_passed':final['passed']}
write('summary.json',summary)
print('COMPLETE',json.dumps({k:v for k,v in summary.items() if k!='runs'}),flush=True)
sys.exit(0 if summary['all_six_passed'] and final['passed'] else 1)
