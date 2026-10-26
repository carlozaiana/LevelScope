#!/usr/bin/env python3
import argparse,pathlib,json,sys,hashlib,platform,subprocess
p=argparse.ArgumentParser();p.add_argument('--candidate',type=pathlib.Path,required=True);p.add_argument('--bundles',type=pathlib.Path,required=True);p.add_argument('--measurement',type=pathlib.Path,required=True);p.add_argument('--library',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);a=p.parse_args()
sys.path.insert(0,str(a.candidate/'Tests'));from verify_m10_decoder import verify_bundle
identity=json.loads((a.candidate/'identity.json').read_text());m=json.loads(a.measurement.read_text());assert m['passed'] and len(m['cases'])==16 and m['build']==identity['measurement_build'] and m['physical_devices_opened'] is False;assert m['architecture']==platform.machine();assert all(x['passed'] for x in m['checks'])
folders=sorted(a.bundles.glob('verified-*'));assert len(folders)==16
result={'passed':True,'source_digest':identity['source_digest'],'build':m['build'],'architecture':m['architecture'],'os':subprocess.check_output(['/usr/bin/sw_vers'],universal_newlines=True),'fresh_measurement_report_sha256':hashlib.sha256(a.measurement.read_bytes()).hexdigest(),'decoder_cases':[dict(folder=f.name,**verify_bundle(f,a.library)) for f in folders],'scope':'Opposite-platform bundles freshly imported and remeasured; saved verdicts are not session authority. Independent recipient-side decoder and manifest/HTML checks.','script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()}
with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print('PASS:',m['architecture'],len(folders),'exchanged bundles')
