"""Review preserved M20 evidence; numerical/UI completion never implies listening approval."""
import argparse
import hashlib
import json
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    return json.loads(path.read_text())


def main():
    p = argparse.ArgumentParser()
    p.add_argument('root', type=Path)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    root, foreign = args.root, args.root / 'golem-receipts-04'
    result = dict(non_audible_complete=False, m20_combined_qualified=False,
                  listening_pending=True, publication_authorized=False, checks=[], hosts={}, comparisons=[])

    def check(ok, name):
        result['checks'].append(dict(name=name, passed=bool(ok)))
        if not ok:
            raise AssertionError(name)

    def receipt(path, build):
        data = read(path)
        check(data.get('passed') and data.get('build') == build and bool(data.get('checks'))
              and all(v['passed'] for v in data['checks']), 'Complete passed receipt: ' + str(path))
        return data

    try:
        identities = {}
        for n in range(1, 5):
            folder = root / ('candidate-%02d' % n)
            info = read(folder / 'identity.json')
            check(sha(folder / 'identity.json') == (folder / 'identity.sha256').read_text().split()[0], 'Frozen identity %d' % n)
            check(all(sha(folder / k) == v for k, v in info['artifacts_sha256'].items()), 'Every frozen artifact %d' % n)
            check(sha(folder / 'identity.json') == sha(foreign / folder.name / 'identity.json'), 'Same package on actual Macs %d' % n)
            identities[n] = info
        baseline, current = identities[1], identities[4]
        source_changes = [k for k in sorted(set(baseline['source_sha256']) | set(current['source_sha256']))
                          if baseline['source_sha256'].get(k) != current['source_sha256'].get(k)]
        allowed_product = {'Source/Standalone/RegionProposalReviewPanel.cpp', 'Source/Standalone/RegionProposalReviewPanel.h',
                           'Source/Standalone/StandaloneMainComponent.cpp', 'Source/Standalone/StandaloneProposalPreferences.cpp',
                           'Source/Standalone/StandaloneProposalPreferences.h'}
        check({k for k in source_changes if k.startswith('Source/')} == allowed_product,
              'Reuse is limited to explicit native review retention/paging changes; all other product inputs identical')
        result['source_changes_after_numerical_candidate'] = source_changes
        result['product_commit'] = current['source_commit']
        result['measurement_build'] = current['measurement_build']
        result['source_manifest_sha256'] = current['source_manifest_sha256']
        groups = {'control-variants', 'core', 'stereo', 'surround', 'heterogeneous', 'failures', 'population', 'audio',
                  'gating', 'cinema', 'devices', 'preparation', 'timeline', 'requirements', 'section-smoke',
                  'legacy-core', 'legacy-stereo', 'legacy-failures', 'legacy-gating'}
        by_host = {}
        for host, numeric, native, exchange, architecture, os_prefix, limit in [
            ('sequoia', root / 'sequoia-01', root / 'sequoia-native-04', root / 'exchange-on-sequoia-04', 'arm64', '15.', 10),
            ('high-sierra', foreign / 'high-sierra-remainder-01', foreign / 'high-sierra-native-04', foreign / 'exchange-on-high-sierra-04', 'x86_64', '10.13.', 20)]:
            summary = read(numeric / 'summary.json')
            check(summary['passed'] and summary['architecture'] == architecture and summary['os'].startswith(os_prefix), 'Actual platform and complete numerical run: ' + host)
            selected = [g for g in summary['groups'] if g['name'] != 'ui']
            check({g['name'] for g in selected} == groups and len(selected) == len(groups), 'All 19 explicit affected groups: ' + host)
            records = {}
            for g in selected:
                path = numeric / g['receipt']
                check(g['passed'] and g['returncode'] == 0 and sha(path) == g['receipt_sha256'], 'Original numerical receipt integrity: ' + host + '/' + g['name'])
                records[g['name']] = receipt(path, baseline['measurement_build'])
            native_run = read(native / 'summary.json')
            check(native_run['passed'] and native_run['architecture'] == architecture and native_run['os'].startswith(os_prefix), 'Actual current native platform: ' + host)
            native_data = receipt(native / 'ui-files/result.json', current['measurement_build'])
            check(not native_data['physical_devices_opened'], 'Native checks used no physical outputs: ' + host)
            heartbeats = [x for x in native_data['cases'] if 'native_heartbeat_p95_ms' in x]
            cancel = next(x['cancel_ms'] for x in native_data['cases'] if 'cancel_ms' in x)
            check(len(heartbeats) == 2 and all(x['samples'] > 20 and x['native_heartbeat_p95_ms'] <= 100 for x in heartbeats) and cancel <= 1000,
                  'Unchanged native response and cancellation criteria: ' + host)
            population = records['population']['cases'][0]
            check(population['64_proxy_seconds'] <= limit and population['analysis_sidecar_bytes'] <= 64 * 1024 * 1024,
                  'Actual cached-planning target and bounded population storage: ' + host)
            exchanged = receipt(exchange / 'result.json', current['measurement_build'])
            check(len(exchanged['cases']) == 2 and {v['layout'] for v in exchanged['cases']} == {'stereo', 'surround'}
                  and all(v['origin'] != v['receiver'] for v in exchanged['cases']), 'Actual stereo/5.1 cross-platform session exchange: ' + host)
            result['hosts'][host] = dict(os=summary['os'], architecture=architecture,
                                        numerical_checks=sum(g['checks'] for g in selected), native_checks=len(native_data['checks']),
                                        exchange_checks=len(exchanged['checks']), heartbeat=heartbeats, cancel_ms=cancel, population=population)
            by_host[host] = records
        for group in ('stereo', 'surround', 'heterogeneous', 'gating', 'cinema'):
            left, right = (by_host[h][group]['cases'] for h in ('sequoia', 'high-sierra'))
            remove_timing = lambda rows: [{k: v for k, v in row.items() if k != 'planning_seconds'} for row in rows]
            check(remove_timing(left) == remove_timing(right), 'Exact reported decisions, scalar results and activity/resource counts: ' + group)
            result['comparisons'].append(dict(group=group, exact_reported_values=True))
        left, right = (by_host[h]['control-variants']['cases'] for h in ('sequoia', 'high-sierra'))
        check(len(left) == len(right) == 15, 'All 15 pre-holdout control variants on both Macs')
        for a, b in zip(left, right):
            check(all(a[k] == b[k] for k in ('variant', 'outcome', 'passes', 'proxy_count', 'retained_attempt'))
                  and [x['lra'] for x in a['actual_attempts']] == [x['lra'] for x in b['actual_attempts']],
                  'Same variant decision and actual LRA: ' + a['variant'])
        held = receipt(root / 'holdouts-01/result.json', current['measurement_build'])
        check(len(held['cases']) == 6, 'Both fixed methods completed on all three independent inputs')
        for info in read(root / 'holdout-inputs-before.json'):
            check(sha(Path(info['path'])) == info['sha256'], 'Original holdout byte-exact: ' + info['path'])
        result['independent_holdout_checks'] = len(held['checks'])
        result['non_audible_complete'] = True
        result['remaining'] = 'User will audition the two independently verified music comparison sessions. No listening/quality conclusion, M20 publication or main integration.'
    except Exception as error:
        result['failure'] = str(error)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('non_audible_complete', 'm20_combined_qualified', 'listening_pending')}, indent=2))
    raise SystemExit(0 if result['non_audible_complete'] else 1)


if __name__ == '__main__':
    main()
