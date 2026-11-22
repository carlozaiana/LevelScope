"""Read-only, checksum-verified session data for exact cross-platform comparison."""
import argparse
import hashlib
import json
from pathlib import Path
import struct


def sha(data): return hashlib.sha256(data).hexdigest()


class Reader:
    def __init__(self, raw): self.raw, self.pos = raw, 0
    def take(self, count):
        assert count >= 0 and self.pos + count <= len(self.raw)
        block = self.raw[self.pos:self.pos+count]; self.pos += count
        return block
    def u(self): return struct.unpack('<Q', self.take(8))[0]
    def f(self): return struct.unpack('<d', self.take(8))[0]
    def optional(self):
        present = self.u(); assert present in (0, 1)
        return self.f() if present else None


def inspect(path):
    raw = path.read_bytes(); assert raw[:8] == b'LSSESS01' and len(raw) <= 128*1024*1024
    size = struct.unpack_from('<Q', raw, 8)[0]; core_bytes = raw[80:80+size]
    assert sha(core_bytes) == raw[16:80].decode('ascii')
    core = json.loads(core_bytes)
    result = dict(path=str(path), sha256=sha(raw), document_schema=core['schema'], programme_schema=core['programme']['schema'],
                  settings=core['regionPlanning'], requirements=core['jobRequirements'], preferences=core['proposalPolicy'],
                  processing=core['current'], layout=core['programme']['layout'], channels=core['programme']['channels'], caches=[])
    for record in core['caches']:
        n, offset = int(record['size']), int(record['offset'])
        if n == 0: result['caches'].append(None); continue
        block = raw[80+size+offset:80+size+offset+n]; assert len(block) == n and sha(block) == record['sha256']
        reader = Reader(block); metadata = json.loads(reader.take(reader.u()))
        detail = dict(payload_sha256=sha(block[reader.pos:]), rate=reader.u(), frames=reader.u())
        for name in ['il', 'lra', 'max_m', 'max_s', 'sample_peak', 'true_peak']:
            detail[name] = dict(status=reader.u(), value=reader.optional())
        for name in ['il_relative_gate', 'lra_relative_gate', 'p10', 'p95']: detail[name] = reader.optional()
        for name in ['il_included', 'lra_included', 'lra_stable']: detail[name] = reader.u()
        for name in ['momentary_windows', 'short_term_windows', 'lra_windows']:
            count = reader.u(); assert count <= 1000000
            detail[name] = dict(count=count, exact_frame_energy_sha256=sha(reader.take(16*count)))
        for name in ['sample_peaks_by_channel', 'true_peaks_by_channel']:
            count = reader.u(); assert count <= 12; detail[name] = [reader.f() for _ in range(count)]
        for name in ['sample_peak_frames', 'true_peak_eighth_frames']:
            count = reader.u(); assert count <= 12
            detail[name] = [struct.unpack('<q', reader.take(8))[0] for _ in range(count)]
        assert reader.pos == len(block)
        result['caches'].append(dict(domain=record['domain'], platform=metadata['platform'], build=metadata['build'],
                                    input_sha256=metadata['content'], identity=metadata['identity'], detail=detail))
    record = core['appliedExplanation']; n, offset = int(record['size']), int(record['offset'])
    block = raw[80+size+offset:80+size+offset+n]; assert len(block) == n and sha(block) == record['sha256']
    result['history'] = dict(bytes=n, sha256=sha(block), schema=struct.unpack_from('<Q', block)[0] if block else None)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('session', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--compare', type=Path, help='Previously captured other-platform JSON')
    args = parser.parse_args(); assert not args.output.exists()
    result = inspect(args.session)
    if args.compare:
        other = json.loads(args.compare.read_text()); checks = {}
        for key in ['document_schema', 'programme_schema', 'settings', 'requirements', 'preferences', 'processing', 'layout', 'channels']:
            checks[key] = result[key] == other[key]
        checks['both_actual_platforms'] = result['caches'][0]['platform'] != other['caches'][0]['platform']
        for n, (a, b) in enumerate(zip(result['caches'], other['caches'])):
            checks['input_'+str(n)] = a is not None and b is not None and a['input_sha256'] == b['input_sha256']
            checks['build_'+str(n)] = a is not None and b is not None and a['build'] == b['build']
            checks['exact_metrics_windows_peaks_positions_'+str(n)] = a is not None and b is not None and a['detail'] == b['detail']
        result['comparison'] = dict(reference=str(args.compare), checks=checks, passed=all(checks.values()),
                                    note='Fresh proposal receipt/history identities may differ by platform/runtime; canonical processing and measurement payload must match exactly.')
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print('Captured checksum-verified session numerics'+(' — exact cross-platform PASS' if args.compare and result['comparison']['passed'] else ''))
    if args.compare and not result['comparison']['passed']: raise SystemExit(1)

if __name__ == '__main__': main()
