"""Prepare a locator-only owned session copy for cross-Mac tests; original stays intact."""
import argparse
import hashlib
import json
from pathlib import Path
import struct


def digest(data):
    return hashlib.sha256(data).hexdigest()


def relocate(original, destination, audio_locator):
    assert not destination.exists(), 'Use a new copy'
    raw = original.read_bytes()
    assert raw[:8] == b'LSSESS01'
    size = struct.unpack_from('<Q', raw, 8)[0]
    core_bytes = raw[80:80+size]
    assert digest(core_bytes) == raw[16:80].decode('ascii')
    core = json.loads(core_bytes)
    assert core['programme']['mode'] == 'interleaved' and len(core['programme']['files']) == 1
    before = core['programme']['files'][0]['locator']
    # The reader requires JUCE's canonical encoding, including its whitespace.
    # Replace exactly one ASCII string token; never reserialize the entire core
    # through Python's different JSON formatter.
    assert all(32 <= ord(c) < 127 and c not in '\\"' for c in before+audio_locator)
    old_token = ('"'+before+'"').encode('ascii')
    new_token = ('"'+audio_locator+'"').encode('ascii')
    assert core_bytes.count(old_token) == 1, 'Locator token must be unique'
    updated = core_bytes.replace(old_token, new_token, 1)
    assert len(updated) <= 4*1024*1024
    result = b'LSSESS01'+struct.pack('<Q', len(updated))+digest(updated).encode('ascii')+updated+raw[80+size:]
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open('xb') as stream: stream.write(result)
    # Reversing just the locator must recover the complete semantic core; all
    # cached measurements and historical-attachment bytes are copied verbatim.
    check = json.loads(updated)
    check['programme']['files'][0]['locator'] = before
    assert check == json.loads(core_bytes)
    assert result[80+len(updated):] == raw[80+size:]
    return dict(original=str(original), original_sha256=digest(raw), copy=str(destination), copy_sha256=digest(result),
                old_locator=before, new_locator=audio_locator, changed_semantic_fields=['programme.files[0].locator'],
                measurement_and_history_bytes_unchanged=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('original', type=Path)
    parser.add_argument('copy', type=Path)
    parser.add_argument('audio_locator')
    parser.add_argument('receipt', type=Path)
    args = parser.parse_args()
    assert not args.receipt.exists()
    result = relocate(args.original, args.copy, args.audio_locator)
    args.receipt.write_text(json.dumps(result, indent=2)+'\n')

if __name__ == '__main__':
    main()
