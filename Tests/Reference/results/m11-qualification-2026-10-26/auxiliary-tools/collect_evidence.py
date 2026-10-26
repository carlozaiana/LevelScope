#!/usr/bin/env python3
"""Collect text-only M11 receipts; never touch frozen packages or raw media."""
import argparse
import hashlib
import json
import pathlib
import shutil


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--work', type=pathlib.Path, required=True)
    parser.add_argument('--output', type=pathlib.Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    records = []

    def copy(source, relative):
        target = args.output / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and digest(target) != digest(source):
            raise ValueError('Refusing to replace historical evidence: ' + str(relative))
        if not target.exists():
            shutil.copy2(source, target)
        records.append(dict(path=str(relative), sha256=digest(source), bytes=source.stat().st_size))

    for candidate in ('candidate-01', 'candidate-02'):
        copy(args.work / candidate / 'identity.json', pathlib.Path(candidate + '-identity.json'))
    disposition = args.work / 'candidate-01-disposition.json'
    if disposition.exists():
        copy(disposition, pathlib.Path(disposition.name))
    results = args.work / 'candidate-02-results'
    for platform in ('sequoia', 'high-sierra'):
        for phase in ('short', 'long', 'exchange'):
            source = results / platform / phase
            if not source.exists():
                continue
            for file in sorted(source.iterdir()):
                if file.is_file() and file.suffix in ('.json', '.log'):
                    copy(file, pathlib.Path(platform) / phase / file.name)
            if phase == 'long':
                for name in ('report.json', 'report.html', 'manifest.json'):
                    file = source / 'long-files/export' / name
                    if file.exists():
                        copy(file, pathlib.Path(platform) / phase / 'long-files/export' / name)
        for file in sorted((results / platform).glob('ui*.json')):
            copy(file, pathlib.Path(platform) / file.name)
    for name in ('short-cross-platform.json', 'complete-cross-platform.json'):
        if (results / name).exists():
            copy(results / name, pathlib.Path(name))
    # Each invocation is separately identified; incomplete collection is not acceptance.
    collection = dict(qualification_claim=False, text_files=records)
    name = 'collection-' + hashlib.sha256(json.dumps(collection, sort_keys=True).encode()).hexdigest()[:12] + '.json'
    (args.output / name).write_text(json.dumps(collection, indent=2) + '\n')
    print(len(records), 'text receipts copied; no media or binaries')


if __name__ == '__main__':
    main()
