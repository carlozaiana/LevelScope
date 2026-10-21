# Qualified ProRes integration into main — 21 October 2026

Accepted integration revision: [`af7ed751ab869f7829c65e7736939e45425c0d54`](https://github.com/carlozaiana/LevelScope/commit/af7ed751ab869f7829c65e7736939e45425c0d54).
The following handover/evidence commit changes no production or test source.
Continue logic-module work from the latest fetched [published main](https://github.com/carlozaiana/LevelScope/tree/main),
which must contain this accepted revision, qualified M5 and the fixed ProRes tag.
Bring main into the existing logic development branch while preserving unmerged work.
This task does not switch or alter that branch or implement another feature.

## Publication verification

Main and the video branch were first published at
`83bc6adc64fc495f5a2ada3e2974e6fe508d583f`. GitHub
[run 39560275560](https://github.com/carlozaiana/LevelScope/actions/runs/39560275560)
passes the universal standalone/test build and all seven CTest suites, then uploads
only `LevelScope-standalone-macOS`. CI artifacts are separate build checks; they do
not replace the frozen endurance-qualified package. The final publication-audit
follow-up changes documentation/evidence only and skips a redundant CI rebuild;
production, test code, CMake and the verified workflow remain byte-identical.

## Integration and scope

Fetched main was `6d304dbecd8fa47f56b0260105aa9d07f96ba1bc`, already an ancestor of
video branch `2e7000af21ce7f6b57b9eca470eafc71ba045077`; there were no newer main
changes or conflicts. Integration uses a checked fast-forward, preserving every
commit. The only integration configuration change restricts GitHub CI to nine named
standalone/test targets and adds seven hardware-independent CTest checks. The old
default build included VST3 and its install rule; CI now packages only standalone.
No plugin target, installation, audio device or system setting is changed.

The qualified source remains `6662091f1dfaff948a4a3eb3f000f47e37cfdc21`, with
production last changed at `3920d08f849d2420509e1d789c5b8208e50d79af`.
All 135 source-manifest hashes match; the complete CMake/Source/test-code diff from
that source is empty. Candidate-04 identity SHA-256 is
`45219de9006b845230c5e1296859bff2c483beb49c4b827b8c33c552f8381285`.
All 286 frozen artifacts match on both Macs before and after these regressions.

Support remains local MOV with ProRes 422 Proxy, LT and standard 422, progressive,
square pixels, full aperture and explicit Rec.709 SDR, decoded to the 8-bit preview.
Required 1080p rates are 24000/1001, 24, 25, 30000/1001 and 30 fps on both actual Macs.
The selected audio clock, exact PTS, offsets, timecode, seeking and bounded queue
remain authoritative. Embedded movie audio stays **undecoded and unrouted**.

## Integration regressions and retained qualification

| Check | Result |
| --- | --- |
| Fresh universal Release build, macOS minimum 10.13 | Nine named standalone/test targets built; no VST3 product |
| Fresh-build CTest on Sequoia | All seven registered suites pass |
| Frozen candidate native ProRes on Sequoia and High Sierra | 52 cases per Mac pass, including actual 51_SyncTest.mov, picture conversion, frame identity/seeks, metadata rejection, lifecycle and nonzero selected-audio isolation |
| Frozen candidate M1–M5, scaffold, H.264 and runner regressions | Both actual platforms pass, including ITU/EBU reference vectors |
| Ordinary cross-platform numerical/performance comparison | Pass, identical candidate/package identity; no reuse-identity exception |
| Existing matrix, UI, large-file loading/cancellation/memory and six 30-minute runs | Reused unchanged candidate-04 evidence; no affected implementation changes |

The fresh build validates CI targets and integration, and has separate executable
hashes. It does not replace the frozen candidate or inherit its endurance claim.
Fresh tests are inaudible and open no physical audio devices. Existing matrix and
endurance results remain internal selected-audio-clock/paint measurements under the
approved memory method. No previous failure is regraded and no new listening or
physical synchronization measurement is claimed.

Raw new reports, identities, comparison and preservation checks are in
[the integration evidence](../../Tests/Video/results/prores-integration-2026-10-21/README.md).
The complete accepted evidence and earlier failures remain in
[the ProRes qualification report](26_ProResQualification.md).

## Checkpoints, source locations and subsequent work

The existing annotated [video-prores-checkpoint-2026-10-21](https://github.com/carlozaiana/LevelScope/tree/video-prores-checkpoint-2026-10-21)
stays at `dcc34aabba6b5b3ca4805e2965d21ffebf2e4d0a`, tag object
`d62c1159eeb2ffbb8b7ffc7c31014b842642d195`. All nine pre-existing tag objects and
peeled commit targets are preserved, including measurement, scaffold, High Sierra,
M2–M5 and both video checkpoints. No tags are added, replaced or force-pushed.

Video source: `/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation`,
`codex/video-implementation`. Original logic checkout:
`/Users/fuzz/Developer/LevelScope`, `codex/logic-module`, remains at
`ae88e187d288c53dbca1287e8341b74b9b5b2619`; tracked, staged and untracked files are
preserved. The final byte inventory detected a change only to the original root
`.DS_Store` Finder metadata; its before/after hashes are recorded in the publication
audit, and its current contents were left untouched. The strict inventory result
is retained as a mismatch, not overwritten with a blanket pass. Registered
worktrees stay at their existing locations.

Preserved tested local package:
`/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/build/prores-20261020/candidate-04`.
Its sibling fixtures, old candidates, raw evidence and historical caches remain.
Golem's tested package remains
`/Users/fuzz/LevelScope-ProRes-qualification-2026-10-20/candidate-04`.
Fresh integration build/logs/results:
`/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/build/prores-integration-20261021`.
Fresh Golem results: `/Users/fuzz/LevelScope-ProRes-integration-2026-10-21`.
Use a new build directory under the video working copy for later builds; do not
reuse or rebuild over preserved qualification directories or historical caches.

The 33 GB fixture requires 184.637 s on Sequoia and 277.960 s on High Sierra before
the first decoded frame is available, predominantly for whole-file hashing. This is
initial attachment time, not an on-screen measurement or a playback sync offset.
It satisfies the declared size-based test budget but remains a loading usability
limitation. See the clarification in document 26; no loading feature is added here.

The next logic milestone should read the measurement contract, M5 handover and
`ProposalEngineSpec_v2.md` and preserve measurement/source identity, video timing
and both macOS versions. Optional embedded soundtrack playback and VST3 parity stay
separate. Physical optical/acoustic synchronization remains unmeasured; 4K/higher
rates, VLC/LibVLC, other codecs, HQ/4444/RAW, HDR/log conversion, interlacing, MXF and
export/transcoding remain outside this accepted milestone.
