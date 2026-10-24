# Video paint diagnostic integration — 24 October 2026

The user authorized integrating the diagnostic tooling and evidence into main,
while keeping M7 frozen. This accepts the investigation tools for development;
it does not resolve the intermittent picture skips or change any video acceptance
limit. M7 remains frozen, unqualified and unpublished.

## Integrated source and branch boundary

Main advances from qualified M6
`b57a1ed23cbe7d74afe98ed16ec8fd615db9d171` to diagnostic implementation
[`c154b469fd8ebcea972d04e284ddde01dbc6df40`](https://github.com/carlozaiana/LevelScope/commit/c154b469fd8ebcea972d04e284ddde01dbc6df40)
by fast-forward. Fetch found no newer main changes, and no conflict resolution or
implementation adjustment was needed. The integration handover/evidence commit
that follows changes documentation and text evidence only. The published
[main branch](https://github.com/carlozaiana/LevelScope/tree/main) and this report's
Git history identify that follow-up revision.

All 121 source inputs still match diagnostic candidate-03, source digest
`0760e080d2e5daa8c89612029236308e134237a142c879003b53c7da1d3727cf`.
The production player's only addition is the read-only probe friend declaration.
The test hooks compile only into `LevelScopeVideoPaintDiagnostic`. Measurement,
audio selection, routing, gain, clock, buffering and embedded-audio policy remain
unchanged. Embedded video audio stays undecoded and unrouted.

The original `/Users/fuzz/Developer/LevelScope` checkout and its pending M7 work
remain on `codex/logic-module`. This integration does not modify that checkout or
resume M7. Subsequent app work should inspect and incorporate the latest fetched
main while retaining its own unmerged changes and the M7 freeze.

## Integration verification

A fresh Release build passed for all twelve explicitly named standalone/test
targets, with arm64 and x86_64 slices and macOS 10.13 as the Intel deployment
target. No default/all or VST3 target was built or installed. All eight native
Sequoia CTest suites passed (62.64 seconds total): video timing, offline
measurement, measured compliance, global-trim proposals, measured diagnostics,
programme sections, verified audition and standalone scaffold.

The bounded trace-buffer tests, six Python runner guards, device-closed diagnostic
startup/export and exclusive-output overwrite guard also passed. Both architecture
slices were inspected to confirm diagnostic symbols are present only in the
diagnostic executable, absent from the normal app and ordinary video harness;
their compile definitions agree. Twelve fresh executable SHA-256 identities are
recorded. The offline export intentionally has no display observer and is not a
valid playback-qualification trace. No audio device was opened for these checks.

The original device-closed checks on actual Sequoia and Intel High Sierra and the
fixed twelve muted 100-second comparisons are reused for the unchanged tested
candidate. The six raw traces remain complete, with no buffer loss. No physical
playback, listening test, new endurance run or Golem rebuild was performed for
this integration. The fresh integration binaries are separately identified;
they are not replacements for the preserved, physically tested packages and are
not claimed to have fresh High Sierra or playback qualification.

See the [integration evidence](../Tests/Video/results/paint-main-integration-2026-10-24),
[original comparison handover](VideoPaintDiagnosticComparison-2026-10-24.md) and
[initial investigation](VideoPaintTimingInvestigation-2026-10-24.md). Historical
reports retain their original paths, hashes, stage/status and acceptance results.

## Open frame-skip finding and synchronization limit

The latest twelve comparisons passed with zero missed painted frames, underruns,
blank paints, frame-bound failures or existing CPU-paint deadline failures. Earlier
intermittent failures remain recorded and unresolved. Passing comparisons do not
establish a scheduling fix, a permanent absence of failures or M7 acceptance.

Earlier evidence is consistent with intermediate pictures being selected but not
painted, followed by pictures selected against the continuing authoritative audio
clock. It did not demonstrate accumulating internal synchronization drift. This
distinguishes a brief picture stutter from an accumulating timing offset; it does
not prove every picture was physically presented on time. Optical/acoustic
synchronization remains unmeasured, and these muted runs provide no listening
judgement.

Other development may continue while this investigation stays open. M7's pending
video qualification and disposition of the intermittent failure remain required;
no acceptance criterion has been relaxed by this merge.

## Using the diagnostic for a recurrence

The normal standalone app and ordinary qualification harness do not automatically
record these traces. The optional diagnostic executable must be selected explicitly,
and a valid trace must be saved. It records frame selection, repaint requests,
paint callbacks, an independent display observer, main-thread vblank observations,
queue/transport identity and window state. These can help locate a delay but do
not directly measure compositor presentation. Trace overflow or export failure
invalidates the diagnostic evidence.

The current diagnostic is bounded to 100 seconds with all routes muted. Any next
physical run needs fresh user coordination and a fixed run budget. A useful next
investigation would reproduce one documented earlier failing condition with the
preserved control and diagnostic, then stop at the declared limit regardless of
outcome. A further passing sample alone does not close the failure. No further
playback or feature work is authorized by this integration record.

## Locations and preservation

Source remains `/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation`.
This chat returns to `codex/video-implementation` after main integration. Build
only named standalone/test targets in a new directory for every new build; never
use default/all or rebuild over a preserved candidate.

Fresh integration build and full build logs:
`/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/build/paint-main-integration-20261024-01`.
Its CMake cache is specific to that source/build location.

Preserved packages:

- M6: `/Users/fuzz/Developer/LevelScope/build/m6-qualification-2026-10-23/candidate-01` (308 manifest artifacts).
- M7: `/Users/fuzz/Developer/LevelScope/build/m7-qualification-2026-10-24/candidate-02` (339 manifest artifacts).
- Diagnostic: `/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/build/paint-diagnostic-20261024-01/candidate-03` (143 manifest artifacts).
- Golem diagnostic: `/Users/fuzz/LevelScope-VideoPaintDiagnostic-2026-10-24/candidate-03`, with its earlier comparison results alongside it.
- Golem M7: `/Users/fuzz/LevelScope-M7-qualification-2026-10-24/candidate-02`.

All ten pre-existing checkpoint tag objects and target commits are preserved.
Their exact identities and package-integrity checks are recorded in the integration
evidence. Golem's preserved copies, installed apps/plugins, Pro40 and system audio
configuration were not modified. VST3 parity, embedded soundtrack playback and
additional features remain separate work.
