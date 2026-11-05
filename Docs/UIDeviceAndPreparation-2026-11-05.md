# UI continuation after M19 — 5 November 2026

Authorised scope: two focused implementation slices on
`codex/ui-device-preparation-m19`, from published main
`5fc8c6ee95e705077cf996cecbc6543dbc20d257`. The M18 branch, M19/Intel checkpoints,
frozen packages and other working copies remain preserved. Publication and main
integration are for subsequent review. Documents 63 and 69 define the roadmap
and independent work; the strategy/audio-quality chat remains design-only.

## Slice 1: remembered output connection

Keep local device preferences outside `.levelscope` documents. Remember only a
successfully opened output, with the backend and CoreAudio UID where safely
available. Key the preference file by the Mac's machine identity so copied user
preferences do not select another Mac's interface. Restore its visible selection;
the new automatic reconnect checkbox defaults to off. No fallback to the default
device or another same-name device is permitted. Ambiguous identification fails
with an explanation and an explicit chooser/retry path.

Explicit selection/Connect is connection intent. If the programme rate is not yet
known, retain that pending intent until a validated programme supplies it. With
automatic reconnection off, simply restoring memory does not connect. Opt-in
reconnection attempts only at explicit lifecycle events, not on every UI timer
tick; failed attempts remain visible without retry loops. A changed programme rate
closes an incompatible existing connection. Automatic reopening requires the
preference; otherwise Connect is explicit. Opening always leaves transport stopped,
routes muted and monitor gain at the existing safe level. Speaker-role confirmation
and programme/measurement identity are independent of device preferences.

Disconnect suppresses automatic reconnection until this application session ends,
including New/Open and preference toggles. Explicit selection/Connect remains
available. This lifetime and the pending-rate state are explained in the UI.
Missing, ambiguous, busy/unavailable and incompatible outputs do not block offline
analysis. Existing output overload and transport safeguards remain unchanged.

Focused validation: injected catalog/device tests for separate Mac stores,
successful-only persistence, renamed/duplicate/missing identities, preference on/off,
pending rate, rate changes, rejected/busy opens, disconnect suppression, no callback
start, routes/monitor safety and session reopen. Reuse M19 measurement/proposal/
history/export and Intel parity evidence because their implementation is unchanged.
Run the affected audition and timeline interaction suites; build a preserved universal
standalone/test candidate with Intel minimum 10.13 and run the focused device checks
on actual High Sierra. No physical opening or audible test without user coordination.

## Slice 2: initial import and analysis feedback

Expose current source validation, preparation/verification and measurement stages
prominently in the main window, with meaningful per-stage progress and cancellation.
Keep worker progress bounded and thread-safe. Preserve existing job/revision acceptance,
preview coverage/gate labels and accepted-result authority. Cancellation during
preparation must retire the job without a restart loop or partial measurement authority.

Focused validation will cover stage transitions, progress mapping, preparation
cancellation, late completions and incomplete/final preview distinctions, with actual
High Sierra checks for affected paths. Reuse unchanged numerical and full-duration
qualification. First-playback feedback, session-wide verified-source reuse, project-owned
audio, early waveform generation during validation and the inherited completion pause
remain separate. No new DSP, proposal method, measurement, history or export behaviour.

Each slice receives its own exact commit, preserved test build, focused results and
manual checklist before publication review.

## Slice 1 review receipt

Product/test commit: `55e5cd1603ed4b6366b244fed12895d48e90099b`.
Measurement build: `935fcf7c38de1ac1bc2240d07e7281d24c6f17f179c6e75c2b0a20615ee42dcc`.
Preserved development app (do not rebuild over it):

`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/ui-device-settings-2026-11-05/candidate-01/LevelScope Standalone.app`

Identical Golem candidate:
`/Users/fuzz/LevelScope-UI-device-2026-11-05/candidate-01/LevelScope Standalone.app`.
Both packages contain the source snapshot, exact binary/source SHA-256 manifest,
compiler configuration and universal test executables. Intel deployment minimum is
10.13; the same frozen executables ran on actual 10.13.6, not only under Rosetta.

| Focused group | Sequoia 15.7.9 arm64 | High Sierra 10.13.6 Intel |
| --- | ---: | ---: |
| Device persistence / reconnection / safety | 170 passed | 171 passed |
| Injected audition | 4,059 passed | 4,059 passed |
| Timeline / roles / interaction | 889 passed | 889 passed |
| Session playback lifecycle / reopen | 1,820 passed | 1,820 passed |

Device counts differ only because read-only native catalog enumeration has one
check per visible output. All opens/playback in these suites use injected devices.
No physical output was opened, no audible test performed, and no system default,
Pro40, installation or driver changed. Hashes were checked before and after testing.
The first sandboxed local audition attempt could not create its existing JUCE-owned
temporary directory; that failed attempt is retained separately. The unchanged
candidate passed with the required temporary-directory access. There was no product
failure or fix between the two attempts.

Evidence: `Tests/Reference/results/ui-device-settings-2026-11-05/`, including the
candidate identity, both host summaries, full focused checks and rendered Audio setup.
Unchanged M19/Intel measurement, proposal, history, export, long-duration and video
qualification is reused; this receipt does not claim new physical routing coverage.

Manual checklist (user-coordinated hardware use):

1. Open Audio setup. Automatic reconnect starts off. Select an interface, import a
   mix and confirm roles. Once validation supplies its rate, the selected output
   connects, transport stays stopped, monitor is −30 dB and all routes are muted.
2. Restart with the preference off: the last successful interface is selected but
   disconnected. Connect explicitly. Enable the preference and restart/import:
   compatible remembered output reconnects, still stopped and muted.
3. Disconnect, then open another project or import another rate. Automatic connection
   remains paused until restart; explicit Connect is available. With automatic
   reconnect off, a changed incompatible rate disconnects and waits for Connect.
4. With the remembered interface absent or an incompatible rate, inspect the reason;
   confirm there is no substituted output and offline analysis remains available.
   Reopen a saved project and verify its speaker roles are retained independently.

If CoreAudio exposes ambiguous duplicate/suffixed names that cannot be joined safely
to a unique device UID through JUCE's public interface, selection is refused with a
request to give the devices unique names in Audio MIDI Setup and rescan. The app
never guesses which same-named device to open.

## Slice 2 implementation notes

The fixed banner is directly below the main toolbar and above Original / Current /
Target. It names source validation, Original or Current preparation/verification,
measurement, and finalisation. Decoded sample checks, FLAC checksum, source identity,
private-copy progress and private-copy identity checks publish their real per-pass
coverage. Percentages are explicitly for the current file/pass (or programme during
measurement), not an estimated percentage of all remaining work. Metadata checks and
finalisation use an indeterminate indicator instead of a fabricated remaining time.
File groups and composition imports expose their existing worker's feedback and
cancellation to this same banner; their separate role/alignment confirmation remains.

`JobFeedback` holds one coherent atomic stage/fraction snapshot. Worker-local children
forward optional feedback to their parent, independently of the existing numeric job
progress. Production retains no stage history or PCM. The banner polls on the existing
UI timer; no feedback is sent to the audio callback. A test receiver observes/cancels
at actual stage boundaries without time-sensitive sleeps or modifying the source.
Feedback is enabled for import and analysis, not first-playback preparation.

Cancel immediately shows that cancellation is pending, then the existing job ticket /
revision checks discard late results. Preparation cancellation removes the owned copy;
measurement cancellation retires incomplete waveform and loudness previews. The
existing once-per-source automatic-analysis policy prevents a cancellation/retry loop.
Manual Measure Source / Re-measure Current State remains explicit. Existing completed
measurements retain their original authority while a repeat measurement is pending.

The waveform budget, per-channel extrema/positions, measured loudness windows, gate
classification and completion acceptance are unchanged. The source verification passes,
measurement engine, proposal strategy, applied explanation, session format and export
logic are unchanged. Normal build/cache identity checks still apply. The inherited
completion pause remains visible as finalisation; it has not been fixed in this slice.

## Slice 2 review receipt

Product/test commit: `0a942b542a4c51824da19940f9500680e5a23c80`.
Measurement build: `2c5187d1c6cd122903fc419347308421c54899a46a5e2175b8610012eb74e991`.
Preserved development app, including both slices:

`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/ui-preparation-feedback-2026-11-05/candidate-01/LevelScope Standalone.app`

Identical Golem candidate:
`/Users/fuzz/LevelScope-UI-preparation-2026-11-05/candidate-01/LevelScope Standalone.app`.
The source/artifact manifest and compiler configuration are beside each app. App and
focused executables are universal arm64/x86_64, with Intel minimum macOS 10.13.

| Focused group | Sequoia 15.7.9 arm64 | High Sierra 10.13.6 Intel |
| --- | ---: | ---: |
| Preparation feedback / cancellation / exact parity / main workflow | 946 passed | 946 passed |
| Progressive waveform / measurement observer parity | 1,305 passed | 1,305 passed |
| Timeline / progressive loudness / interaction | 889 passed | 889 passed |

These are executed results from the frozen candidate on both actual operating
systems. Preparation checks cancel inside copy, copy-hash, decoded validation,
source-hash, FLAC checksum, measurement and finalisation. They verify owned-copy
cleanup, untouched originals, no cancelled-result authority, no automatic retry loop,
manual retry, stale completion rejection and initial Current's identical-Source origin.
Stereo and 5.1 feedback-on/off comparisons exactly match all metrics and availability,
per-channel sample/true peaks and positions, gates, percentiles, populations and every
M/S/LRA window. The existing waveform suite also covers Current processing/transitions,
assemblies/compositions, supported rates/channels, memory and display reduction.

UI images verify the stage banner layout with staged progress; live worker preview
coverage/authority is tested separately. Visual inspection on both Macs confirmed
readable stage, per-pass percentage and cancellation. The finalisation indicator
remains indeterminate. Candidate artifact hashes match before and after the runs.
Evidence: `Tests/Reference/results/ui-preparation-feedback-2026-11-05/`.

Manual checklist:

1. Open Mix and choose an interleaved file or channel-file group. Watch the banner
   immediately below the toolbar identify source validation and the current pass.
   Channel roles and separate-file alignment still require their existing confirmation.
2. Confirm roles and analyse. Watch preparation/copy/identity stages change to measurement.
   Waveform and M/S curves appear progressively; Original and initial Current stay pending
   until the accepted completed Source result exists. The banner identifies finalisation
   separately from completion.
3. On a sufficiently long file, Cancel during preparation. Confirm it visibly enters
   Cancelling, returns without a completed measurement and does not restart itself.
   Use Measure Source to retry; repeat cancellation during progressive measurement.
4. Let a retry finish. Confirm both Original and untouched Current become ready and
   Current retains its inspectable identical-Source measurement origin. Inspection and
   timeline navigation do not start audition. Playback preparation feedback is unchanged.

Review state: both slices are local on `codex/ui-device-preparation-m19`; neither is
published or integrated into main. Prior M18/M19/Intel checkpoints and the original UI
branch remain intact. No physical output opening/listening, new video qualification or
long-duration rerun was performed. Reuse unchanged qualification; retain the separate
source-reuse, project-owned audio, early-waveform, first-playback and completion-pause
work items, and the later stereo monitoring design. Proposal strategy/audio quality
remains the other chat's design-only assignment.
