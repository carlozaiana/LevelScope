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
