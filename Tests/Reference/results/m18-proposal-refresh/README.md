# M18 proposal UI hang: development regression

A five-second sample of the user's running candidate-04 app on macOS 15.7.9
showed the process at 99.7% CPU. The main thread had 448 sampled stacks in the
event loop servicing timer callbacks. Proposal panel refresh, gain-lane sync
and session dirty-state checks repeatedly reconstructed the Source and Current
snapshot digests. Those digests include every measured loudness window. The raw
process sample stays local; it contains host-specific information.

The proposal path now reuses the original SHA-256 digest of each published,
immutable measurement. Each thread has two weakly owned entries. Replaced
snapshots calculate a fresh digest, mutable/raw inputs remain uncached, and
every request still checks its current context, processing, target, preferences,
annotation geometry and applicable source file stamps. The fix grants no
measurement or Apply authority.

The regression fixture uses 36,000 momentary, 36,000 short-term and 3,600 LRA
window entries in each domain. It compares both identities to fresh uncached
digests, measures repeated panel/lane/actionable-status refresh, replaces one
window energy by the adjacent representable value, verifies stale preview and
Apply refusal, changes expected session identity, and checks that old window
arrays can be released. This synthetic population exercises UI cost without
running a full-duration media job.

The user narrowed the trigger to Propose followed by **Proposal limits…**,
then a second click when the editor appeared behind the timeline. Each click
previously replaced the entire native window. The action now raises and focuses
the same open editor, restores it if minimized, and preserves unconfirmed text.
After closing the window, reopening it starts from the current session's limits.
Native two-window checks exercise the first open, bringing the timeline forward,
the second click, typing an unconfirmed duration and reopening after close.

The separate build directory is `build/m18-refresh-fix` so compilation and
validation leave the user's running executable and session untouched. Named
standalone, section proposal workflow/preferences, global proposal and measured compliance
targets only. Physical hardware/video, actual High Sierra, cross-Mac and
full-duration/RSS qualification stay deferred.

The original app became responsive after the proposal was discarded. Its CPU
fell to 4%, the window again answered accessibility requests, and the proposal
panel reported no proposal to discard. The original instance was not terminated.

Final [candidate 02](candidate-02/summary.json) passed **1,136 checks**:
395 proposal workflow, 256 preferences, 303 global proposal and 182 measured
compliance. The large-window refresh took 1.02 ms p95 on the
development Mac. A single uncached measurement pair took
59.34 ms on the same fixture. The native peers verified
foreground focus and draft preservation. Opening the actual corrected app's
limits window also showed its interactive controls in the foreground; the window
was then closed without changing preferences. Candidate 01 is an earlier passing
digest-only pilot before the window-opening fix.

Development app: `/Users/fuzz/.codex/worktrees/a411/LevelScope/build/m18-refresh-fix/LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app`.
Build `8e566ab112b5632fb395c46b5a9d7e49dd7b685f97768d888f622fc644ec72ac`; the app and three test executables are universal arm64/x86_64
with Intel minimum 10.13. This does not qualify actual High Sierra.
