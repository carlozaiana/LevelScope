# Remembered monitoring UI slice — 5 November 2026

**Publication update, 6 November 2026:** all three convenience slices are now
published and integrated into main at checkpoint `ui-convenience-2026-11-06`.
See the [publication handover](UIConveniencePublication-2026-11-06.md) for exact
identities and preservation checks. The review-time status below is historical;
the frozen packages and qualification evidence remain unchanged.

User-authorised continuation on `codex/ui-device-preparation-m19`, after the two
completed device/preparation slices. Their commits and frozen packages remain
preserved. Publication and main integration remain for review.

Local monitoring memory is separate from the connection-preference file and portable
sessions. Its namespace uses the same machine identity as connection preferences.
Within it, backend + stable interface UID (or exact name when UID is unavailable)
select the interface. Stereo and 5.1 routes are independent, keyed by speaker role;
the existing supported side-surround aliases map to Ls/Rs. Confirmation of actual
programme roles is still required. The complete ordered interface output-name list
must match, including its count and unused outputs. No routing is inferred from a
channel count, another interface, or another programme layout. The app cannot detect
physical cable/speaker changes that do not change the interface's exposed metadata.

A changed/unavailable configuration keeps routes muted without overwriting stored
assignments. Audio setup shows the reason and the saved/available lists (full text
also in the tooltip). Edit the outputs or explicitly choose Remember routing after
review; this also saves a deliberately all-muted configuration. Valid routing edits
save immediately. Existing duplicate-output, role, output-range and no-downmix
safeguards remain. Other supported programme layouts retain manual routing.

Monitor gain remembers only explicit user notifications, separately per Mac/interface.
Internal restoration uses no-notification updates and does not write preferences.
The existing -60 through 0 dB range, half-dB UI steps, gain smoothing and overload
protection remain. No valid saved value means the initial -30 dB default. Neither
routing nor volume changes affect programme dirtiness, processing, measurement
identity or export. File storage is bounded and atomically replaced.

Read-only output discovery may construct a JUCE device metadata wrapper; it never
calls open/start or changes hardware settings. Matching monitoring can therefore be
shown while disconnected, with reconnection off or explicitly suppressed. Connection
intent remains owned by the prior slice. Restoration occurs stopped and never invokes
Play. Routine refresh only checks context/output changes, not preference writes.

Focused checks: expanded device suite for new/reopened/restarted sessions, both layouts,
reordered 5.1 channels, assigned/muted/all-muted roles, independent interfaces/Macs,
changed/missing output configurations, volume defaults/restoration/range, reconnection
on/off and Disconnect suppression. Retain affected injected audition/smoothing and
session-playback tests plus the preparation workflow integration checks. Run the same
frozen universal candidate on Sequoia and actual High Sierra 10.13.6. No physical
opening or audible tests without user coordination. Reuse unchanged M19/Intel,
long-duration, video, proposal, measurement/history/export qualification.

## Review receipt

Exact product/test commit: `96159ebd79f762bea147b41ad696fd80f449dfc9`.
Measurement build: `4a44fc52dd01be5947cc5627041391dac42ac5e2ee76044c60dc2d8f0faa3b9a`.
Standalone executable SHA-256: `62c76abdbaa2ec0ed7a299d56eec4374de34910b2b7eeca9b8a4232bfde0e55a`.

Fresh preserved development app, containing all three convenience slices:

`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/ui-monitoring-memory-2026-11-05/candidate-01/LevelScope Standalone.app`

Identical preserved Golem app:

`/Users/fuzz/LevelScope-UI-monitoring-2026-11-05/candidate-01/LevelScope Standalone.app`

Both packages include exact source/artifact manifests, configuration and focused test
executables. App and test binaries are universal arm64/x86_64; Intel minimum is 10.13.
The frozen binaries ran on actual High Sierra 10.13.6 and Sequoia 15.7.9. Artifact
hashes match before/after each run. Both earlier device/preparation packages were
also checked against their original manifests and remain unchanged.

| Focused suite | Sequoia arm64 | High Sierra Intel |
| --- | ---: | ---: |
| Device intent / routing / volume / native read-only metadata | 445 passed | 445 passed |
| Injected audition / smoothing / overload safeguards | 4,059 passed | 4,059 passed |
| Initial import / preparation / main-window workflow | 946 passed | 946 passed |
| Session reopen / playback lifecycle | 1,820 passed | 1,820 passed |

The routing suite covers actual save/read/verify/install, application-object restart,
New and fresh imports, confirmation gates, alternating stereo/5.1, reordered six-channel
files and side-surround aliases, assigned/muted/all-muted configurations, two interfaces
and two Mac namespaces, output-name/count changes, missing devices, live configuration
changes, invalid saved volume, range endpoints, preference on/off and explicit Disconnect.
It also covers a pending user route choice during routine refresh and deferred JUCE
notifications after internal restoration. Volume notifications are captured synchronously;
internal no-notification changes never become saved user edits.

An explicit Play through an injected callback verifies the actual software output:
reordered speaker destinations, muted LFE and the restored -24 dB monitor gain after
the existing startup ramp. Restoration itself starts no callback. Physical speakers,
interfaces, system defaults, Pro40 and drivers were not opened or modified. Native
metadata discovery only reads device/output properties; it does not play or connect.

Curated full checks, summaries, logs and rendered routing UI on both Macs are under
`Tests/Reference/results/ui-monitoring-memory-2026-11-05/`. During development, adding
multiple test engines exceeded the test process's stack limit (`___chkstk_darwin`).
The engines in that harness were moved to heap ownership; the correction changed no
production code. Both frozen-candidate runs subsequently passed. No frozen failed
candidate was replaced or overwritten.

Manual checklist (physical connection/listening remains user-coordinated):

1. In Audio setup, configure stereo destinations, deliberately mute a role, and set a
   monitor level. New, Open and application restart should restore those settings on
   the same Mac/interface; with reconnect off, the interface stays disconnected.
2. Configure 5.1 separately, then import a differently ordered 5.1 file and confirm its
   roles. Check L/R/C/LFE/Ls/Rs against their intended outputs. Alternate with stereo;
   each layout keeps its own assignments and mutes. Both use that interface's volume.
3. Mute every role and use Remember routing if they were already initially muted.
   Reopen/restart and verify no channel is assigned automatically. Change interfaces
   and return: each interface retains its own volume and routing; an unsaved one uses
   muted routing and the initial -30 dB default.
4. Use Disconnect with reconnect enabled, then New/Open: monitoring restores but the
   interface does not reconnect until explicit selection/Connect or an app restart.
   Playback stays stopped until Play. With missing/changed outputs, Scan devices and
   inspect the reason and old/new lists; review assignments before saving the new map.

No publication or main integration has occurred. Measurement, proposal, history and
export implementation is unchanged. Reuse their M19/Intel and full-duration/video
qualification. The earlier progressive-waveform/loudness evidence remains applicable;
this slice changes local monitoring only. The separate first-playback, source-reuse,
project-owned audio, early-waveform, completion-pause and stereo-fold-down work remains
outside this slice.
