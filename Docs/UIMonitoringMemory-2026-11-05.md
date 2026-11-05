# Remembered monitoring UI slice — 5 November 2026

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
