# Six local tests with Screen Sharing connected, no SSH — 24 October 2026

**All six fixed tests passed. M7 remains frozen and unqualified overall.**

At the user's request, this repeats the six local stereo reference tests with
Screen Sharing connected to Golem's idle desktop and no SSH session. The user
confirmed a live normal-window connection and agreed to leave both Macs untouched.
No remote commands, additional test groups, source changes or builds ran.

The same Mac Studio boot session, arm64 Sequoia 15.7.9 (24G830), Apple M1 Max,
PB248 1920×1200 at 60 Hz, existing power settings, binaries, source manifest,
reference movie and stereo clock WAV match the disconnected baseline. Runs
alternated preserved M6 and frozen M7 three times in the same order, with
100 seconds of playback each through Mac Studio Speakers at 48 kHz, every route
muted. Testing ran approximately 17:27–17:38 Europe/Rome. Summed executable elapsed
time was 638.462 seconds (10 minutes 38 seconds), excluding boundary checks.

| Run | Elapsed seconds including preparation/shutdown | Result | Skipped UI frames |
| --- | ---: | --- | ---: |
| m6-stereo-1 | 106.867 | Pass | 0 |
| m7-stereo-1 | 106.345 | Pass | 0 |
| m6-stereo-2 | 106.287 | Pass | 0 |
| m7-stereo-2 | 106.401 | Pass | 0 |
| m6-stereo-3 | 106.249 | Pass | 0 |
| m7-stereo-3 | 106.313 | Pass | 0 |

Every report also records zero audio underruns, blank paints, frame-selection
bound failures and CPU-paint deadline failures. Hardware audio timestamps were
active, and existing attachment and memory limits passed. Maximum observed paint
intervals ranged from 26.175 to 36.520 ms across these six runs.

Fourteen local observations (preflight/final and before/after each run) show an
established Screen Sharing TCP connection with no SSH client or other candidate
app. These boundary observations are not a continuous network trace. Addresses
and connection details are omitted from the recorded counts. Screen Sharing was
left connected after the tests; the test apps exited and closed their audio
device before writing reports.

The fixed disconnected baseline passed 6/6 and this connected condition passed
6/6. The matched comparison is in `comparison-to-disconnected.json`. This sample
does not support attributing the failure to an idle Screen Sharing connection
alone. It does not test active remote-screen playback, SSH, concurrent remote
work, or ordinary desktop interaction, and does not establish the original cause
or a permanent fix. Earlier failures remain preserved and the remaining physical
M7 matrix/endurance qualification is unchanged.

`preflight.json`, the six invocation/outcome records, raw reports/logs and
`final-integrity.json` retain exact identities and observed conditions. All 155
M7 source inputs and 339 artifacts, both executable hashes, branch HEAD and all
existing local tags remained unchanged. The preserved baseline is unmodified.
No M7 commit, tag, push or merge was made. This addendum has its own evidence hash
manifest; earlier evidence manifests remain valid for their original file sets.
