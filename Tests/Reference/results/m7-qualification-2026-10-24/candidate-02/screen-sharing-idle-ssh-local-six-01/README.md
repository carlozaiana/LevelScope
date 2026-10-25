# Six local tests with Screen Sharing and idle SSH — 24 October 2026

**All six fixed tests passed. M7 remains frozen and unqualified overall.**

The user connected SSH to Golem, left it at its idle shell prompt, kept Screen
Sharing connected to the idle desktop, confirmed Ready and agreed to leave both
Macs untouched. The runner issued no remote commands. No additional qualification
groups, builds, source changes or installations ran.

This used the same Mac Studio boot session, arm64 Sequoia 15.7.9 (24G830), Apple
M1 Max, PB248 1920×1200 at 60 Hz, power configuration, binaries, source manifest,
stereo clock WAV, reference movie and acceptance limits as the two earlier sets.
Three preserved M6 and three frozen M7 runs alternated in the same order, each
with 100 seconds of Source stereo playback through Mac Studio Speakers at 48 kHz,
all routes muted. Testing ran approximately 17:49–18:00 Europe/Rome. Summed
executable elapsed time was 638.336 seconds (10 minutes 38 seconds), excluding
local boundary checks.

| Run | Elapsed seconds including preparation/shutdown | Result | Skipped UI frames |
| --- | ---: | --- | ---: |
| m6-stereo-1 | 106.792 | Pass | 0 |
| m7-stereo-1 | 106.315 | Pass | 0 |
| m6-stereo-2 | 106.316 | Pass | 0 |
| m7-stereo-2 | 106.334 | Pass | 0 |
| m6-stereo-3 | 106.274 | Pass | 0 |
| m7-stereo-3 | 106.306 | Pass | 0 |

All six reports also record zero audio underruns, blank paints, frame-selection
bound failures and CPU-paint deadline failures. Hardware audio timestamps were
active; attachment and memory limits passed. Maximum observed paint intervals
ranged from 27.514 to 38.280 ms across the six runs.

Fourteen local observations (preflight/final and before/after each run) found
established Screen Sharing and SSH TCP connections with matching destination
hosts and no other candidate app. Hosts were compared in memory; addresses and
credentials are not recorded. These are boundary observations, not a continuous
network trace or verification of remote CPU activity. The idle shell and desktop
are user-confirmed conditions. Both connections were left open after the tests.
Every test child exited and closed its audio device before writing its report.

## Comparison

| Fixed post-reboot condition | M6 | M7 | Skipped UI frames |
| --- | --- | --- | ---: |
| Disconnected | 3/3 passed | 3/3 passed | 0 |
| Screen Sharing connected, no SSH | 3/3 passed | 3/3 passed | 0 |
| Screen Sharing plus idle SSH | 3/3 passed | 3/3 passed | 0 |

`comparison-three-conditions.json` verifies the common recorded configuration and
links the two earlier summaries by hash. All eighteen 100-second runs passed.
This does not support an idle Screen Sharing connection, with or without idle
SSH, as the sole explanation under these post-reboot conditions. Active remote
work, changing remote-screen content, ordinary desktop interaction and other
session/load conditions were not tested here. The cause of the original skips
remains unresolved; these samples do not establish a permanent fix.

`preflight.json`, each invocation/outcome, the raw reports/logs and
`final-integrity.json` retain exact identities and observations. All 155 M7 source
inputs and 339 packaged artifacts, both executable hashes, branch HEAD and all
existing local tags remained unchanged. Earlier passing and failed evidence is
preserved. No M7 commit, tag, push or merge occurred. Remaining required physical
M7 checks are still deferred; full physical 7.1 listening remains separately
pending. This addendum has its own evidence hash manifest.
