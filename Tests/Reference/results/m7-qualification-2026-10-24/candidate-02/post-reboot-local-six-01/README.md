# Six local reference-movie tests after reboot — 24 October 2026

**All six predeclared tests passed. M7 remains frozen and unqualified overall.**
These are additional diagnostic observations; they do not overwrite the earlier
failed attempts, complete the deferred physical matrix/endurance checks, or
establish physical 7.1 listening qualification.

The user rebooted the Mac Studio, disconnected Golem, closed Screen Sharing and
agreed to leave the Mac untouched. The local boot record is 17:00:55 Europe/Rome;
testing ran approximately 17:06:22–17:17:00. Preflight and final local process
checks found no SSH client, Screen Sharing client or other candidate app. No
remote command or connection was used. The test app was brought forward by the
unchanged harness; no computer-use inspection or parallel builds/tests ran.

Actual platform: arm64 Sequoia 15.7.9 (24G830), Apple M1 Max. The read-only display
snapshot records PB248 at 1920×1200, 60 Hz. Mac Studio Speakers operated at 48 kHz
with every route muted. Existing power/display/audio configuration was preserved.

The fixed sequence alternated preserved M6 and frozen M7 stereo, three runs of
each. Every run used the same local `51_SyncTest.mov`, stereo clock WAV and
100-second duration as the earlier reference comparison. No binary rebuild,
threshold change or retry-until-pass occurred. Binary/media hashes and exact
arguments are in `preflight.json` and each `*-invocation.json`; original raw
reports and logs are preserved beside this file.

| Run | Elapsed seconds including preparation/shutdown | Result | Skipped UI frames |
| --- | ---: | --- | ---: |
| m6-stereo-1 | 106.864 | Pass | 0 |
| m7-stereo-1 | 106.326 | Pass | 0 |
| m6-stereo-2 | 106.292 | Pass | 0 |
| m7-stereo-2 | 106.321 | Pass | 0 |
| m6-stereo-3 | 106.286 | Pass | 0 |
| m7-stereo-3 | 106.316 | Pass | 0 |

Total elapsed test time was 638.405 seconds (10 minutes 38 seconds). Every report
also recorded zero audio underruns, blank paints, frame-selection bound failures
and CPU-paint deadline failures; hardware audio timestamps were active. Existing
attachment and memory limits passed. Maximum observed paint intervals ranged
from 29.828 to 38.317 ms across the six runs.

`final-integrity.json` confirms unchanged M6/M7 executable hashes, all 155 M7
source inputs and 339 packaged artifacts, branch HEAD and all existing local tag
objects/targets. Every child process exited; the harness disconnects its device
before writing the report. No new tests are scheduled. No M7 commit/tag/push/merge
was made.

## Interpretation

The prior failure did not reproduce in this fixed post-reboot, disconnected
six-run sample. This supports investigating session/workload conditions before
changing the video code. Rebooting and disconnecting Screen Sharing/SSH changed
multiple conditions together; these results do not identify Screen Sharing,
SSH, Golem's age/OS, or any particular process as the cause. Three passes per
binary also do not prove that an intermittent failure has been eliminated.

A separately coordinated comparison that changes one condition at a time can
help isolate the cause; ordinary-interaction robustness should remain a separate
case. M7 stays frozen pending interpretation and its remaining required checks.
The prior evidence inventory and qualification summary remain historical records;
this addendum has its own hashes in `evidence-manifest.json`.
