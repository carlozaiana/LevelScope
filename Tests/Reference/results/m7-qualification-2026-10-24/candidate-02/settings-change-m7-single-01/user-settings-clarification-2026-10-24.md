# User clarification of settings — 24 October 2026

This note was added after the completed single M7 test. The original reports,
README, comparison context and evidence manifest retain what was known at their
collection time; their earlier “description outstanding” status is resolved by
this later user clarification.

The user reports that, before requesting this test, they:

- Granted Firefox Local Network access on the Mac Studio.
- Opened some browser tabs.

The user also reports that Sequoia blocks this access again after every reboot,
and that they normally need it to check signal strength through their local
router/receiver. The recurring permission behavior is a user observation; it was
not independently reproduced or diagnosed in this qualification work.

These changes describe the environment of the passing 100-second M7 run at
approximately 18:23:52–18:25:39 Europe/Rome on 24 October 2026. The report retained
zero skipped frames, underruns, blank paints, clock-bound failures and CPU-paint
deadline failures. Screen Sharing and idle SSH remained connected. The recorded
display/power snapshots being unchanged is compatible with these browser and
permission changes, which those snapshots do not inventory.

Exact tab count/content, browser resource use and network activity during the
run were not measured. Do not assume Firefox was absent from earlier runs, infer
that the permission reset caused the original paint skips, or treat this single
pass as proof of a fix. It records a passing run with the user-reported setup.

No additional tests, remote commands, permission changes or source edits were
performed for this clarification. M7 remains frozen and unqualified overall;
its remaining required checks and original failed evidence are unchanged.
