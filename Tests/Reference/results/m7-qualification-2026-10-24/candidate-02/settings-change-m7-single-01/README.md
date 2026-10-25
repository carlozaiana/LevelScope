# One M7 reference test after user-reported settings changes — 24 October 2026

**The single 100-second muted M7 stereo test passed. M7 remains frozen and
unqualified overall.** This is separate from the three earlier six-run samples.

The user requested one further test after changing settings on the Mac Studio.
The assistant selected frozen M7 candidate-02, stated this scope, and requested
the exact settings for annotation. That description had not yet arrived at
evidence collection. The user was told to wait until playback finished before
replying so the run could remain untouched.

Playback ran approximately 18:23:52–18:25:39 Europe/Rome; executable elapsed time
was 106.928 seconds. The same preserved movie and stereo clock WAV were used,
with Mac Studio Speakers at 48 kHz and every route muted. Actual execution was
arm64 Sequoia 15.7.9 (24G830). Screen Sharing and SSH connections to matching
hosts were observed at preflight and run boundaries; no remote command was sent.

The report records zero skipped UI frames, audio underruns, blank paints,
frame-selection bound failures and CPU-paint deadline failures. Attachment and
memory limits passed, hardware audio timestamps were active, and the maximum
observed paint interval was 25.227 ms. The test process exited and closed its
audio device; remote connections were left open.

Recorded boot, OS/CPU, display and power settings match the preceding idle-SSH
comparison. These observations do not cover every macOS/application setting and
do not contradict the user's reported change. `comparison-and-context.json`
marks the exact change description as outstanding. This single pass does not
establish an effect of the changed settings or resolve the original failure.

`preflight.json`, the invocation, raw report/log, outcome and final integrity
record retain exact source/media/binary identities. All 155 frozen source inputs,
339 packaged artifacts, both preserved executable hashes, branch HEAD and local
tags remained unchanged. No rebuild, code edit, extra test, commit, tag, push or
merge occurred. Earlier evidence remains preserved. A later description of the
settings should be added as a separate dated note without replacing these records.
