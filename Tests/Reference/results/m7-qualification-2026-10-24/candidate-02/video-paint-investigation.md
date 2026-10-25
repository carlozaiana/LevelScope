# Candidate-02 intermittent paint skips — reproduced; handed over

This is a qualification gap, not a passing video result. No threshold or candidate
implementation was changed during this investigation.

The first Sequoia all-muted Source ProRes matrix completed 16 cases: 14 passed,
while Proxy 30000/1001 and 30 fps failed the zero-skipped-frame gate with two and
one skipped UI frames respectively. All cases had zero audio underruns, blank
paints, frame-selection bound failures and CPU-paint deadline failures. Each failed
case had queued decoded pictures (7–8) at the skip events. The observed UI paint
gaps were approximately 49–50 ms; frame PTS selection matched the current audio
presentation clock when painting resumed. This does not excuse a failed gate.

The user could not confirm whether desktop/window interaction occurred during
that first group. Its cause is unknown; user interference is not established.
The user correctly emphasized robustness during ordinary interaction. These
records distinguish audio continuity and clock synchronization from every-frame
visual presentation. Neither a controlled test nor a stereo listening check is
an unrestricted guarantee under arbitrary desktop or system load.

The runner stopped after the failed Source matrix, so its Current, H.264 and
endurance groups did not run. Original evidence is retained in
`sequoia-muted-video-first-attempt`.

One bounded diagnostic comparison followed, with no retry-until-pass loop:

| Preserved binary / selected source | Movie | Outcome |
| --- | --- | --- |
| Qualified M6 / stereo clock WAV | Proxy 30000/1001 | Pass, zero skipped frames |
| Candidate-02 / stereo clock WAV | Proxy 30000/1001 | Pass, zero skipped frames |
| Candidate-02 / eight-channel clock FLAC | Proxy 30000/1001 | Pass, zero skipped frames |
| Candidate-02 / eight-channel clock FLAC | Proxy 30 | Pass, zero skipped frames |

All used the approved Mac Studio Speakers at 48 kHz with every route muted. Exact
binary hashes, commands and reports are retained in `paint-investigation-01`.
These passing comparisons demonstrate intermittency, not a root cause or fix.

At approximately 14:12 UTC the user was given an explicit hands-off notice before
one complete controlled Sequoia repeat, leaving the window visible and using the
same frozen candidate and unchanged gates. All fifteen generated Source cases
passed. The reference movie failed with one skipped UI frame: 1419 → 1421 at
audio 56.842292 s, a 50.023709 ms paint gap, seven pictures queued and no worker
selection skips. Audio underruns, blank paints and clock/CPU-paint bound failures
remained zero. The completed failed attempt is retained separately as
`sequoia-muted-video-controlled` in this evidence directory.

One final bounded 100-second reference-movie comparison reproduced the failure
with the unchanged qualified M6 binary and stereo audio (two UI frames skipped)
and M7 with stereo audio (one). These reports and exact commands are in
`paint-reference-investigation-01`. Thus eight-channel processing is not required
to reproduce the issue. A main-thread scheduling explanation is a hypothesis,
not a proven cause. Earlier M6 qualification reports remain unchanged.

The downstream physical Current/H.264/endurance groups were deferred after this
reproduction. High Sierra's already-started Source matrix finished normally:
all sixteen cases passed with zero skipped frames, including the reference movie.
Its reports are retained in `high-sierra-muted-video/prores-matrix-source`.
The user explicitly chose
to keep M7 frozen and prepare a Video Implementation handover. No further video
code changes or qualification groups are started here. M7 remains unqualified.

The concrete handover is
[ProjectIndex 33](../../../../../Docs/ProjectIndex/33_M7VideoPaintTimingHandover.md).
