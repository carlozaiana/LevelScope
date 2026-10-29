# M14 pre-implementation impact and qualification plan

Recorded 29 October 2026, before source implementation. Approved contract:
[50](50_M14MonoProcessingAndVerifiedExport.md). Base main/logic commit
`a9be070ab9eb5cb4f591ff9d991d8999c3a016d0`. Main stays there.

| Affected behavior / dependency | Fresh Sequoia evidence required |
| --- | --- |
| Mono section-layout gate; existing M12/M13 processor and callers | Four-rate L/R/C mono manual and proposal workflow; exact independent samples; transition/seek/EOF; no-op, policy failures, stale/apply/cancel safeguards. Focused M12/M13 stereo/5.1 regressions. |
| PCM24 format/writer/strict inspector; export mapping, reports | Odd/even mono frames, masks, padding, exact PCM+dither oracle, strict artifact rejection, boundary/sparse >4 GiB cases, disk/write/cancel/publication failures; existing stereo/5.1 export suite. |
| Separate RF64 input policy; reader/identity/preparation/session callers | Synthetic extended ds64/zero count/trailing metadata, malformed sizes/alignment/masks/counts/chunks, complete read/seek samples; original cinema versus M13 copy at bounded positions past 4 GiB. Existing grouped/composed input smoke, pinned EBU/ITU and measurement/scaffold checks. |
| Session/diagnostic/audition UI integration and identity | Mono roundtrip, old M13 reader rejection, cache recovery, stale jobs, monitor-only changes, audio-only injected playback and overload, silent Sequoia UI; video alignment data only, no video decoding. |
| Named target/build changes | Universal standalone and explicit test targets, arm64 execution, Intel slice/deployment inspection; no VST3/default-all/install. |
| Sustained changed workload | One 90-minute 192 kHz mono workflow: sections/proposal/Apply/measurement/diagnostics/session/audition captures/verified RIFF export/resources; <=120 min execution. Independent short PCM oracle, bounded long oracle + complete hash. |

Runner phases must enumerate required checks, record argv/exit status/log hashes,
OS/CPU/architecture and candidate/build identity, and fail missing or duplicated
evidence. No recursive historical runners. Re-run affected checks after fixes;
do not substitute development or a rebuilt binary for frozen evidence.

Reuse sources: M13 candidate-02 source manifest
`50eba19a2d66096607ec878c8eaa05e610580c7f5e8854d7b375d846a89c5f27`,
build `5fd55bb1a81c80b4b556715936b6873bccae29a5d04c00b0890a9b43296a20d2`;
its unchanged arithmetic/policy, broad representation, full video/long and cinema
measurement evidence remain historical. Each reused receipt must identify the
original candidate and unchanged dependencies. New mono and RF64 input claims
require fresh tests. Cinema equivalence may reuse the archived PCM payload hash
and source measurements after verifying original identity and bounded reader parity;
no new full cinema proposal/export run.

Budgets and estimates are in document 50: >=20 GiB working headroom plus measured
size preflight; <=256 MiB additional working RSS, <=100 ms UI p95, <=1 s cancellation,
<=5 s save/parse, <=128 MiB session and <=16 MiB reports. Estimate 30-70 minutes for
the one sustained workflow and 45-90 minutes short/reference/UI work, excluding
build/fixes/operator delays. Track real bytes separately from sparse logical size.

Explicitly deferred: Golem/actual High Sierra runtime and UI, cross-Mac exchange,
all video load/decode/playback, new three-hour mono/RF64 end-to-end qualification.
Combined target window is after M18, no later than M19, using one frozen combined
candidate. Earlier macOS API/dependency, measurement-arithmetic, buffering/clock/
latency changes or concrete failures require an explained impact-based revision.
No routine Golem access or physical-device test is authorized by this plan.

Completion means **Sequoia phase complete; High Sierra, video and full-duration
mono/RF64 qualification pending**, never fully qualified. Stop for readiness review
before publication/integration/next milestone. Preserve all frozen packages/tags,
the video working copy, original media and installations/audio configuration.

Subsequent review: the user authorized candidate-02 development-checkpoint
publication on 29 October. Qualification deferrals and the prohibition on main
integration remain in force; see [52](52_M14QualificationAndHandover.md).
