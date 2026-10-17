# Separate milestone: VST3 measurement parity

Recorded 17 October 2026. This is future work; M1 does not update, build or install
the plugin, and its standalone numerical qualification does not qualify plugin meters.
Use `measurement-milestone-2026-10-17` as the preserved standalone reference; see
[17_M1QualificationResults.md](17_M1QualificationResults.md) for its exact scope.

## Objective and preserved behavior

Bring VST3 measurement semantics into agreement with the qualified M1 kernel for the
same samples, roles, sample rate and explicit measurement interval. Preserve parameter
IDs/ranges/defaults, APVTS serialization, existing saved-state chunks and migrations,
automation interpretation, bus negotiation, bypass, latency reporting and processing
behavior. Metering changes must not alter output audio or dynamics decisions.

## Integration boundaries

- Audit `PluginProcessor` taps and `ProcessContext`: distinguish pre-processing input,
  post-processing output, monitor-only gain and transport ramps. Name each metered
  domain and use the same sample domain when comparing with M1.
- Keep detector/control filters and dynamics sidechains separate from compliance
  metering. Replacing plugin K weighting or statistics must not change a module's
  envelope, gain computer, limiter or LFE application policy.
- Adapt ordered host buses to explicit measurement roles, not channel-count guesses.
  M1 supports conventional mono/stereo/5.0/5.1; extended/height layouts require their
  own standards mapping and reference qualification before parity can be claimed.
- The offline kernel currently allocates growing history vectors and performs final
  sorting/draining. Do not place it unmodified in the audio callback. Design bounded,
  allocation-free capture and worker processing, with explicit overflow/discontinuity
  handling, or qualify a separate real-time accumulator against the same numerical core.
- Define start/pause/resume/reset, host play/stop/seek/loop, bypass, rate/block-size and
  bus-layout changes. An offline EOF drain must not be mistaken for a live pause or
  inserted into audible output. Reconcile gate populations, exact M/S windows, LRA
  cadence/percentiles, true-peak boundary state and metric-unavailable reasons.
- Keep measured values distinct from predictions and UI history decimation. Introduce
  meter version/interval/generation identity without reinterpreting old preset state.

## Acceptance tests

1. Feed identical pinned ITU/EBU and generated multirate samples through the plugin
   measurement adapter and M1. Require published per-test tolerances and paired
   differences no larger than 0.01 LU/dB (I/M/S/peaks) and 0.1 LU (LRA), with equivalent
   interval/boundary semantics and reasons for missing metrics.
2. Randomized block sizes, host offline/realtime rendering, repeated runs, channel
   permutations, LFE-only, silence, short inputs, intersample peaks and tail boundaries.
   Qualify additional layouts separately; never silently map heights/rears as 5.1.
3. Host transport discontinuities, loop ranges, pause/resume/reset, sample-rate changes,
   editor close/reopen and multiple instances. Stale UI/worker generations never publish.
4. Realtime instrumentation: no allocation, locks, filesystem/network access or
   unbounded processing on the audio callback; bounded memory and declared overflow
   behavior. Stress the slow worker/closed editor paths without affecting playback.
5. Compare rendered output before/after with fixed plugin states and automation,
   including dynamics/LFE/bypass modes. Require unchanged processing samples, latency
   and parameter behavior; investigate any difference rather than accepting meter fixes
   as a reason for changed sound.
6. Golden parameter-ID inventory, old preset/session round trips, state chunk retention,
   bus-layout restoration, automation and host reopening checks on High Sierra/Sequoia.
   Use isolated test plugin locations and obtain authorization before any installation.

Completion requires separate plugin/host evidence and a dedicated checkpoint. M1's
tag is a preserved standalone reference, not an assertion of VST3 parity.
