# Parameter Capabilities Matrix

Purpose: baseline capability map for automation safety, playback-lock behavior, latency impact, and state scope.

Legend:
- Playback-lock behavior: **Locked** (disabled during effective playback) / **Live** (can change during playback) / **Hybrid** (mixed by sub-control).
- Scope: **Preset** = included in `.lscsettings` behavior (APVS + MODG), **Snapshot** = broader project/session state context.

| Parameter ID / Group | Automatable | Playback-lock behavior | Latency-affecting | Scope (preset vs snapshot) | Compatibility risk notes |
|---|---|---|---|---|---|
| `lvlr.enabled`, `mtdm.enabled` (module enable/bypass families) | Yes (typical enable controls) | Live (except if routed to structural/latency-changing behavior) | Usually No | Preset | Medium: behavior coupling with audition/routing can change perceived gain staging |
| `lvlr.targetLufs`, `lvlr.maxBoostDb`, `lvlr.maxCutDb`, `lvlr.rate*` | Yes | Live | No | Preset | Low/Medium: automation abruptness may cause audible transitions if smoothing assumptions change |
| `lvlr.measChoice`, `lvlr.modeChoice`, `lvlr.learn01`, `lvlr.controlModeChoice`, `lvlr.captureToHost01` | No (explicit non-automatable controls documented for these families) | Hybrid (mode toggles may be playback-guarded by UI policy) | No | Preset | Medium: changing enum semantics breaks preset portability if value mapping shifts |
| `lvlr.hostGainDb` | Yes | Live | No | Preset | Low: standard continuous parameter; ensure range/default stability |
| MTDM thresholds `t0Lufs..t3Lufs` | Yes | Live | No | Preset + Snapshot-visible effect | Medium: ordering/min-gap invariants must remain stable for old sessions |
| Upward stage dynamics controls (amount/maxBoost/curve/knees/attack/release) | Mixed (core continuous controls expected automatable; some advanced choices may be non-automatable) | Hybrid | Potentially (if implementation switches structure/FFT mode) | Preset | High: advanced mode/curve definition changes can alter session recall sonics |
| Upward advanced analysis structure (`fftSizeChoice`, `bandsPerOctChoice`, min/max freq, calibration trim, curve type) | Mixed | Locked for structural controls during effective playback | **Yes** for structure-dependent controls | Preset | High: structural latency/response changes require strict versioned migration notes |
| Downward stage controls (`down*`) | Yes (continuous controls), bypass/enable bools yes | Live | No | Preset | Medium: ratio/knee law updates can change backwards sonic equivalence |
| Limiter controls (`limCeilingDb`, `limLookaheadMs`, `limAttackMs`, `limReleaseMs`, `limDriveDb`, `limOversamplingChoice`, enable/bypass) | Mixed | **Locked** for latency/structure-sensitive edits during effective playback | **Yes** (lookahead / oversampling choices) | Preset | High: latency reporting must remain synced with parameter changes |
| Multichannel policy (`mcPolicyChoice`, `dialogDetectorChoice`, `dialogApplyChoice`, `lfeInDetector`, `lfeInApply`) | Mixed (policy enums often non-continuous and may be non-automatable) | Hybrid (policy changes may be guarded) | Usually No (unless future structure switch) | Preset + Snapshot-context implications | High: policy meaning drift can invalidate compliance comparisons across versions |
| Zone audition controls (`zoneAud.*`, solo/mute combinations) | Mixed | Live | No | Preset (if parameterized) + Snapshot workflow impact | Medium: audition logic coupling risk with bypass/threshold interaction |
| UI-only layout/history controls (non-APVTS `UIST` state) | N/A (not plugin automation parameters) | N/A | No | Snapshot/UI state, not settings preset | Medium: keep isolated from APVS to avoid accidental preset contamination |

## Operational notes

- Exact per-parameter automatable flags are defined at registration time in `createParameterLayout()` and should be treated as source-of-truth when this matrix is revised.
- Playback-lock behavior is based on effective playback policy (transport playing + callback freshness), not only host play flag.
- Any change to parameter IDs, enum ordering, or latency behavior should trigger updates to this matrix and compatibility notes.
