# Compliance Matrix (Foundation)

Purpose: map key loudness/compliance requirements to current implementation status and validation coverage so DSP/proposal work can be checked systematically.

Status legend: **Implemented** = present in current code/docs, **Partial** = some support but not full profile compliance, **Planned** = identified but not implemented.

## EBU R128 / ATSC A/85 Requirement Mapping

| Requirement ID | Standard/Profile | Current status | Validating method | Existing test/check reference | Gaps / next action |
|---|---|---|---|---|---|
| CM-LOUD-001 | EBU R128 + ATSC A/85: program loudness measurement path aligned to BS.1770-family method | **Implemented** (measurement primitives available) | Manual | Architecture snapshot loudness primitives + current input-side analysis note; ProjectIndex standards coverage | Add deterministic numeric regression fixtures (44.1/48/96/192 kHz) and pass/fail tolerances |
| CM-LOUD-002 | EBU R128: Momentary (M), Short-term (S), Integrated (I) visibility | **Implemented** | Manual | Architecture snapshot: loudness/features computation and running stats | Add automated assertions for windowing correctness and continuity across transport seeks |
| CM-LOUD-003 | EBU R128: Loudness Range (LRA) support | **Implemented** (running/LRA views available) | Manual | Architecture snapshot running LRA + rolling LRA controls | Add offline reference-vector verification for LRA algorithm against known material |
| CM-TP-004 | EBU R128 / ATSC A/85 operational practice: output peak protection (true-peak-oriented limiter path) | **Partial** (limiter path + oversampling option, no full true-peak compliance claim yet) | Manual | ProjectIndex standards coverage: limiter true-peak protection path | Define explicit true-peak conformance target and add oversampling validation vectors |
| CM-PROFILE-005 | Profile targeting: EBU vs ATSC preset/profile behavior | **Planned** | Manual | ProjectIndex standards roadmap (profiles prioritized, not yet productized) | Add formal profile schema (targets, tolerances, gating/report outputs) |
| CM-REPORT-006 | Compliance reporting/export artifacts | **Planned** | Manual | No explicit report pipeline documented | Specify report contract (measured I/LRA/TP, pass/fail, source metadata, render settings) |
| CM-TRANSPORT-007 | Measurement robustness across start/stop/seek/loop discontinuities | **Implemented** (guard rails in plugin transport logic) | Manual | Architecture snapshot transport ramp guard, stop freeze, discontinuity detection | Add scenario checklist and automated host-transport simulation harness |
| CM-SR-008 | Timing/sample-rate invariance for measurement/control behavior | **Partial** (sample-rate aware but mixed fallback defaults) | Manual | Sample rate & timing audit (48k defaults + 44.1 fallback noted) | Standardize fallback policy and add multi-rate conformance test matrix |
| CM-MULTI-009 | Multichannel policy consistency (channel roles, LFE handling) | **Partial** (up to 7.1.4 support + current LFE policy choices) | Manual | Architecture snapshot multichannel scope and default LFE policy | Split “detector policy” vs “compliance profile” defaults and add profile-specific checks |
| CM-STATE-010 | Reproducible state/preset behavior for compliance workflows | **Partial** (additive chunks + settings preset separation implemented) | Manual | Architecture snapshot persistence strategy + APVS/MODG/UIST contracts | Add round-trip state compatibility test set and profile-lock metadata |

## Validation Foundation Notes

- **15 October 2026 standalone update:** a separate offline kernel now passes 39
  official ITU BS.2217 integrated-loudness vectors (±0.1 LU) and 184 generated/integration
  checks on arm64 Sequoia; the original 205 scaffold checks also pass. This evidence
  applies to the new standalone measurement path, not the existing plugin algorithms.
  The supplied official EBU v05 suite now also passes: 66/66 minimum-requirements files
  (72 rules, 1,052 readings), plus 4/4 supplementary alignment sanity checks. Eight
  runner failure-detection tests pass. No DSP correction was required. All official
  EBU files are 48 kHz; live-signal cases test the offline kernel, not a live GUI.
  Evidence: `../Tests/Reference/results/ebu-v05-sequoia-arm64-2026-10-15.json` and the
  pinned manifest/runner in `../Tests/Reference` / `../Tests/run_ebu_measurements.py`.
  Intel High Sierra measurement execution and interactive UI measurement testing
  remain outstanding. See `ProjectIndex/13_StandaloneScaffoldHandover.md`
  and `14_OfflineMeasurementContract.md`. No EBU/ATSC delivery-compliance verdict is
  implemented. Earlier implemented/manual rows above describe historical primitives,
  not completed numerical qualification of those algorithms.

- Historical plugin evidence remains documentation- and implementation-contract-based.
  The standalone M1 kernel now has the automated numeric evidence above. Scaffold tests
  separately validate media/routing/lifecycle plumbing. No standalone result qualifies
  the older plugin's measurement algorithms or enables delivery-compliance verdicts.
- Short-term recommendation: keep this matrix updated whenever DSP behavior, transport semantics, sample-rate policy, or preset/state behavior changes.
- Medium-term recommendation: add machine-readable mirror (YAML/JSON) of this table to drive automated checks and proposal gating.
