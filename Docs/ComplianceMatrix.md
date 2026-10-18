# Compliance Matrix (Foundation)

Purpose: map key loudness/compliance requirements to current implementation status and validation coverage so DSP/proposal work can be checked systematically.

Status legend: **Implemented** = present in current code/docs, **Partial** = some support but not full profile compliance, **Planned** = identified but not implemented.

## Historical plugin primitives / roadmap (not standalone M1 qualification)

| Requirement ID | Standard/Profile | Current status | Validating method | Existing test/check reference | Gaps / next action |
|---|---|---|---|---|---|
| CM-LOUD-001 | EBU R128 + ATSC A/85: program loudness measurement path aligned to BS.1770-family method | **Implemented** (measurement primitives available) | Manual | Architecture snapshot loudness primitives + current input-side analysis note; ProjectIndex standards coverage | Add deterministic numeric regression fixtures (44.1/48/96/192 kHz) and pass/fail tolerances |
| CM-LOUD-002 | EBU R128: Momentary (M), Short-term (S), Integrated (I) visibility | **Implemented** | Manual | Architecture snapshot: loudness/features computation and running stats | Add automated assertions for windowing correctness and continuity across transport seeks |
| CM-LOUD-003 | EBU R128: Loudness Range (LRA) support | **Implemented** (running/LRA views available) | Manual | Architecture snapshot running LRA + rolling LRA controls | Add offline reference-vector verification for LRA algorithm against known material |
| CM-TP-004 | EBU R128 / ATSC A/85 operational practice: output peak protection (true-peak-oriented limiter path) | **Partial** (limiter path + oversampling option, no full true-peak compliance claim yet) | Manual | ProjectIndex standards coverage: limiter true-peak protection path | Define explicit true-peak conformance target and add oversampling validation vectors |
| CM-PROFILE-005 | Standalone selected-rule evaluation | **Partial** | Automated + silent UI qualification | `ProjectIndex/21_M2MeasuredCompliance.md` | M2 implements one bounded EBU programme QC profile; ATSC, general editing and full delivery/export compliance remain outside scope |
| CM-REPORT-006 | Compliance reporting/export artifacts | **Planned** | Manual | No explicit report pipeline documented | Specify report contract (measured I/LRA/TP, pass/fail, source metadata, render settings) |
| CM-TRANSPORT-007 | Measurement robustness across start/stop/seek/loop discontinuities | **Implemented** (guard rails in plugin transport logic) | Manual | Architecture snapshot transport ramp guard, stop freeze, discontinuity detection | Add scenario checklist and automated host-transport simulation harness |
| CM-SR-008 | Timing/sample-rate invariance for measurement/control behavior | **Partial** (sample-rate aware but mixed fallback defaults) | Manual | Sample rate & timing audit (48k defaults + 44.1 fallback noted) | Standardize fallback policy and add multi-rate conformance test matrix |
| CM-MULTI-009 | Multichannel policy consistency (channel roles, LFE handling) | **Partial** (up to 7.1.4 support + current LFE policy choices) | Manual | Architecture snapshot multichannel scope and default LFE policy | Split “detector policy” vs “compliance profile” defaults and add profile-specific checks |
| CM-STATE-010 | Reproducible state/preset behavior for compliance workflows | **Partial** (additive chunks + settings preset separation implemented) | Manual | Architecture snapshot persistence strategy + APVS/MODG/UIST contracts | Add round-trip state compatibility test set and profile-lock metadata |

## Validation Foundation Notes

- **M1 qualification completed, 17 October 2026:** actual Intel High Sierra 10.13.6
  and arm64 Sequoia 15.7.9 both passed 184 measurement checks, 205 scaffold checks,
  eight runner tests, all 39 pinned ITU files, 66 EBU normative files (72 rules,
  1,052 readings), four supplementary files and the silent interactive UI checklist.
  Cross-platform comparisons passed. See `ProjectIndex/17_M1QualificationResults.md`
  for exact identities and `ProjectIndex/20_MainIntegrationAcceptance.md` for the
  accepted video integration and unchanged-measurement regressions. This evidence
  qualifies the standalone offline kernel, not historical plugin meters.
- **M2 qualification completed, 18 October 2026:** versioned EBU R128-2023 whole-programme
  linear-audio QC and measured Source/Current State rule evaluation pass on actual
  Intel High Sierra and arm64 Sequoia, including 128 M2 checks, existing regressions,
  both official reference suites and silent UI checks. The final display correction
  was rechecked on both platforms; all automated suites were rerun.
  See `ProjectIndex/21_M2MeasuredCompliance.md` for applicability, identity, exact
  comparisons and the current qualification status. A scalar rule pass is not a
  broadcaster delivery, metadata or encoded-export compliance claim.

- Historical plugin evidence remains documentation- and implementation-contract-based.
  The standalone M1 kernel now has the automated numeric evidence above. Scaffold tests
  separately validate media/routing/lifecycle plumbing. No standalone result qualifies
  the older plugin's measurement algorithms or enables delivery-compliance verdicts.
- Short-term recommendation: keep this matrix updated whenever DSP behavior, transport semantics, sample-rate policy, or preset/state behavior changes.
- Medium-term recommendation: add machine-readable mirror (YAML/JSON) of this table to drive automated checks and proposal gating.
