# Standards Coverage & Roadmap Matrix

M13 qualification update, 28 October 2026: bounded candidate-02 is implemented
and qualified on actual Sequoia/High Sierra; [ProjectIndex 47](47_M13BoundedSectionTrimProposals.md) is authoritative for this scope. The broader
no-change rule covers all applicable explicitly selected, method-defined objectives;
M13 still evaluates only supported programme rules. The intended cinema comparison
is user-reported Nuendo statistics supplied 28 October, while the historical screenshot
is Cubase 10 supplementary evidence. No exact Nuendo version is claimed. The user authorized checkpoint publication
and combined M12/M13 main integration on 28 October. The preserved M12 qualification
status remains unchanged; [49](49_M13QualificationAndHandover.md) records the completed combined
M12/M13 qualification. The [published M13 checkpoint](https://github.com/carlozaiana/LevelScope/tree/m13-bounded-section-proposals-2026-10-28)
identifies the combined main revision at publication.

Design/status update: 28 October 2026. This document distinguishes qualified
standalone behavior from historical plugin capabilities and future proposals.

## Current standalone baseline

Main receives the qualified combined M12/M13 [checkpoint](https://github.com/carlozaiana/LevelScope/tree/m13-bounded-section-proposals-2026-10-28)
by fast-forward from M11 `cfb1033aa7fb65a2655db276bee696161196bdeb`.
M1 measurement, M2 evaluation, M3 global trim and M4–M11 diagnostics, sections,
audition, sessions, assembly, verified PCM export and composition have their own
published evidence; see the [index](README.md). The authoritative EBU R128-2023
whole-programme linear-audio QC profile retains -23 LUFS intended target, ±0.2 LU
allowance and strict -1 dBTP ceiling. LRA is informative, with no universal ceiling.
No engine-test or meter tolerance silently widens these profile limits.

Published M12 candidate-03 is **“Sequoia tested; High Sierra and video qualification
pending.”** That historical tag remains unchanged; M12 processing is included in the combined M13 main checkpoint. [46](46_M12QualificationAndHandover.md)
records its actual evidence and preserves the frozen candidate. The approved
bounded method in [47](47_M13BoundedSectionTrimProposals.md) is implemented in the
qualified M13 checkpoint. [49](49_M13QualificationAndHandover.md) records its
completed combined qualification and the passing final required-evidence gate.

Historical plugin loudness/processing code, including upward/downward processing
and a limiter, does not establish qualified standalone DSP or VST3 measurement
parity. Standalone integration/measurement/audition/export qualification must precede
any new processor's compliance claim. No plugin installation changes are implied.

## Shared strategy and standards boundaries

The [proposal specification](../ProposalEngineSpec_v2.md) defines dynamics-led,
section-led and mixed processing strategies with separate algorithms/preferences.
They share profiles, exact measurement, processing, identity and verification/export
gates. All preserve an already-passing result and consider global-only correction
first; subsequent preferences depend on strategy. Planning reserves and preservation
budgets never relax required rules. Predictions, private measured candidates,
Current measurements and verified artifact measurements remain separate.

ATSC A/85, other distribution profiles, dialogue-gated methods and intelligibility
objectives remain separately defined/validated work. A loudness-region boundary,
measurement gate and compressor threshold are different concepts, with explicit
mapping needed for future processing. No broad broadcaster-delivery claim follows
from passing the selected programme loudness/true-peak rules.

## Dependency order after the bounded M13 method

1. Mono section processing/export and optional basic marker import, each separately
   reviewed; bounded EDL import follows later without becoming a proposal prerequisite.
2. Qualify standalone limiter, downward and upward processing in small stages.
3. Dynamics-led proposals, then enhanced section-led proposals, then mixed strategy.
4. Later component-relative dialogue/M&E processing and validated dialogue/intelligibility
   objectives. Keep analysis-only references independently identified and non-audible.

Generate only the selected strategy; explicit comparisons reuse compatible analysis
and share a frozen measured baseline with bounded resources. No detailed future
milestone numbering is committed. M13 development remains Sequoia-first, followed
by combined M12/M13 actual High Sierra and selected deferred video qualification.
Accepted occasional picture skips retain all existing audio/drift/lifecycle gates.
Both combined qualification and the subsequent publication/integration approval are complete.

Detailed requirements, provenance and evidence remain in
[the measurement contract](14_OfflineMeasurementContract.md),
[the M2 profile contract](21_M2MeasuredCompliance.md),
[the M10 export handover](40_M10QualificationAndHandover.md) and
[the compliance matrix](../ComplianceMatrix.md). No tests are run by this design update.
