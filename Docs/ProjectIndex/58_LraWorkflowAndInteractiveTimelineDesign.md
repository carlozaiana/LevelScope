# Selective LRA workflow and interactive timeline — design priorities

30 October 2026 user direction. Architectural preparation and roadmap, not authorization to implement these features in M16.

Intended workflow: **Open mix and optional video → automatic analysis → select delivery requirements → Propose → inspect, audition and refine → export the verified result.** Selective LRA correction is a central product objective. Preserve dynamics-led, section-led and mixed strategies, but a complete upward/downward compressor is not a prerequisite for the section-led LRA/timeline path.

## Proposed dependency order after M16 (each requires review)

1. **M17 approved: explicit job LRA requirements and reusable named presets.** [ProjectIndex 60](60_M17JobLraRequirementsAndPresets.md) specifies the approved bounded implementation and acceptance tests; implementation/Sequoia qualification approved 30 October 2026; [candidate-04 Sequoia phase complete](62_M17QualificationAndHandover.md), development checkpoint published as `m17-job-lra-presets-sequoia-2026-10-30`; main remains M13. Plain EBU R128 retains no universal LRA ceiling. Ask inline whether the job requires a maximum LRA, without recurring popups. Unset means “no LRA acceptance limit specified.” An added ceiling has its own provenance and base-profile identity. All applicable selected objectives must pass: IL/TP alone cannot yield already compliant with an unmet LRA requirement. Define method, short-programme interpretability, exact thresholds, identity and acceptance before implementation. Configuring a requirement does not invent a correction strategy.
2. **M18 proposed for review: useful interactive audio timeline and progressive analysis.** [ProjectIndex 63](63_M18InteractiveTimelineAndProgressiveAnalysis.md) defines the bounded first slice and acceptance tests; no implementation is authorized. Deliver a bounded slice using existing processing: waveform/measurement curves, markers/boundaries, editable gain envelopes, playhead/seek, audition and original/processed comparison. Keep edits/drafts/applied/verified states explicit. Start analysis automatically after programme and roles are established; ambiguous roles require confirmation. Progressive waveform/curves are unfinished until final authority exists. Distinguish rolling LRA (LU), cumulative IL (LUFS), final programme IL/LRA and time scales. The concrete M18 proposal includes progressive waveform/M/S previews but deliberately stages rolling LRA and cumulative IL; final programme IL/LRA remain authoritative. Editing uses existing numeric section-gain drafts rather than freehand automation. Combined actual High Sierra/video/three-hour qualification follows M18 and must occur no later than M19; these candidate numbers do not authorize implementation.
3. **Separately bounded region suggestions and LRA-directed proposals.** Reuse full-resolution LRA distribution/tail diagnostics; group windows into editable investigation/processing regions, not claimed editorial scenes or automatically faulty passages. Tail membership alone does not prove causal contribution. Begin selective attenuation; selective raising requires separate qualification. Introduce neighbour/transition/preservation controls only under an approved method and remeasure LRA, IL and TP together. A true-peak limiter is supporting separately qualified work, not a mandatory prerequisite to early timeline interaction.

## Preservation and future lanes

Meet the selected ceiling with justified interventions, rather than minimise LRA or equalise scenes. Preserve dialogue clarity, ambience and dramatic contrast. Whole-mix gain cannot independently improve dialogue/background balance; component-relative dialogue/M&E processing and validated intelligibility analysis remain separate.

Plan DSP parameter/activity lanes and a change-overview lane highlighting large, rapid, sustained or dense alterations, each with a readable explanation and source-time location. Descriptive highlights do not become undeclared compliance criteria. Define freeze semantics before implementation: locking a processing setting differs from preserving a passage’s audible result, which also constrains global gain and overlapping processing. Replanning respects protections, explains conflicts and updates whole-programme verification.

M16 only implements approved proposal limits/inspection. It does not implement custom LRA rules, automatic analysis, waveform/timeline editing, region suggestions, freeze controls, new DSP or intelligibility methods. Keep original/applied audio identity and analysis-only references distinct.

## Historical applied-proposal explanations — future timeline requirement

On reopening, eventually restore reasons behind the exact applied proposal, its
relevant measured baseline and requirement/policy/method provenance as historical
information. Private candidate measurements and old verdicts must never become
current verification. Distinguish explanations for the still-applied plan from
those associated with subsequently edited processing; preserve that association
through immutable processing identity. M17 records this extension point only.
Keep the first interactive timeline/progressive-analysis implementation bounded;
stage nonessential editing/detail rather than delay useful selective LRA proposals.

The M18 proposal includes one versioned last-applied historical explanation record,
not an undo/history browser. It also plans the combined deferred qualification directly
after M18, retaining the no-later-than-M19 limit and subsequent main-integration review.
