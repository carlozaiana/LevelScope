# LevelScope — Proposal Engine Specification v2

Status: Proposed  
Scope: Standalone-first, reusable in Core, plugin-preview-compatible later  
Depends on: `Docs/ArchitectureSnapshot.md`, `Docs/ComplianceMatrix.md`, `Docs/ProjectIndex/*`

---

## 1. Purpose

The Proposal Engine generates explainable, minimally destructive processing plans that move a source mix toward a selected broadcaster target profile.

Primary success condition:

1. final verified output must be inside the target compliance window
2. inside that compliance window, the engine should minimize audible alteration
3. the user must be able to inspect, override, freeze, and re-measure the current state before export

This engine is intended first for the future standalone workflow, while keeping the implementation reusable in `LevelScopeCore`.

---

## 2. Product Philosophy

The engine is not a black box. It should answer:

- what is wrong
- where the problem comes from
- which module is chosen to solve it
- how much processing is proposed
- what trade-off mode is active
- what changed after user edits

Core philosophy:

- meet the broadcaster requirements
- preserve cinema intent as much as possible
- prefer the least destructive intervention first
- solve LRA by tail steering, not flattening everything
- support user-in-the-loop refinement

---

## 3. Compliance Model

## 3.1 Pass conditions

The engine succeeds only if final verified output is inside the profile limits.

Typical default interpretation:

- Integrated Loudness: pass if within ±0.5 LU of target
- LRA: pass if not above target, unless the selected profile explicitly allows overage
- Max Peak / True Peak: pass if at or below target ceiling

Profile-specific rules must be data-driven and owned by the compliance/profile layer, not hard-coded inside proposal methods.

## 3.2 Constraint types

### Hard constraints
- true peak / max peak limits
- technical/legal broadcaster limits

### Near-hard constraints
- integrated loudness target window

### Soft constraints
- exact LRA hit when lower-than-target is acceptable
- intelligibility target details
- tonal preservation
- scene preservation
- automation smoothness
- limiter activity budget

---

## 4. Displayed States

The UI should expose three user-facing value groups:

### A1. Source
Measured values of the original/imported source program.

### A2. Target
The selected profile target values and allowed tolerances.

### A3. Current State
Measured values of the current editable state after proposal, user edits, freezes, and re-measurement.

Current State is the operational truth for export readiness.

## 4.1 Internal distinction
Internally, the system must distinguish:

- Predicted Proposal State
- Measured Current State

The UI may choose to show only Source / Target / Current State, but the engine must keep predicted and measured states separate.

## 4.2 Re-measurement rule
If the user changes parameters or freezes regions/scenes, the app must support re-measuring the Current State before generating a new proposal.

This may be:
- automatic after edit completion
- user-triggered via explicit “Re-measure Current State”
- or both, depending on workflow cost

A new proposal must start from the actual measured Current State, not stale predicted values.

---

## 5. Gate Visibility

The app should expose gate values numerically and visually.

Recommended visible gates:

- Integrated Loudness absolute gate
- Integrated Loudness relative gate
- LRA effective lower inclusion threshold / gate

## 5.1 UI representation
On the loudness timeline:
- show gates as horizontal lines

On histogram/distribution views:
- show gates as vertical markers on the loudness axis

Purpose:
- explain what content is counted
- explain why some sections contribute to IL and LRA and others do not
- increase trust and reduce trial-and-error

---

## 6. Time Scales

The engine must reason at multiple time scales.

### Program scale
Used for:
- integrated loudness
- target profile strategy
- global compliance planning

### Scene scale
Used for:
- preserving dramatic arc
- scene classification
- scene trims
- threshold adaptation at coarse scale

### Local window scale
Used for:
- short-term loudness shaping
- LRA tail management
- local dialog support
- compressor activity planning

### Event scale
Used for:
- transient peaks
- isolated true-peak violations
- short intelligibility risks
- limiter cleanup planning

---

## 7. Intervention Hierarchy

The engine should prefer the least destructive tool first.

Recommended intervention order:

1. Global trim
2. Scene/section trim
3. Broadband downward control for upper-tail loudness
4. Broadband upward support for lower-tail loudness
5. Spectral/dialog-focused support only where broadband would cause collateral damage
6. Limiter for residual peak containment only

This hierarchy must be encoded in policy and proposal logic.

---

## 8. Architecture

## 8.1 Shared foundation vs separate methods

All proposal methods must share the same:

- compliance/measurement truth
- analysis outputs
- attribution maps
- policy definitions
- current-state model
- freeze model
- verification pipeline

Alternative proposal methods should live in separate files/classes and operate on shared data structures.

This allows:
- apples-to-apples comparison
- listening evaluation between methods
- experimentation without duplicating compliance logic

---

## 8.2 Core modules

### `ProposalEngineFacade`
Orchestrates the full proposal workflow.

Responsibilities:
- receives Source + Target + Current State
- selects method + mode
- runs diagnosis, proposal, simulation, refinement, verification
- returns `ProposalReport`

### `MeasurementEngine`
Authoritative compliance measurement.

Responsibilities:
- integrated loudness
- short-term loudness
- momentary loudness
- loudness range
- true peak / max peak
- later: dialog/intelligibility proxies

Hard requirement:
- final verification must use this same truth model

### `AnalysisEngine`
Builds the descriptive program model.

Responsibilities:
- loudness curves
- running / rolling stats
- gated distribution views
- histogram and CDF
- peak map
- segmentation
- scene classification
- headroom and dynamic profile

### `AttributionEngine`
Computes metric contribution maps.

Responsibilities:
- IL contribution map
- LRA lower-tail / upper-tail contribution maps
- peak contribution map
- later intelligibility contribution map

### `PolicyEngine`
Owns user-facing proposal modes and weighting strategy.

Responsibilities:
- mode selection
- hard vs soft constraint weighting
- dialog protection strength
- preservation weighting
- limiter aversion
- scene policy weighting

### `ProxySimulationEngine`
Fast predicted-output simulator.

Responsibilities:
- simulate metric effects without full render
- support fixed-point / iterative LRA gate handling
- estimate candidate proposal quality quickly

### `RefinementEngine`
Performs bounded local refinement around an initial proposal.

Responsibilities:
- threshold nudging
- scene-trim balancing
- upward/downward rebalance
- limiter load reduction
- micro-correction planning

### `VerificationEngine`
Runs exact post-processing measurement after actual render/current-state measurement.

Responsibilities:
- measure final truth
- compare against profile
- emit pass/fail
- allow micro-correction if policy permits

### `FreezeEngine`
Maintains user-fixed constraints.

Responsibilities:
- freeze thresholds
- freeze regions
- freeze scenes
- freeze module choices
- protect frozen state during recalculation
- smooth transitions near freeze boundaries

### `ExplanationEngine`
Turns the result into a user-readable rationale.

Responsibilities:
- explain what was wrong
- explain where issues were found
- explain which modules were chosen
- explain what changed
- explain remaining risk or trade-off

---

## 8.3 Data structures

### `ProgramAnalysis`
Whole-program measured and derived state.

Fields:
- source metrics
- current metrics
- target profile
- loudness curves
- gates
- histogram
- CDF
- peak map
- scene list
- global headroom
- running stats

### `SceneAnalysis`
Per-scene or per-segment model.

Fields:
- time range
- scene class
- median / mean loudness descriptors
- local range
- peak density
- dialog confidence
- preservation weight
- protected/user-marked flags

### `ContributionMap`
Metric attribution structure.

Fields:
- metric type
- ranked contributing windows
- lower-tail vs upper-tail tags for LRA
- efficiency score
- alteration-cost estimate

### `ProposalIntent`
High-level intent layer.

Intent lanes:
- global loudness correction intent
- upper-tail range control intent
- lower-tail support intent
- peak safety intent
- intelligibility support intent

### `ModulePlan`
Derived low-level processing plan.

Fields:
- thresholds T0–T3
- global trim
- scene trims
- upward settings
- downward settings
- limiter settings
- optional automation
- final micro-correction plan

### `ProposalReport`
Top-level returned proposal result.

Fields:
- diagnosis summary
- predicted metrics
- measured metrics if available
- module plan
- explanation
- pass/fail status
- warnings
- next recommended action

---

## 9. Workflow States

Recommended top-level workflow state machine:

### `Idle`
No source analyzed yet.

### `SourceMeasured`
Source metrics available.

### `TargetSelected`
Profile selected and validated.

### `CurrentStateMeasured`
Current editable state measured and available.

### `Diagnosed`
Root-cause diagnosis complete.

### `ProposalDrafted`
Initial deterministic proposal available.

### `ProposalPredicted`
Fast proxy prediction complete.

### `ProposalRefined`
Bounded local refinement complete.

### `UserReview`
User inspects curves, gates, thresholds, and explanation.

### `CurrentStateEdited`
User changed something manually.

### `CurrentStateFrozen`
User froze one or more regions / scenes / parameters.

### `CurrentStateRemeasured`
Edited current state re-measured.

### `ProposalRecalculated`
New proposal generated around frozen constraints.

### `Verified`
Exact measured output available.

### `ExportReady`
Current state passes compliance window and is user-approved.

### `Exported`
Final export completed.

Invalid transitions should be blocked or clearly explained.

---

## 10. User-Facing Modes vs Internal Methods

These must remain separate concepts.

## 10.1 User-facing modes
These are policy presets:

- `StrictCompliance`
- `Balanced`
- `MaximumPreservation`

These modify:
- penalty weights
- aggressiveness
- limiter aversion
- scene protection
- allowable softness on non-hard metrics

## 10.2 Internal methods
These are implementation strategies:

- `RuleBasedMethod`
- `DistributionReshapingMethod`
- `HybridMethod`
- later `LocalOptimizerMethod`
- later `MlAssistedInitializer`

Methods are not the same as user modes.

---

## 11. Measurement and Analysis Outputs

Required measurement outputs:

- integrated loudness
- short-term loudness
- momentary loudness
- LRA
- true peak / max peak

Required analysis outputs:

- loudness-over-time curves
- histogram of short-term loudness
- CDF of short-term loudness
- current P10 / P95
- effective gates
- cumulative / rolling LRA
- peak map
- scene boundaries
- scene classification
- headroom map

Note:
- histogram is best for human understanding
- CDF is best for percentile math, threshold derivation, and LRA steering

---

## 12. Threshold Logic (T0–T3)

Thresholds are central and must be explainable.

### T1–T2: Preserve zone
This should bracket the healthy middle of the mix:
- normal dialogue
- stable narrative passages
- material best left alone

### T0–T1: Lower-tail support zone
This should catch:
- meaningful low-level content
- weak dialog passages
- relevant quiet material

It must avoid:
- silence
- dead air
- irrelevant room tone

### T2–T3: Upper-tail control zone
This should catch:
- upper-tail LRA contributors
- sustained loud sections
- likely TP-risk passages after correction

### T3
Marks stronger control / hard-control / limiter-risk territory.

### Automation rule
Thresholds may adapt at scene scale, but should not wobble rapidly at fine granularity.

---

## 13. Deterministic Expert Proposal

The initial proposal should be deterministic.

Meaning:
- same input + same settings => same result
- fixed reasoning chain
- explainable and debuggable

Recommended deterministic flow:

1. measure source and current state
2. compute compliance deltas
3. test whether global trim can solve IL safely
4. diagnose whether LRA excess is lower-tail, upper-tail, or split
5. place thresholds T0–T3
6. choose least destructive intervention by hierarchy
7. draft `ProposalIntent`
8. derive `ModulePlan`
9. run proxy simulation
10. refine only if required

The deterministic proposal is the anchor. Optimization is allowed only as refinement around this anchor.

---

## 14. Proposal Methods

## 14.1 `RuleBasedMethod`
Purpose:
- instant, explainable first proposal

Strengths:
- deterministic
- easy to debug
- fast
- good baseline

## 14.2 `DistributionReshapingMethod`
Purpose:
- improve LRA handling using loudness-distribution logic

Strengths:
- good for percentile steering
- good for threshold derivation
- good for analytical tail decisions

## 14.3 `HybridMethod`
Recommended primary method.

Composition:
- rule-based diagnosis and first draft
- distribution-aware LRA refinement
- proxy simulation
- bounded local refinement
- exact verification

This should be the default production method.

## 14.4 Later methods
- `LocalOptimizerMethod`
- `MlAssistedInitializer`

These should remain optional and isolated from the shared truth model.

---

## 15. Attribution Maps

The engine should compute at least:

### `ILContributionMap`
Shows which windows materially affect Integrated Loudness.

### `LRALowerTailContributionMap`
Shows which windows drive the lower tail.

### `LRAUpperTailContributionMap`
Shows which windows drive the upper tail.

### `PeakContributionMap`
Shows which windows create TP / peak violations.

### Later: `IntelligibilityContributionMap`
Shows which speech-bearing windows have support/masking risk.

This is the core of explainability and efficient correction.

---

## 16. Decision Strategy

## 16.1 Integrated Loudness
Preferred solving order:
1. global trim
2. scene trim
3. dynamic modules only if needed because IL conflicts with LRA/TP/intelligibility

## 16.2 True Peak / Max Peak
Preferred solving order:
1. preserve headroom during planning
2. local downward control for sustained problems
3. limiter for residual peaks

## 16.3 Loudness Range
Solve by tail steering:
- upper-tail issue => reduce upper-tail passages
- lower-tail issue => support meaningful lower-tail passages
- both => split strategy

Do not flatten the whole mix if localized tail steering can solve it.

## 16.4 Dialog intelligibility (later)
Preferred solving order:
1. dialog-aware levelling
2. center/dialog-specific support
3. spectral support
4. broad gain only when justified

---

## 17. Penalty Model

The engine needs explicit alteration-cost budgets.

## 17.1 Base penalties
- total gain movement
- total compression activity
- limiter activity
- spectral deviation
- automation roughness
- dialog disturbance

## 17.2 Content-aware penalties
- stronger penalty in quiet/simple material
- stronger penalty in protected scenes
- stronger penalty in dialog scenes
- lower penalty in dense/loud passages where processing is less audible

## 17.3 Practical tuning model
For each intervention type, define:
- preferred range
- caution range
- excessive range

These are soft budgets for tuning, not hard bans, except where compliance requires hard limits.

---

## 18. Scene Intelligence

Automatic scene recognition is useful but incomplete.

The engine should support both:

- automatic scene classification
- user-marked protected/key scenes

User protection marks should override automatic guesses.

Suggested flags:
- protect dramatic impact
- do not flatten
- do not lift ambience
- preserve dialog perspective
- preserve music phrasing

---

## 19. Final Correction Stage

For v1, do not add a second full leveller at the end of the chain.

Instead, implement a small `FinalCorrectionStage` before the limiter.

Responsibilities:
- tiny final global trim
- later optional slow section trim if needed
- feed final limiter
- re-measure exact result

Reason:
- simpler
- less risk of masking upstream design mistakes
- better for keeping the chain understandable

---

## 20. Standalone-First Strategy

This workflow should be implemented first in the standalone application.

Reasons:
- whole-program analysis
- offline iterative verification
- current-state re-measurement
- freeze-aware recalculation
- export readiness workflow

The plugin may later reuse:
- measurement pieces
- proposal preview logic
- explanation snippets
- selected module planning components

But the full proposal/export workflow belongs first in the standalone product.

---

## 21. Edge Cases

The engine must explicitly handle:

1. small IL miss with available headroom
2. IL correction that would create TP violation
3. LRA excess driven almost entirely by upper tail
4. LRA excess driven almost entirely by lower tail
5. both tails extreme
6. isolated outlier event such as one explosion
7. mostly quiet program with a few loud events
8. sparse ambience / room tone regions
9. music-driven sections where speech logic would be harmful
10. dialog masking without overall loudness deficit
11. upward support pushing content into later downward zone
12. proxy passes but exact verification misses slightly
13. frozen user regions that conflict with ideal compliance
14. no acceptable solution inside alteration budget

In the last case, the engine should be able to present alternatives, for example:
- stricter hit with more alteration
- better preservation with near-miss or reduced ambition on soft constraints

---

## 22. Recommended First Release

### v1
- exact measurement engine
- source/target/current-state model
- segmentation
- contribution maps
- deterministic rule-based proposal
- hybrid default method scaffold
- fast proxy simulation
- freeze-aware re-measurement and recalculation
- final exact verification
- explanation output

### v2
- stronger percentile steering for LRA
- better passage efficiency ranking
- scene-policy weighting refinement

### v3
- dialog-aware support
- intelligibility contribution map
- spectral/dialog-specific proposal lanes

### v4
- bounded local optimizer
- optional ML-assisted initializer

---

## 23. Recommended File Placement

Recommended documentation placement:

- main spec file: `Docs/ProposalEngineSpec_v2.md`
- add an index entry in `Docs/ProjectIndex/README.md`
- add a summary line in `Docs/ArchitectureSnapshot.md` under planned proposal/policy hooks
- later add code-level design notes under `Docs/ProjectIndex/` if the implementation splits into multiple tracked milestones

Reason:
- this keeps one canonical proposal-engine spec at the top `Docs` level
- keeps `ProjectIndex` as the navigation/index layer
- keeps `ArchitectureSnapshot.md` as the contract summary, not the full design document

---

## 24. Recommended Code Namespace / Folder Direction

Suggested future code grouping in Core:

- `Source/Core/Proposal/`
  - `ProposalEngineFacade.h/.cpp`
  - `MeasurementEngineBridge.h/.cpp`
  - `AnalysisEngine.h/.cpp`
  - `AttributionEngine.h/.cpp`
  - `PolicyEngine.h/.cpp`
  - `ProxySimulationEngine.h/.cpp`
  - `RefinementEngine.h/.cpp`
  - `VerificationEngine.h/.cpp`
  - `FreezeEngine.h/.cpp`
  - `ExplanationEngine.h/.cpp`

Suggested separate methods:

- `Source/Core/Proposal/Methods/RuleBasedMethod.h/.cpp`
- `Source/Core/Proposal/Methods/DistributionReshapingMethod.h/.cpp`
- `Source/Core/Proposal/Methods/HybridMethod.h/.cpp`

Suggested shared model types:

- `Source/Core/Proposal/Model/ProgramAnalysis.h`
- `Source/Core/Proposal/Model/SceneAnalysis.h`
- `Source/Core/Proposal/Model/ContributionMap.h`
- `Source/Core/Proposal/Model/ProposalIntent.h`
- `Source/Core/Proposal/Model/ModulePlan.h`
- `Source/Core/Proposal/Model/ProposalReport.h`

This keeps methods isolated while preserving one shared truth model.

---

## 25. Summary

The recommended v2 direction is:

An offline, explainable, freeze-aware hybrid proposal engine that uses exact measurement, attribution maps, deterministic first-pass rules, fast proxy prediction, bounded local refinement, and exact verification to meet compliance targets with minimum audible alteration.

This should be implemented standalone-first, with reusable Core components and method isolation for listening-based comparison.
