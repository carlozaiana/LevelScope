# Standalone Requirements (Deterministic Conforming Workflow)

## 1) Purpose and scope
- Define v1 requirements for a future standalone app that performs deterministic loudness-conforming workflows.
- Align standalone behavior with existing Core contracts and standards roadmap.
- Cover ingest, analysis, proposal, apply/audition, and export responsibilities.
- Exclude implementation details not yet frozen in current architecture docs.

## 2) Supported input/output formats and channel layouts

### Initial (v1)
- **Input audio:** TODO (exact container/codec list not yet defined).
- **Output audio:** TODO (exact container/codec list not yet defined).
- **Channel layouts:**
  - Non-ambisonic layouts up to 12 channels (7.1.4) should be supported, consistent with current scope assumptions.
  - Main output layout must match main input layout for v1.
  - LFE handling must remain explicitly policy-driven (detector/apply/exceptions).

### Future (post-v1)
- Ambisonics/object workflows: TODO.
- Expanded format matrix (broadcast interchange/delivery variants): TODO.
- Stem-aware scene ingest/export: TODO.

**Assumption:** Standalone I/O matrix should be finalized only after compliance/reporting output requirements are frozen.

## 3) Standalone workflow contract
Required pipeline contract:
1. **Ingest**
   - Load source audio and metadata required for processing/reporting.
   - Validate channel layout and project/profile settings before processing.
2. **Analyze**
   - Compute loudness/features using Core-compatible analysis semantics.
   - Record baseline analysis snapshot used for proposal decisions.
3. **Propose**
   - Generate conforming adjustment proposal against selected profile/target.
   - Proposal must be reproducible from saved inputs + settings.
4. **Apply/Audition**
   - Apply proposal through deterministic processing chain.
   - Allow audition pass without mutating original source.
5. **Export**
   - Render conformed output and required report artifacts.
   - Emit sufficient metadata to re-run and verify the same job.

## 4) Determinism requirements
- Same input media + same settings + same profile + same engine version + same platform + same floating-point mode must produce the same analysis/proposal/export results.
- Determinism applies to both:
  - Interactive/offline standalone run.
  - Batch job run.
- Offline and batch parity requirements:
  - No algorithmic divergence between interactive offline render path and batch render path.
  - Any unavoidable numeric tolerance must be documented as TODO with explicit thresholds.
- Versioning requirement:
  - Reports must include engine/profile version identifiers to support reproducibility checks.
- TODO: Define cross-platform tolerance policy (numeric acceptance criteria across CPU/OS/toolchain combinations).

## 5) Transport/analysis semantics in standalone context
- Standalone must define analysis timeline semantics independent of DAW transport quirks.
- Required semantics:
  - Explicit run boundaries (start/end) for analysis windowing.
  - Explicit behavior for seeks/restarts/reanalysis requests.
  - Stable frame/timing basis for reproducible metrics.
- Reanalysis semantics (v1):
  - Segment-level reanalysis is **allowed** for operator iteration.
  - Compliance pass/fail decisions must use full-program analysis unless a profile explicitly defines segment scope.
  - Segment reanalysis reports must record segment bounds and state that surrounding context may affect integrated metrics.
- **Assumption:** Existing plugin transport guard concepts inform behavior, but standalone should not depend on host callback freshness models.
- TODO: Define canonical context padding/window policy when segment-level reanalysis is run.
- **Assumption:** Existing plugin transport guard concepts inform behavior, but standalone should not depend on host callback freshness models.
- TODO: Define canonical handling of partial re-runs and segmented jobs.

## 6) Profile/compliance requirements linkage (EBU R128 / ATSC A/85)
- v1 must support selectable compliance profiles at minimum:
  - EBU R128
  - ATSC A/85
- For each profile, standalone requirements must map:
  - target metrics,
  - gating/measurement interpretation,
  - pass/fail criteria,
  - required report fields.
- Profile parameter policy must be split into three classes:
  - **Hard-locked by profile:** cannot be edited.
  - **Soft defaults:** editable with an explicit warning.
  - **Free parameters:** user-editable without profile-compliance warning.
- **Assumption:** Detailed requirement-to-control mapping will be maintained in a dedicated compliance matrix document.
- TODO: Finalize per-profile parameter classification and warning text policy.

## 7) Reporting/export artifact requirements
Each export job must record:
- Input identity:
  - content hash (SHA-256 or equivalent cryptographic hash),
  - channel layout metadata,
  - sample rate metadata,
  - source path/URI (informational, non-authoritative identity field).
- Output identity and render configuration.
- Selected compliance profile and target settings.
- Key analysis metrics (pre and post).
- Proposal summary (what adjustments were recommended/applied).
- Processing chain configuration relevant to reproducibility.
- Engine/build/profile version identifiers.
- Timestamp and job mode (interactive offline vs batch).
- Pass/fail compliance result.

Minimum required fields for reproducible pass/fail verification:
- Profile ID and profile version.
- Metric values pre and post processing.
- Thresholds/targets used for pass/fail evaluation.
- Processing graph definition and module versions.

Artifacts to emit:
- Conformed audio output.
- Machine-readable report: TODO (exact schema/format).
- Human-readable summary report: TODO (exact format).

### 7.1) Batch failure handling requirements
- Retry policy:
  - Deterministic retry count/order must be defined in job configuration.
  - Retries must not mutate analysis/proposal settings between attempts.
- Partial job behavior:
  - Each item result must be recorded independently (success/failure/skipped).
  - Batch summary must include counts and stable identifiers for failed items.
- Deterministic error codes:
  - Failures must emit stable, documented error codes suitable for automation.
  - Error reports must include stage-of-failure (ingest/analyze/propose/apply/export).
- TODO: Finalize canonical error code registry and retry backoff policy.

## 8) Non-goals for v1
- Real-time DAW-host transport integration behavior.
- Full automatic dialog intelligibility optimization loop.
- Stem + marker-track aware scene logic.
- Ambisonics/object-based workflows.
- Region-specific compliance profiles beyond EBU R128 and ATSC A/85.
- Cloud/distributed render orchestration.

## 9) Open questions / decisions needed
- TODO: Final v1 input/output format list (containers/codecs/metadata support).
- TODO: Canonical report schemas (machine-readable + human-readable).
- TODO: Numeric tolerance policy for determinism validation across platforms.
- TODO: Batch job manifest format and retry/error semantics.
- TODO: Final per-profile lock/default/free parameter tables.
- TODO: Profile override policy (what is locked vs editable per profile).
- TODO: Scope of reanalysis primitives (full-file only vs segment-level).