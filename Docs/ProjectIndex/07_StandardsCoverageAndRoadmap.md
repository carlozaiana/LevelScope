# Standards Coverage & Roadmap Matrix

## Implemented/visible now
- Loudness tracking primitives aligned with BS.1770 style processing (momentary/short-term/integrated support in architecture and code paths).
- Running LRA visualization/metrics in current plugin workflow.
- True-peak protection path via limiter (with optional detector oversampling).

## Not yet fully realized (product-level)
- Full automatic conform/proposal engine for target matching with minimal intervention.
- Dialog intelligibility scoring/optimization loop.
- Standalone app ingest/analyze/apply/export pipeline.
- Stem + marker-track aware scene logic.

## Standards focus for next phases
1. EBU R128 compliance profile
2. ATSC A/85 compliance profile
3. Region/profile extensions (future)

## Recommended future deliverable
Add a machine-readable compliance matrix file (e.g., `Docs/ComplianceMatrix.md`) that maps each standard requirement to:
- implementation status,
- validating tests,
- user-visible controls,
- export/report outputs.