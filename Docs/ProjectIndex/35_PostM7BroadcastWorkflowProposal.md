# After M7 — practical stereo/5.1 workflow

Proposal only; no implementation is authorized by this document. The user deferred
7.1.4 and further immersive expansion after M7. Preserve the extensible channel
interfaces and qualified 7.1 behavior, but optimize subsequent work for practical
stereo/5.1 broadcast use. Read with `ProposalEngineSpec_v2.md`; its broad future
engine design does not expand the existing versioned profile or imply a universal
LRA ceiling.

## Recommended M8: reproducible sessions and safe reopening

Save the work already supported by M1–M7 so a long programme can be resumed without
reconstructing its channel interpretation, global trim, sections, target and video
offset. This establishes a stable document and identity boundary for subsequent
verified exports and reports. It does not add processing.

Bounded scope:

- One versioned local session document, explicit Save / Save As / Open, atomic
  replacement through an owned temporary file, and clear unsaved-change state.
  No autosave, cloud sync, project-media copying or embedded credentials.
- Store selected source content identity and locator, ordered confirmed roles and
  provenance, programme range, full-precision global trim and processing identity,
  profile identity/definition digest and applicability confirmation, M5 section
  annotations, optional reference-video identity and exact rational signed offset.
- Preserve measured snapshots only with their exact method/build, source, role,
  range, processing and profile bindings. Reopened records begin untrusted. Verify
  the external media content and supported schema/method before marking any
  measurement current or reconstructing an M2 evaluation. Cached predictions never
  satisfy measured rules. Derived diagnostics and M3 proposals can be rebuilt from
  validated measurements rather than saved as independent truths.
- Missing, changed or relinked media has an explicit stopped unresolved state.
  Relinking must verify content; a same-named replacement cannot inherit a pass.
  Unsupported/future schemas or corrupt documents fail clearly and preserve the
  previous session. Old measurements remain unavailable/stale when exact reuse
  cannot be established; offer explicit remeasurement.
- Reopening never opens a device or starts playback. Monitoring starts muted;
  no machine-specific physical output assignment is silently restored. Video
  remains subordinate to selected audio, with embedded sound undecoded.

Dependencies: published M7 identities and limitations, stable serialization of
the existing in-memory models, current profile validation and qualified measurement
worker, atomic filesystem operations and bounded cancellable media verification.
No eight-output monitoring equipment is needed.

Acceptance tests on arm64 Sequoia, then the same candidate on actual High Sierra:

1. Round-trip stereo/5.1 sessions exactly: non-round gain, original role order,
   rational ranges/offsets, section IDs/labels and profile/applicability. Retain a
   supported 7.1 regression without broadening listening claims.
2. Unchanged media with matching identities permits measured reuse; changed bytes,
   roles, range, processing, method/build or profile mismatch never inherit an
   invalid verdict. Profile-only compatible reevaluation remains distinct from
   audio remeasurement. Predictions and diagnostic displays cannot authorize reuse.
3. Cross-Mac reopen and deliberate relink to the same content, including a moved
   file and a same-name different-content negative case. Missing video must not
   replace or reinterpret the authoritative audio.
4. Invalid/oversized documents, unknown versions, duplicate IDs, non-finite gain,
   invalid roles/ranges and malformed cache records fail safely. Inject disk-full,
   create/write/flush/replace failures and cancellation; retain the last saved file
   and previous valid session, and remove owned partials.
5. A representative three-hour snapshot remains bounded and responsive. Freeze
   specific open/save, memory and cancellation budgets before implementation,
   accounting separately for unavoidable full-media hash time.
6. Silent UI Save/Open/unsaved-change/cancel/relink and stopped lifecycle on both
   Macs; M1–M7 regressions and relevant H.264/ProRes identity/seek/offset checks.
   Coordinate any device use separately. No physical listening is needed to qualify
   session persistence.

Decisions to settle at approval: session filename/extension; whether retaining
validated measurement caches is in the first M8 candidate (recommended for long
programmes); exact resource bounds. Existing build identity is deliberately
conservative: an incompatible build invalidates cached measured reuse.

## Following dependency order

1. **Verified PCM export and an artifact-bound QC report.** Render explicit Current
   State through the shared full-precision path, pin output format/channel order,
   remeasure the actual written artifact and evaluate that artifact. Never export
   monitor gain/routing. Define full-scale/quantization/dither and file-size limits
   before implementation; write failures must not publish a misleading result.
   Keep report wording as narrow as M2 and identify exact exported hashes.
2. **Workflow reporting and handoff improvements.** A concise human-readable report
   can compare measured Source, Current and exported audio with diagnostic context,
   profile provenance and unresolved limitations. Reports require stable session
   and artifact identity; they must not imply complete broadcaster compliance.
3. **Import/container work driven by actual delivery media.** Add separate-file
   grouping with explicit ordering/alignment when needed. Add RF64/other large-file
   support when selected export rates/layouts exceed existing container limits,
   potentially as a prerequisite to a broad long-programme export claim.

Section processing, dynamics, automatic correction, additional video features,
7.1.4, Atmos and VST3 parity remain separate proposals. The existing optional paint
diagnostic executable remains an investigation aid for a later reproduction;
occasional accepted picture skips alone do not start another investigation cycle.
