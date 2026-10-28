# After M7 — practical stereo/5.1 workflow

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

Historical M8 contract and post-M7 sequence below are retained. As of 28 October
2026, M8–M11 are published on main, M12 is a Sequoia-tested development checkpoint,
and combined M13 candidate-02 is qualified and authorized for checkpoint publication/main integration. The current staged order and optional-marker,
three-strategy design are in [47](47_M13BoundedSectionTrimProposals.md) and
[the shared specification](../ProposalEngineSpec_v2.md#22-staged-development-roadmap).
The bounded implementation follows the explicit impact plan in [48](48_M13ImpactQualificationPlan.md); broader strategies remain proposals.

M8 approved 25 October 2026, including the session, replaceable-cache and coherent-save clarifications below. Implementation and qualification are authorized; publication remains a separate checkpoint action. The user deferred
7.1.4 and further immersive expansion after M7. Preserve the extensible channel
interfaces and qualified 7.1 behavior, but optimize subsequent work for practical
stereo/5.1 broadcast use. Read with `ProposalEngineSpec_v2.md`; its broad future
engine design does not expand the existing versioned profile or imply a universal
LRA ceiling.

## Approved M8: reproducible sessions and safe reopening

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

Approved defaults: `.levelscope`, measurement caches included, exact conservative
build/platform compatibility. Budgets on each Mac: three-hour Source + Current
with 1,000 sections; save and core/cache parsing each <=5 seconds on local disk,
additional persistence memory <=256 MiB, UI interaction p95 <=100 ms, cancellation
acknowledgement <=1 second, document <=128 MiB. Full-media verification and video
indexing are separately timed, asynchronous, bounded and cancellable.

## M8 implementation contract

One programme is distinct from its physical files. Version 1 records a programme
ID, a file table (stable ID, locator, content SHA-256 and format description), and
ordered channel assignments (file ID, channel index, explicit role). M8 supports
only `interleaved` assembly with one file and identity channel ordering. Unsupported
assembly versions/modes fail safely; multiple-file importing is not implemented.
Programme identity uses content, channel mapping, alignment and range, not pathname.
M8 retains the qualified single-file measurement binding behind this representation;
future assembly support must version that binding and cannot reuse it by accident.

The `.levelscope` container has a bounded, independently checksummed core and
separately framed/checksummed replaceable measurement records. Invalid core data
fails atomically. Absent, truncated, damaged or incompatible caches are discarded
with an explanation when core framing is intact. They cannot prevent recovery of
valid settings, annotations and references. Keep original measurement provenance;
recompute evaluation from validated records. Missing media opens explicitly
unresolved; no measurement, audition or inherited pass is available until verified.

Save captures one immutable message-thread snapshot. Persistent document revision
covers settings, annotations, video references and accepted measurement changes;
audio revision remains separate. Save completion records only the captured revision:
new edits or measurement completions stay dirty. A persistent document UUID does
not replace the unique runtime opening/job identity. Old completions cannot install
results after New/Open, including reopening the same document.

Save uses an owned same-directory temporary and checked atomic replacement. It
must not overwrite media (including aliases), leave partials on cancellation, or
replace the previous document on create/write/flush/rename failure. No autosave,
media copying or automatic device opening. Open/New requires stopped transport;
installed sessions start device-closed and routes-muted. Measurements remain
independent of monitor settings and reference picture offsets.

## Following dependency order

1. **Separate-file stereo/5.1 import (next after M8).** Two mono files form stereo,
   or six mono files form 5.1. Explicit ordered roles, matching sample rates and
   frame counts, confirmed alignment, per-file verification, filenames beside roles,
   grouped save/reopen/relink and an explicitly incomplete group when any file is
   missing. No silent channel reassignment, padding, resampling or mix summing.
   This importer is not implemented by M8.

2. **Verified PCM export and an artifact-bound QC report.** Render explicit Current
   State through the shared full-precision path, pin output format/channel order,
   remeasure the actual written artifact and evaluate that artifact. Never export
   monitor gain/routing. Define full-scale/quantization/dither and file-size limits
   before implementation; write failures must not publish a misleading result.
   Keep report wording as narrow as M2 and identify exact exported hashes.
3. **Workflow reporting and handoff improvements.** A concise human-readable report
   can compare measured Source, Current and exported audio with diagnostic context,
   profile provenance and unresolved limitations. Reports require stable session
   and artifact identity; they must not imply complete broadcaster compliance.
4. **Container work driven by actual delivery media.** Add RF64/other large-file
   support when selected export rates/layouts exceed existing container limits,
   potentially as a prerequisite to a broad long-programme export claim.

Section processing, dynamics, automatic correction, additional video features,
7.1.4, Atmos and VST3 parity remain separate proposals. The existing optional paint
diagnostic executable remains an investigation aid for a later reproduction;
occasional accepted picture skips alone do not start another investigation cycle.

## Future component purposes and independent invalidation

Channel assembly places channels of one mix; it does not sum stems. Stem composition
would explicitly sum/route selected components, such as M&E and one language dialogue
stem. Analysis references assist analysis of an existing mix and contribute no audio
to playback, processing or export. Future component records need stable ID, purpose,
language/version, per-file identities, exact rational/source-frame alignment, channel
mapping and explicit audible/analysis-only participation. Unknown purposes remain
unsupported, never silently audible.

Audible programme identity includes every audible component's bytes, alignment,
selection, channel mapping and composition operation. Changes invalidate ordinary
programme measurements, processing verification, proposals and audition preparation.
Analysis-reference identity is separate: changes invalidate only dependent dialogue
analysis and its reports; unchanged audible programme measurements remain current.
Paths and relinking alone do not change either identity after bytes are verified.
Profile annotations remain separately bound to the applicable measured analysis.

Dialogue-stem loudness, speech-gated measurements of a final mix, and dialogue/
background diagnostics are distinct products with distinct input/method identities.
No speech detection, stem mixing or dialogue analysis is authorized in M8. Any
intelligibility score or delivery-specific dialogue-gated method needs its own
specification, primary-reference review and validation. These follow separate-file
import and verified export/reporting as bounded proposals. 7.1.4 remains deferred.
