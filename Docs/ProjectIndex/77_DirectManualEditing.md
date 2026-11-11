# Direct manual editing — 10 November 2026

Continuation: [Timeline selection and boundary editing](78_TimelineSelectionBoundaries.md).

Latest correction: [Re-measure after gain editing](../RemeasurePreviewCrash-2026-11-10.md).
The original candidate-02 remains preserved but is superseded for use by that fix.
Corrected app source: `5ac6ac57530004ab28167e65e63fda7475b09b68`; 2,154 affected
checks passed on each actual host. The crash was reproduced before the fix with
the same null-pointer paint failure as the user report.

See [handover](../DirectManualEditing-2026-11-10.md) for immediate section gains,
session-local disk-backed Undo/Redo, full-section marquee selection, reset/merge
semantics, provisional waveform feedback, default peak warnings and native menus.

All Current metrics remain unavailable after processing changes until explicit
measurement. Proposal acceptance and core DSP remain unchanged. Adjustable transitions
and custom key assignments follow separately. Publication/main integration awaits review.

Preserved `candidate-02` is built from `082d45ff7d88768f44bd9b557f16da561f2b7d93`.
Each actual host passed 2,625 focused checks: Sequoia 15.7.9 arm64 and High Sierra
10.13.6 x86_64. The handover contains the exact app paths, manual checklist and
acceptance receipt. Earlier packages and the running user app are preserved.
