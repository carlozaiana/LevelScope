# M19 qualification evidence — 5 November 2026

Current implementation: candidate-02, source
`cf275afb133118dd2a06e7bca59f82f98898b3cd`, branch `codex/m19-lra-passages`.
The approved impact-selected qualification is complete; ready for checkpoint review.
No M19 publication, tag or main merge is authorized or performed here.
See [the handover](../../../../Docs/ProjectIndex/68_M19QualificationAndHandover.md)
for scope, identities, limits, commands and readiness status.

- `candidate-01-identity.json`: original frozen implementation, all eleven direct
  automated groups on actual Sequoia and High Sierra.
- `candidate-02-identity.json`: UI legend correction, fresh native UI and measured
  workflow checks on both platforms. Original candidate-01 remains unchanged.
- `local-integrity-and-reuse.json`, `golem-integrity.json`: all 377/380 artifacts
  intact on both Macs, only the two legend UI files differ in product sources.
  Numerical reuse retains candidate-01's original source/build/binary identities.
- `sequoia-01`, `high-sierra-01`: exact explicit commands, executable hashes,
  platform, exit codes, raw reports, timings and peak RSS. Raw generated audio
  and disposable binary sessions are deliberately outside Git.
- `sequoia-native-02`, `high-sierra-native-02`: real native product components,
  controlled picker completions/events, device-closed inspection, assessment,
  cancellation, Apply, independent Current and history. Screenshots are fresh
  and visually reviewed. These are not claims of manual OS-picker interaction.
- `numerical-ui-review.json`, `compared-session-cores.json`: cross-platform
  numerical/processing identity comparisons and checksums of local sessions.
- `resources-*.json/log`, `supplemental-provenance.json`: explicit grouping plus
  envelope-selection resource probes. No three-hour audio was decoded.
- `*-surround-rates.json/log`, `surround-rates-comparison.json`: supplemental
  candidate-02 5.1 workflows at 44.1/96/192 kHz on both Macs; exact numerical,
  envelope/global-gain and policy/requirement comparisons. Original eleven-group
  evidence remains unchanged.
- `listening-*.json/log`: approved Saffire sequences, exact EOF/zero underruns,
  separate perceptual observations, supplied-file hashes, and the investigation
  of quiet clip halves. The final user inspection reports all four stereo and
  corresponding 5.1 captures flawless. Inspection hardware/method was not specified;
  full physical multichannel routing is not thereby qualified.
- `cinema-excerpt-provenance.json`: unchanged original identity and exact
  byte-identical eight-minute extraction. The cinema attempt fails its explicit
  15 LU ceiling honestly, without Apply or changed preferences.
- `current-status.json`: overall readiness, separate from automated completion.
  `review.py` deliberately does not qualify the candidate from numerical/native
  receipts alone. Coordinated listening observations and final review are separate.
- `completion-review.json`: all final required checks, candidate source/build
  identity and readiness decision. Reproduce read-only with
  `python3 Tests/Qualification/M19/complete_review.py build/m19-qualification-2026-11-05 --output /tmp/m19-completion-review.json`.
  `high-sierra-final-integrity.json` and `remote-final-refs.txt` retain the final
  package/ref verification. Frozen identity files retain their original
  pre-qualification flags; they were not rewritten after testing.

Development failures, the first unreadable legend image, and the sandbox `time`
wrapper error remain in the raw evidence; later receipts do not rewrite them.
The first completion-auditor receipt records an assertion-schema handling error
(the inherited measurement test reports a count, newer tests report named checks).
Only that read-only auditor was corrected; no product change or test rerun was needed.
Frozen apps/media are under the preserved local and Golem qualification folders.
Raw JUCE receipts retain their original CRLF and trailing spaces so their bytes
and recorded hashes remain unchanged. Whitespace checks apply to maintained
source/documentation; raw historical evidence is intentionally not normalized.
No installation, video execution, main merge, tag or M19 push is included here.
The 23 existing tag objects/targets remain unchanged.
