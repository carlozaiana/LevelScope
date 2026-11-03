# M18 first implementation slice — development validation

This is short Sequoia/arm64 development evidence, not full M18 qualification.
The logic handover `9935ad725670958971c1aaed28d8e84f9b3d3e9a` was cherry-picked
as `cdac0c3` before these changes. Main and the other worktrees were not changed.

`candidate-01/summary.json` accounts for all nine required runs. All passed
against measurement build
`d2b478709ba907c2eede7d9d4b6360205e2fa4f84b1c4f3d71d3104b6ea68b4b`.
`source-manifest.json`, `validation-manifest.json` and `binaries.json` retain
source/test inputs and binary hashes. All nine binaries have arm64 and x86_64
slices; every Intel slice declares macOS 10.13 minimum. CMake configuration and
build logs are retained. JUCE emitted its existing dependency warnings.

The working directory was `/Users/fuzz/.codex/worktrees/a411/LevelScope` and the
build directory was `build/m18-dev-01`. Configuration used Release, JUCE 7.0.12
from the existing dependency checkout, arm64/x86_64 and deployment target 10.13.
Only the standalone and explicitly selected test targets were built. The plugin
and default/all targets were not built.

The following commands use executables under
`build/m18-dev-01/<Target>_artefacts/Release/<Target>`. Outputs are in
`candidate-01`. Fixture roots were new owned directories under `/tmp`.

| Target | Arguments |
| --- | --- |
| LevelScopeMeasurementTests | `--report <measurement.json>` |
| LevelScopeComplianceTests | `--report <compliance.json>` |
| LevelScopeProposalTests | `--report <proposals.json>` |
| LevelScopeExportTests | `--artifacts /tmp/levelscope-m18-export-final-01 --output <export.json>` |
| LevelScopeSectionProposalTests | `--smoke --root /tmp/levelscope-m18-sections-final-01 --output <section-proposals.json>` |
| LevelScopeSectionProposalTests | `--preferences --root /tmp/levelscope-m18-preferences-final-01 --output <preferences.json>` |
| LevelScopeJobRequirementsTests | `--root /tmp/levelscope-m18-targets-final-01 --output <targets.json>` |
| LevelScopeSessionTests | `--short --output <sessions.json>` |
| LevelScopeDiagnosticsTests | `--report <diagnostics.json>` |

JUCE's macOS temporary files live in `~/Library/Caches/<executable>`, so execution
needed the sandbox escalation for those generated fixtures. No device was opened.

`iteration-01` retains initial failures and superseded reports. The initial
sandboxed runs could not create their JUCE cache files. Development testing then
caught the programme/core schema mismatch during the schema-8 change and an old
proposal test that assumed Current stayed unmeasured after Source. The format was
corrected and the test now exercises automatic initial Current plus an explicit
reset to unmeasured Current. Only `candidate-01` is final evidence for this slice.

The native window was launched and visually inspected with devices closed. A
label-encoding issue was corrected. Computer-use actions did not reliably operate
JUCE child controls, so those attempts do not count as interactive native
qualification. Component tests cover direct field edits, validation, unchanged
edits and compressed fixed-font geometry. Full native interaction, progressive
PCM, gain/inspector integration, saved explanation history, actual High Sierra,
cross-Mac, video and three-hour checks remain pending.
