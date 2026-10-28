# M13 combined qualification evidence

**Candidate-02 qualified on actual arm64 Sequoia and Intel High Sierra, 28 October
2026. Ready for checkpoint/publication review; nothing has been published or merged.**
Main remains M11 `cfb1033aa7fb65a2655db276bee696161196bdeb`.

Source manifest: `50eba19a2d66096607ec878c8eaa05e610580c7f5e8854d7b375d846a89c5f27`.
Production build: `5fd55bb1a81c80b4b556715936b6873bccae29a5d04c00b0890a9b43296a20d2`.
Universal app executable: `6f843857838c64accc92536f6ccba2c9b0f36113dacca56dc2beef55c96233c5`.

See [the handover](../../../../Docs/ProjectIndex/49_M13QualificationAndHandover.md),
[completion gate](qualification-complete.json),
[cross-platform comparison](cross-platform-comparison.json) and
[current status](status.json). The frozen identity files retain their original
pre-qualification status fields; the separate completion gate records acceptance.

Each actual Mac passed the 22 required short/reference checks, 40 M13 cases,
32 M12 processing cases, 17 silent UI checks, six foreign-session exchanges,
80 injected video receipts plus timing regressions, four coordinated muted physical
video cases, and one consolidated three-hour composition/section/RF64 workflow.
The long phases took 117.3 minutes on Sequoia and 202.1 minutes on High Sierra,
including their independent decoder checks, and passed every resource budget.
Both 9,331,200,104-byte RF64 artifacts have SHA-256
`953b2bcecbf068f5fe153461fb1ea0999939ddb40adf607052b05f4ad2146ade`.
The final comparison found zero float32 ULP difference in the retained captures
and maximum metric difference `3.552713678800501e-15`, within the documented
comparison-only tolerances. Profile evaluation bounds remain unchanged.

Candidate-01 and candidate-02 remain preserved under
`build/m13-qualification-2026-10-28`. Candidate-02 repairs only five qualification-tool
Python files; every production source and C++ test/standalone binary is identical.
Both source manifests and the explicit equivalence record are retained. Candidate-01
Sequoia short/UI/injected-video/cinema/long evidence keeps its original identity.
Fresh candidate-02 evidence includes both repaired-tool sets, both physical video
sets, both session exchanges and all High Sierra phases. No historical long runner
was invoked recursively. The final gate requires this exact reuse relationship.

The first physical-runner preflight failed on a composition report field before any
device was opened. That failure remains under `historical`; it is not a passed video
run. The fixed runner validates the verified programme-definition digest. The cinema
qualification completed safely with an honest policy-limited failure: one whole-programme
section exceeded the 20% local-duration budget. No plan was applied, exported or played.
Its original media remains unchanged; a payload-identical RF64 rewrap is documented.
User-reported Nuendo observations remain external context, with the Cubase 10 screenshot
separate. They are not exact reference vectors or proof of an engine discrepancy.

No binaries, raw audio/video, raw captures, credentials or connection setup are archived
here. Reproduction and the final gate use the preserved local candidate/results trees;
this reduced Git archive does not replace those test packages. `evidence-manifest.json`
hashes every archived file except itself. Local and High Sierra post-test integrity
records confirm preserved binaries; all 16 earlier tag objects/targets and all branch
and worktree heads are unchanged. Main remains at M11 pending separate approval.

Physical video evidence is distinct from injected timing tests. All eight fresh muted
physical cases recorded zero skips. Occasional picture skips remain an accepted
limitation only under the retained audio continuity, authoritative-clock drift, identity,
memory and lifecycle conditions. Earlier strict failures remain historical failures.
Full physical 7.1 listening is separately pending and outside M13. Fixed policy/manual
sections, no limiter/dynamics and the inherited M12 report caption remain documented
limitations; the exact report build identifies M13.
