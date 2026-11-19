# M20 UI and independent numerical evidence — 6 November 2026

**Non-audible acceptance complete; user listening review pending.** No M20
publication, checkpoint tag or main merge is authorized. No physical device was
opened. The user will audition the two independently verified music comparison
sessions personally. `acceptance-04.json` deliberately leaves combined M20
qualification false.

Accepted development app:
`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/m20-ui-integration/candidate-04/LevelScope Standalone.app`.
Product/test input commit: `18d87f9a3246ccbfa0c98963119aaff2e30f25aa`.
Build: `2bc868fb5198536720a2fa41c9b359c7c931753e25e32e66bc2429c69072ee18`.
App SHA-256: `08da73ad862f7beed937ade0df05c126191d01718aa4eff4708cf36ab2312ecf`.

## Evidence and provenance

- `identities/`: the four preserved package identities. All source/test/configuration
  and binary hashes are retained; packages themselves remain outside Git.
- `sequoia-numerical/` and `high-sierra-numerical/`: original candidate-01 receipts for
  19 selected affected groups, 4,248/4,247 checks respectively. The difference is an
  arm64-only cached-planning assertion; the audit also enforces Intel's 20 s target.
  Per-case JSON files are byte-exact aliases named by case. Original summary paths
  still refer to the owned original run directories, not these aliases.
- `sequoia-native/`, `high-sierra-native/`: fresh candidate-04 UI controls, actual
  native peers, all-attempt review, diagnostic/stale guards, dense geometry/history,
  explicit Apply/independent Current, save/reopen and cancellation: 93 checks each.
- `*-exchange.json`: fresh candidate-04 stereo/5.1 exchange in both directions,
  36 checks per receiving Mac. Exact historical attachment/settings/roles/processing
  survive; foreign caches do not confer verification. Independent Current passes.
- `holdouts/`: 87 checks and all six complete M19/M20 outcomes on the three original
  user-supplied inputs. The protocol in `Docs/M20IndependentReviewProtocol-2026-11-06.md`
  was committed before measurement. The music case passes both methods with exact
  independent Current reproduction; dialogue/ambience and 5.1 have no applicable
  candidate under the fixed limits. This is not a delivery or listening-quality claim.
- `development-failures/`: original failed row-notification/history-reader/native
  performance observations, preserved as failures. None is reclassified as a pass.
- `acceptance-04.json`: explicit completeness, artifact/source equivalence, unchanged
  numerical behavior, actual platforms, timing/storage and original-input guards.
  `candidate-02-reuse-proof.json` is intermediate history; the final acceptance file
  lists the final allowed product differences from candidate 01.
- `preserved-packages.json`: unchanged prior UI convenience and Logic package hashes.
  `holdout-inputs-before.json`: immutable original input identities/role order.
- `evidence-sha256.json`: byte-preserving inventory of the copied evidence.

Only five product files change after candidate 01: the M20 review view header/source,
the shared page-view header/source (old default unchanged), and the main window's
historical-page constructor. They retain immutable report views and select smaller
M20 pages. Numerical/measurement/policy/processing/session/device sources are otherwise
identical. Candidate-01 numerical/convenience evidence is explicitly reused;
candidate 04 receives fresh actual-platform native and exchange checks. The complete
acceptance guard is `Tests/Qualification/M20/review_ui.py`, run against the preserved
owned `build/m20-ui-integration` tree; missing evidence cannot produce acceptance.

No generated audio is committed here. Original media, all frozen packages, main,
all 25 tags and other working copies are preserved. See
[the full handover](../../../../Docs/M20UIIntegration-2026-11-06.md) and
[ProjectIndex 73](../../../../Docs/ProjectIndex/73_M20UIIntegrationAndIndependentReview.md).
