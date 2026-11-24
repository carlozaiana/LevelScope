# Graph boundary candidate 01

The acceptance receipt records two affected groups on actual Sequoia arm64.
High Sierra's new graph run is pending; the user deferred reconnection. No
physical device or audible test is authorized or required by these commands.

Preserved package:
`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/graph-region-boundaries-2026-11-24/candidate-01`

The package's own frozen runner verifies every artifact before and after testing:

```sh
python3 ABS_CANDIDATE/source/Tests/Qualification/M20/run_ui.py ABS_CANDIDATE ABS_FRESH_RESULTS --case graph-boundaries --case graph-ui --inner-session ABS_INPUTS/supplied.levelscope
```

For Golem, transfer the identical candidate and the byte-original input pair
from `build/applied-region-boundaries-2026-11-22/inputs/`: `supplied.levelscope`
and `supplied-audio.wav`. Keep those filenames together; the native test copies
and relinks its own session/audio before editing. Preserve macOS 10.13-compatible
archive metadata, as documented by earlier qualification. Run with normal user
cache access in the logged-in native desktop context. No broad suites or timing
benchmarks are requested.

`freeze.py` and `verify_preservation.py` here are exact archived copies of the
scripts used in the task's build directory, where their path-relative roots
resolve. They are provenance, not scripts to execute from this evidence folder.
The candidate/source runner above is the reproducible test entry point.

The failed initial development receipt is retained. Its graph-mode fixture
wait was corrected for accepted background preparation; the final product
Source is the same in both successful development groups and the frozen run.
