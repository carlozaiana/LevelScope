#!/bin/sh
set -eu
R=/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation
B="$R/build/prores-integration-20261021"
ctest --test-dir "$B/build" --build-config Release --output-on-failure --output-junit "$B/ctest.xml" > "$B/ctest.log" 2>&1
python3 -B /private/tmp/levelscope-prores-integration-20261021/run-regressions.py sequoia
