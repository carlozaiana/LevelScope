#!/bin/sh
# Run from the isolated bundle root. No build, installation, GUI or audio playback.
set -eu
if [ "$#" -ne 3 ]; then
    echo 'Usage: sh Tests/run_m1_bundle.sh /absolute/python3 high-sierra|sequoia results/fresh-run'
    exit 2
fi
m1_python="$1"
m1_platform="$2"
m1_results="$3"
case "$m1_platform" in high-sierra|sequoia) ;; *) echo 'Unknown platform'; exit 2 ;; esac
if [ -e "$m1_results" ]; then echo 'Results path exists; choose a fresh run directory'; exit 2; fi
mkdir -p "$m1_results"
"$m1_python" -B Tests/collect_m1_environment.py --output "$m1_results/environment.json" \
    --identity identity.json --measurement bin/LevelScopeMeasurementTests \
    --scaffold bin/LevelScopeScaffoldTests \
    --app 'LevelScope Standalone.app/Contents/MacOS/LevelScope Standalone' \
    --require-os "$m1_platform"
run_check() {
    m1_log="$1"
    shift
    if "$@" > "$m1_results/$m1_log" 2>&1; then
        cat "$m1_results/$m1_log"
    else
        cat "$m1_results/$m1_log"
        exit 1
    fi
}
run_check measurement.log bin/LevelScopeMeasurementTests
run_check scaffold.log bin/LevelScopeScaffoldTests
run_check runner.log "$m1_python" -B Tests/EbuReferenceRunnerTests.py
run_check itu.log "$m1_python" -B Tests/run_reference_measurements.py \
    --binary bin/LevelScopeMeasurementTests --vectors vectors/itu \
    --report "$m1_results/itu.json" --compare baseline/itu.json
run_check ebu.log "$m1_python" -B Tests/run_ebu_measurements.py \
    --binary bin/LevelScopeMeasurementTests --vectors vectors/ebu \
    --report "$m1_results/ebu.json" --compare baseline/ebu.json
echo 'Hardware-independent M1 qualification and baseline comparisons PASS. UI qualification is separate.'
