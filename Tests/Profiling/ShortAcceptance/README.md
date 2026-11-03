# Repeating the short M18 acceptance pass

Use an existing enabled Release build on the development Mac. The harness links the application objects and uses its original timer and controls. Only the diagnostic driver is compiled with access to private controls. File selections are returned to the production chooser completions. The injected output uses the existing test device, implements no physical I/O, and is first ticked as soon as the app starts its callback. Per-run role memory, audio and exported files are isolated from user projects. The native window and verified-copy storage require the same local execution privileges as the standalone.

Run from the repository root with a new working directory:

```sh
work=$(mktemp -d /private/tmp/levelscope-short-acceptance.XXXXXX)
python3 Tests/Profiling/ShortAcceptance/fixtures.py "$work/audio"
python3 Tests/Profiling/ShortAcceptance/build.py build/m18-peak-dev "$work/bin"
"$work/bin/short-acceptance" "$work/audio/programme-2.wav" "$work/stereo" 2
python3 Tests/Profiling/ShortAcceptance/check_export.py "$work/stereo" "$work/audio/programme-2.wav"
"$work/bin/short-acceptance" "$work/audio/programme-6.wav" "$work/surround" 6
python3 Tests/Profiling/ShortAcceptance/check_export.py "$work/surround" "$work/audio/programme-6.wav"
```

The source WAVs are deterministic. For repeated timings, keep the same inputs and use a new destination for every process. Run jobs sequentially on the same Mac. The driver writes stage timings, individual checks, waveform and proposal captures, the session, and the published verified audio/report bundle. The separate export check uses Python's standard library and requires assertions enabled.
