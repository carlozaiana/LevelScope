# M18 timeline section gains and playback keyboard focus

The main timeline now displays Current's all-channel section gain **before global
trim**, irrespective of the inspected measurement or audition domain. The solid
line is the applied audio envelope; the dashed line is a fresh numerical draft.
Both use M12's canonical transition frames, including shortened ramps. Annotation
boundaries do not move existing audio ramps. At crowded zoom levels the graph keeps
the first, last, smallest and largest gain knots per visible column in time order
(maximum 4,096 columns). Clicking a section is inspection only.

The UI stages finite values in the existing −12 to +12 dB range; **Apply gains**
and **Discard gains** act on the whole section draft. Changing Apply invalidates
Current only. Identical Apply and Discard preserve the accepted measurement.
Unsupported layouts stay inspectable and report why manual editing is unavailable.
The draft edit serial is runtime-only; there is no DSP or session schema change.

Playback preparation still makes and verifies its own protected source copy.
The status now says “Verifying audio for playback...”. It compares the copied
file content and confirms that the audio decodes as imported; it does not perform
another loudness measurement. The audition service reuses a matching prepared
copy. Seek and playback toggles now put focus on a stable timeline component
before the Seek button disables, leaving text and modal input alone.

### Manual checklist

1. Import a supported programme, confirm its roles and let Original and Current
   become Ready. Open **Sections / gain…** if you want to split the default
   whole-programme annotation. Click the new gain lane or use its arrow buttons.
   Check the section name, frame range and applied gain. **Fit selection** fits
   the section on the same waveform/seconds axis without moving audition.
2. Enter a decimal such as `-1.23456789` and press Return or **Set draft**.
   Select a second section to stage another value. Only the dashed envelope
   changes. Try **Discard gains**: the solid curve and Current measurement stay
   as they were.
3. Stage a change and choose **Apply gains** while stopped. Current becomes
   Unverified, Original remains Ready, and Current audition uses the solid
   envelope. **Re-measure Current State** verifies it independently. Zoom into
   a transition and inspect the applied dB at a frame versus the section setting.
4. With an output configured, enter seconds and click **Seek**. Press Space to
   play, then Space again to stop without clicking elsewhere. On first play,
   look for the verification status. With a numeric field focused, Space should
   remain text input.

Run `python3 -B Tests/Reference/results/m18-gain-lane/run_slice.py <new-directory>`
after building only LevelScopeTimelineProgressTests, LevelScopeStandalone and
LevelScopeSectionProcessingTests in the existing universal Release build. It records
source/test/binary hashes, architecture/Intel minimum, exact commands, component
captures and JSON reports. Fixtures use a new private /tmp directory; injected
audio and a JUCE test window exercise keyboard focus. No physical output or video
is opened. The standalone app path is
`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/m18-dev-01/LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app`.

The early iteration in `pilot-01` exposed a slow full-path dashed renderer on
10,000 sections. Its timing is not final. The renderer now bounds the displayed
knots and draws dashes directly in screen space. Only the final run qualifies
this development slice. Previous waveform, loudness arithmetic, session/export
and proposal checks are reused for unchanged paths. Broader native interaction,
actual High Sierra/cross-Mac, physical device/video and full-duration/RSS
qualification remain deferred. Saved applied-proposal explanations are still
separate M18 work.

### Final development result

[Candidate 01](candidate-01/summary.json) passed **715 timeline checks and
334 section-processing checks** (1,049 total). Its standalone app and both test
executables contain arm64/x86_64 and declare Intel minimum 10.13; tests ran on
the development Mac. Measurement build identity:
`42315ef049da2fc4465d5840652fac210f48b7bd70e4c5024ba21196e812a4f6`.
The same inputs were verified after the run. The single dense 10,000-section
preparation and software paint sample was **6.041 ms**; it is a bounded fixture,
not a full-duration memory or display qualification.

The [timeline report](candidate-01/timeline.json) includes stopped numeric seek,
Space/start/Space/stop through a real JUCE component peer with injected audio,
text input, section selection, numerical precision across both editors, a stale
draft after a boundary edit, invalid input, no-op/discard/change, cancelled and
independent verification, unsupported layout and dense drawing. The
[processing smoke report](candidate-01/section-processing.json) checks the M12
audio contracts. [Draft](candidate-01/gain-draft.png) and
[applied](candidate-01/gain-applied.png) captures show the gain lane aligned with
actual Original waveform and loudness data. The later native/device/High Sierra
qualification block remains unchanged.
