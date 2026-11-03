# M18 manual sections, role confirmation and preset dialogs

The main and Sections windows use one annotation document and audition service.
A waveform drag defines the interval on release. Either drag direction works.
Click to inspect, drag a violet boundary to move it, or Alt-drag the waveform
to pan. Ruler drag, wheel and scrollbar still pan. Each section is a contiguous
part of the programme. Draw within an existing section to retain surrounding
annotations. A range spanning an existing internal boundary reports the
conflict and preserves the map. Moving a boundary keeps both neighbours
nonempty. Zoom in or use the numeric frame editor for finer placement.

A fresh stereo or 5.1 import shows the roles in file-channel order. Supported
WAVE speaker masks, native FLAC order or an explicit FLAC mask, and supported
Mac AIFF CHAN layouts provide file evidence. When labels disagree with memory
the review shows both orders. The fallback is the last explicitly confirmed
order for that format on this Mac, then a standard unconfirmed suggestion.
Confirmed stereo and 5.1 orders have separate local preference files. Confirm
roles and analyse accepts the displayed assignment and starts Source analysis.
A separate-file import displays each filename and requires alignment confirmation.
The saved session owns the assignments on reopening. Routing is configured in
the output column independently.

### Manual checks

1. Save a named Target preset. The name window should close and the saved name
   should appear in the dropdown. Open Save as preset again and Cancel.
2. Drag across part of the waveform, then release. Try drawing in the reverse
   direction. Open Sections / gain and move the shared violet boundary. Zoom in
   and try another adjustment. Inspect Original and Current after these edits.
3. In Sections, click the upper timecode ruler while stopped and try Play, Stop
   and Space. The main window should show the same audition position. Click the
   waveform to inspect a passage independently.
4. Confirm a custom stereo order, then a different 5.1 order. Alternate fresh
   unlabelled imports and check the suggested roles before confirmation. Import
   a file with a conflicting speaker mask and inspect both orders. Reopen a
   saved project and verify its own roles. For separate files confirm alignment
   before completing the import.

Only the named standalone, timeline interaction and Target test targets are
built. The focused runner records the source and binary identities, component
renders, commands and reports. It uses an injected audio device and isolated
role and preset directories. Actual High Sierra, cross-Mac, physical audio/video
and full-duration/RSS qualification are deferred to the agreed combined block.
Device reconnection and the monitor-only stereo fold-down are planned in the
[playback roadmap](../../../../Docs/ProjectIndex/63_M18InteractiveTimelineAndProgressiveAnalysis.md#playback-roadmap-following-m18).

[Candidate 01](candidate-01/summary.json) passed **1,378 focused checks**, including
889 timeline, import and interaction checks and 489 Target/preset checks.
The standalone and both test binaries share build identity
`5811184d944d23b6b660dd6a4da55cfa24f81ffaa318414805d6c473fed92067`.
All contain arm64 and x86_64 code with Intel deployment minimum 10.13.
Verification ran on the development Mac. The source and test hashes were
checked again after the run.

The [role conflict](candidate-01/role-conflict.png), [separate-file review](candidate-01/separate-role-review.png)
and [drawn section](candidate-01/drawn-section.png) captures accompany the
interaction report. The early development run in build/file-labels.json exposed
a queued layout initialization notification that could clear a reviewed group.
The final selection covers that separate-file flow and alignment gate.

Development app
`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/m18-dev-01/LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app`
