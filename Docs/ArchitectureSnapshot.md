# LevelScope — Architecture Snapshot (Frozen Contracts + Implemented State)

**Repo/plugin name (current):** LevelScope  
**Intended suite/product name (later):** LevelFlow  
**Vendor (later):** Studio Paradiso  
**Rule:** keep “LevelScope” naming in code/CMake/UI until release.

This document is a **contract** for parallel development chats.  
Anything marked **Implemented** must be preserved.  
Anything marked **Planned** must not be assumed to exist until integrated.

Last updated: **2026-03-11**

---

## 1. Layers

### 1) LevelScopeCore (static lib; no GUI)

**Implemented (baseline, preserved):**
- Timeline truth storage (60 Hz) + persistence (chunked gzip state)
- Loudness/features computation (momentary/short-term/integrated/LRA/rolling LRA)
- Running program stats support (I running, LRA running, rolling LRA selector/numeric)

**Implemented (current):**
- Module DSP chain scaffolding (ProcessorCore + module interfaces)
- Module 1 DSP lives in Core and is reusable by plugin/standalone targets (no GUI dependencies)

**Planned / next milestones:**
- IN snapshot vs OUT live analysis tracks (separate histories/curves, separate visibility)
- Offline render/analysis hooks for standalone determinism
- Policy/proposal engine hooks

Hard rule: Core must be usable by plugin + future standalone + headless tests.

### 2) Plugin wrapper (JUCE AudioProcessor)

**Implemented (baseline, preserved):**
- DAW timebase/transport/discontinuity logic
- transport start ramp guard
- stop-time silence freeze (do not corrupt analysis)
- discontinuity detection (seek/loop)
- smooth zoom/pan/follow playhead

**Implemented (current):**
- Hosts ProcessorCore and runs module graph inside `processBlock()`
- Owns APVTS and binds raw atomics into Core modules RT-safely
- Host latency is reported and updated when latency-affecting params change
- RT-safe input/output metering snapshots for UI

**Planned / next milestones:**
- Capture-pass controls for IN snapshot analysis
- OUT live analysis track refresh workflows

### 3) UI (JUCE Components)

**Implemented (baseline, preserved):**
- Timeline drawing + navigation (performance critical)
- Multi-resolution LOD history pyramid for fast drawing at all zoom levels

**Implemented (current):**
- Threshold overlay (T0–T3) over loudness history with ordered push-drag (min gap 0.1 LU)
- Mission Control top strip (targets/current readouts, curve toggles, follow, rLRA window, multichannel policy controls + routing graphic)
- Right-side strip next to timeline: LUFS scale + stage meters (In/Up/Dn/Lim/Out)
- Resizable rolling LRA lane (height divider)
- Bottom MTDM panel (card-based controls with internal resizers)

**Planned / next milestones:**
- IN vs OUT curve sets with toggles/labels
- Module strip UI (order/enable/bypass for multiple modules)
- Advanced overlays (hotspots/gates/etc.)
- UI-only state persistence chunk (splitter positions, lane heights, strip widths)

---

## 2. Analysis Truth Model: IN snapshot vs OUT live

**Planned (not implemented yet):**
- IN (Input) snapshot captured from **unprocessed input** during a defined capture pass; frozen until re-captured
- OUT (Output) analysis computed from **post-processing output**; refreshable via capture/offline render

Current behavior note:
- Loudness/history analysis in the plugin is currently computed on the **input-side** (pre-module-chain).
- I/O meters reflect pre- and post-module-chain audio.

---

## 3. Threading / RT Safety Rules (Hard)

Audio thread must:
- No locks (no mutex/spinlock), no blocking waits
- No heap allocations (no new/delete, no growing std::vector, no juce::Array growth)
- No file I/O, no logging
- No ValueTree access or mutation
- Only read atomics / prebuilt structures

Module graph swaps:
- Use an **RCU-style** graph swap (read-copy-update).
- In C++17, implemented via `std::atomic_load/std::atomic_store` free functions on `std::shared_ptr<T>` (C++11+), or equivalent.
- Audio thread is read-only + lock-free; non-audio thread builds and swaps.

APVTS listeners:
- `AudioProcessorValueTreeState::Listener::parameterChanged()` may be called from non-message threads.
- UI must not call `repaint()` or touch Components directly from `parameterChanged()`.
- Correct pattern: set an atomic/flag + use `AsyncUpdater` or a message-thread `Timer` to repaint.

---

## 4. Multichannel Scope (Hard)

**Implemented (current):**
- Plugin supports non-ambisonic layouts up to **12 channels** (7.1.4) where main input layout == main output layout.
- Channel role mapping uses `juce::AudioChannelSet::getTypeOfChannel (channelIndex)`; no channel-order assumptions.

Loudness analysis:
- Follows BS.1770 principles; LFE excluded from loudness by default.

Processing:
- Channel-role-aware via detector/apply masks (bitmask/list based).
- Default LFE policy:
  - excluded from detector and apply in dynamics stages unless user opts in
  - limiter is output protection and applies to **all channels including LFE** by default (current design choice)

Development staging rule: **5.1-first** (stereo/mono remain quick test paths).

---

## 5. Module Chain Contracts (Frozen)

### ProcessContext (per block) — Implemented
Contains:
- `juce::AudioBuffer<float>& audio`, optional `juce::MidiBuffer*`
- `sampleRate`, `numSamples`
- `juce::AudioChannelSet channelSet`
- transport flags + discontinuity flag
- absolute sample index and optional 60 Hz frame index

### IAudioModule — Implemented
Required methods:
- stable `getModuleID()` (persistence ID; must not change once released)
- `getDisplayName()`
- `prepare(spec)`, `reset()`
- `process(ctx) noexcept` (**RT-safe; no allocations, no locks**)
- bypass set/get (RT-safe)
- `saveState()/loadState()` (non-audio-thread only)

### ProcessorCore — Implemented + integrated
- Owns active module graph snapshot
- Audio thread reads active graph snapshot and iterates modules
- Non-audio thread builds a new graph and swaps it in (RCU-style)

---

## 6. Persistence Strategy (Additive, Versioned)

**Implemented (baseline, preserved):**
- Chunked gzip plugin state container with stable IDs.
- Existing chunk IDs (baseline; do not break):
  - `LSCP` container magic/version
  - `HIST` (energy/history)
  - `LRAG` (LRA gate curve)
  - `TCOF` (user time offset)

**Implemented (current additive chunks):**
- `APVS` — APVTS state (ValueTree binary), schema v1
- `MODG` — module graph (module IDs/order; bypass byte stored but currently not relied upon for restore), schema v1

Rules:
- Never break old states: new chunks are optional with defaults
- Unknown chunks must be ignored on load (forward compatibility)

Presets (current UI behavior):
- `.lscpreset` save/load uses the **full plugin state blob** (equivalent to getStateInformation/setStateInformation).
- Planned refinement: “Settings preset” (APVS+MODG only) vs “Session snapshot” (full blob including history).

Planned:
- `UIST` (or similar) UI-only preferences chunk (splitter positions, lane heights, strip widths, card heights)

---

## 7. Module 1 (MTDM) — Implemented State

Module 1 is a multistage dynamics/leveling prototype hosted as a Core module.

### Processing order (inside MTDM)
1) Upward stage (Spectral or Broadband; optional)
2) Downward stage (Broadband compressor; optional)
3) Limiter stage (lookahead; optional; output protection)
4) Post-chain **zone audition gate** (time-membership gating aligned with chain latency)

### Upward stage
- Spectral (STFT-based) or Broadband (time-domain)
- Supports Linked / Dialog-mask / Unlinked policies via detector/apply masks
- Spectral mode contributes FFT latency; provides delay-preserving audition bypass (unity through pipeline)

### Downward stage
- Broadband downward compressor with T2–T3 engagement zone
- Supports Linked / Dialog-mask / Unlinked via masks
- T2/T3 are clamped/reordered at runtime to preserve T1–T2 untouched semantics

### Limiter stage
- Lookahead limiter with drive, attack ramping (back-scheduled), release smoothing
- Optional FIR oversampling on detector path (true-peak-ish)
- Applies to all channels by default (output protection)
- Delay-preserving audition bypass (unity through delay pipeline when bypassed but enabled)

### Multichannel policies
- Linked: detector = all non-LFE; apply = all non-LFE (default)
- Dialog-mask: detector/apply selectable as C or LCR (fallback to non-LFE if not present)
- Unlinked: per-channel detector/gain states (advanced)

### Zone audition (important)
- Implemented via `mtdm.zoneAud.*` toggles (belowT0 / t0t1 / t1t2 / t2t3 / aboveT3)
- This is “solo/mute zones between thresholds” by timeline membership.
- Legacy “stage solo/mute” params exist but are **deprecated** (not used by DSP/UI); plan is to keep IDs for backward session compatibility but not expose them.

### Latency reporting & updates
- Host latency computed from latency-affecting params and updated via APVTS listener + 10 Hz message-thread timer.
- Latency-affecting params include MTDM enabled, upward enabled/mode/FFT size, limiter enabled/lookahead/oversampling.

Known watch item:
- Ensure any internal delay alignment (e.g. zone audition gate) accounts for all latency components (including oversampling FIR detector delay), consistent with host-reported latency.

### Metering
- RT-safe snapshots for:
  - I/O peak+RMS (pre and post processing)
  - Upward boost
  - Downward gain reduction (excluding makeup)
  - Limiter gain reduction
- UI polls snapshots on message thread; audio thread writes atomics once per block.

---

## 8. Coding & Patch Workflow (Hard)

- Large files must not be replaced wholesale.
- Add edit anchors:
  - `// [BEGIN <BLOCK-ID>]` ... `// [END <BLOCK-ID>]`
- Prefer patch-style edits: replace only marked blocks or specific functions.

---

## 9. Roadmap Chats (Parallel Work Split)

### Chat A: Modules / DSP
- Extend MTDM correctness (latency alignment edge cases, presets semantics, RT audits)
- Add future modules (Module 2+)

### Chat B: UI Overhaul Milestone (later)
- IN vs OUT curves + toggles
- module strip (order/enable/bypass)
- overlays/hotspots/policy visualizations
- UI state persistence chunk

### Chat C: Logic/Policy Proposal Engine
- Presets/specs + proposal outputs (static params + automation curves)
- Operates on analysis snapshots (IN and/or OUT)
- Deterministic standalone/offline behavior