# LevelScope — Architecture Snapshot (Frozen Contracts + Implemented State)

**Repo/plugin name (current):** LevelScope  
**Intended suite/product name (later):** LevelFlow  
**Vendor (later):** Studio Paradiso  
**Rule:** keep “LevelScope” naming in code/CMake/UI until release.

This document is a **contract** for parallel development chats.  
Anything marked **Implemented** must be preserved.  
Anything marked **Planned** must not be assumed to exist until integrated.

Last updated: **2026-03-22**

---

## 1. Layers

### 1) LevelScopeCore (static lib; no GUI)

**Implemented (baseline, preserved):**
- Timeline truth storage (60 Hz) + persistence (chunked gzip state)
- Loudness/features computation (momentary/short-term/integrated/LRA/rolling LRA)
- Running program stats support (I running, LRA running, rolling LRA selector/numeric)

**Implemented (current):**
- Module chain scaffolding (ProcessorCore + module interfaces)
- Module 1 DSP: MTDM (multi-stage dynamics) in Core (reusable)
- Module 2 DSP: Leveler in Core (reusable)

**Planned / next milestones:**
- IN snapshot vs OUT live analysis tracks (separate histories/curves, separate visibility)
- Offline render/analysis hooks for standalone determinism
- Policy/proposal engine hooks
- Additional modules as needed (OutputTrim later, etc.)

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
- Host latency is reported and updated when latency-affecting params change (APVTS listener + message-thread timer)
- RT-safe input/output metering snapshots for UI
- Exposes UI-safe “effective playback” state (transport playing + callback freshness)

### 3) UI (JUCE Components)

**Implemented (baseline, preserved):**
- Timeline drawing + navigation (performance critical)
- Multi-resolution LOD history pyramid for fast drawing at all zoom levels

**Implemented (current):**
- Threshold overlay (T0–T3) over loudness history with ordered push-drag (min gap **0.1 LU**)
- Mission Control top strip:
  - target/current readouts
  - curve toggles + follow
  - rLRA window selector
  - multichannel policy controls + routing graphic (Detector vs Apply)
  - preset save/load (settings vs snapshot; see Persistence section)
- Right-side strip next to timeline:
  - LUFS scale
  - stage meters (at least In/Up/Dn/Lim/Out; Leveler metering available via snapshot API)
- Resizable rolling LRA lane (height divider)
- Bottom processing UI:
  - MTDM cards/panels
  - Leveler controls integrated
- Playback-lock policy in UI:
  - latency/structure controls disabled during **effective playback**
  - tooltips explain “Stop playback to change (changes latency)” etc.
  - preset load is disabled and guarded during effective playback
- Meter display behavior on stop/stale callbacks:
  - when callbacks go stale, UI meters decay toward rest (do not freeze forever)

**Implemented (current): UI view state persistence**
- UI layout + history viewport state is persisted via UIST chunk (see Persistence).

**Planned / next milestones:**
- IN vs OUT curve sets with toggles/labels
- Module strip UI (order/enable/bypass for multiple modules)
- 2D transfer/curve displays for processors (compressors/leveler)
- Advanced overlays (hotspots/gates/etc.)

---

## 2. Analysis Truth Model: IN snapshot vs OUT live

**Planned (not implemented yet):**
- IN (Input) snapshot captured from **unprocessed input** during a defined capture pass; frozen until re-captured
- OUT (Output) analysis computed from **post-processing output**; refreshable via capture/offline render

**Current behavior note:**
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
- Correct pattern: set a flag + use `AsyncUpdater` or a message-thread `Timer` to repaint.

UI “effective playback” gating (implemented):
- UI edit-lock decisions must use **effective playback**, not only last-known transport:
  - transport is playing AND audio callback is “fresh/recent”
- This avoids controls getting stuck disabled if a host stops calling `processBlock()` on stop.

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
- `MODG` — module graph (module IDs/order), schema v1
  - backward-compatible graph upgrade behavior: older sessions that only had MTDM can be upgraded to include Leveler before MTDM
- `UIST` — UI-only persisted state (layout + history viewport/toggles), schema v2
  - loader accepts schema v1 and v2; unknown schemas ignored

Container API:
- Core state container supports passing additive chunks via helper structs:
  - `ExtraStateChunksIn`
  - `ExtraStateChunksOut`

Rules:
- Never break old states: new chunks are optional with defaults
- Unknown chunks must be ignored on load (forward compatibility)

### Presets vs Snapshots (Implemented)
Two file semantics exist:

**Settings Preset**
- file extension: `.lscsettings`
- contains: `APVS + MODG`
- excludes: `UIST` and baseline history chunks
- loading settings presets must not change UI layout/view state

**Snapshot / DAW session**
- file extension: `.lscpreset` (snapshot)
- DAW session state is equivalent to full plugin state
- contains: baseline state + `APVS + MODG + UIST`

---

## 7. Module 1 (MTDM) — Implemented State

MTDM is a multistage dynamics/leveling module hosted as a Core module.

### Processing order (inside MTDM)
1) Upward stage (Spectral or Broadband; optional)
2) Downward stage (Broadband compressor; optional)
3) Limiter stage (lookahead; optional; output protection)
4) Post-chain **zone audition gate** (time-membership gating aligned with chain latency)

### Multichannel policies
- Linked: detector = all non-LFE; apply = all non-LFE (default)
- Dialog-mask: detector/apply selectable as C or LCR (fallback to non-LFE if not present)
- Unlinked: per-channel detector/gain states (advanced)

### Zone audition (implemented, important)
- Implemented via `mtdm.zoneAud.*` toggles (belowT0 / t0t1 / t1t2 / t2t3 / aboveT3)
- This is “solo/mute zones between thresholds” by timeline membership.
- Legacy “stage solo/mute” params are deprecated/unneeded; keep IDs for backward session compatibility but do not expose.

### Latency reporting & updates (implemented)
- Host latency computed from latency-affecting params and updated via APVTS listener + message-thread timer.
- MTDM internal delayed features (zone audition) are aligned to the **true effective chain latency**, including limiter FIR detector delay when oversampling is enabled.

### Stop-playback edit policy (implemented)
Certain params are treated as “structural/latency/quality” and are:
- **Non-automatable** (backend)
- **Playback-locked** in UI using effective playback state

Structural enable params marked non-automatable:
- `mtdm.enabled`
- `mtdm.up.enabled01`
- `mtdm.lim.enabled`

Structural/quality params marked non-automatable:
- `mtdm.upwardModeChoice`
- `mtdm.suc.fftSizeChoice`
- `mtdm.suc.bandsPerOctChoice`
- `mtdm.suc.minFreqHz`
- `mtdm.suc.maxFreqHz`
- `mtdm.lim.lookaheadMs`
- `mtdm.lim.oversamplingChoice`

### Metering (implemented)
- RT-safe snapshots for:
  - I/O peak+RMS (pre and post processing)
  - Upward boost
  - Downward gain reduction (excluding makeup)
  - Limiter gain reduction

---

## 8. Module 2 (Leveler) — Implemented State

Leveler is a zero-latency Core module inserted **before MTDM** by default.

### Graph position (default)
1) Leveler  
2) MTDM

### Leveler parameters (lvlr.*)

**Core**
- `lvlr.enabled` (default OFF)
- `lvlr.targetLufs` (default -27.0)
- `lvlr.maxBoostDb` (default 12.0)
- `lvlr.maxCutDb` (default 12.0 magnitude)

**Measurement / algorithm**
- `lvlr.measChoice` = Auto / Momentary / Short-term (non-automatable)
- `lvlr.modeChoice` = Adaptive / Learn-Hold (non-automatable)
- `lvlr.learn01` (non-automatable)
- `lvlr.rateUpDbPerSec` (default 1.0)
- `lvlr.rateDownDbPerSec` (default 3.0)

**Control source + host lane capture**
- `lvlr.controlModeChoice` = Internal / Host Gain (non-automatable)
- `lvlr.hostGainDb` (automatable, default 0.0 dB, range approx -24..+24 dB)
- `lvlr.captureToHost01` (non-automatable; explicit arm/toggle; default OFF)

Important: no `lvlr.appliedGainDb` APVTS parameter exists. “Applied gain” is exposed via metering snapshots only (telemetry), not as a host-controlled parameter.

### Detector / measurement (implemented)
- Self-contained detector (does not depend on timeline/history model)
- BS.1770-oriented K-weighted path; LFE excluded by default unless routing override includes it
- Internal detector update cadence: 60 Hz
- Momentary window: 0.4 s (24 frames)
- Short-term window: 3.0 s (180 frames)
- Auto measurement rule: use Short-term when valid; otherwise use Momentary

### Modes (implemented)
- Adaptive: continuously measures and applies bounded/rate-limited correction
- Learn-Hold:
  - while learn01 true: measure and compute candidate
  - on learn01 false: commit one held gain and apply until next learn pass
  - commit also occurs on playing→stopped edge when learning (best-effort)

### Control modes (implemented)
- **Internal**
  - Leveler computes and applies internal gain (Adaptive/Learn-Hold behavior)
  - `lvlr.hostGainDb` is not used as the audio gain source
- **Host Gain**
  - Leveler ignores internal gain decisions for the applied gain
  - applied gain comes from `lvlr.hostGainDb` and is smoothed with a per-block ramp to avoid zippering

### Capture: Adaptive → Host Gain lane (implemented)
- When:
  - `lvlr.controlModeChoice == Internal`
  - `lvlr.captureToHost01 == true`
  - transport is **effectively playing**
- The processor periodically writes the Leveler’s **actual applied gain** into `lvlr.hostGainDb` using proper host gesture semantics:
  - beginChangeGesture → periodic setValueNotifyingHost → endChangeGesture
- Capture cadence: 30 Hz, with a write threshold (≈ 0.1 dB) to limit automation density.
- Host recording of plugin-driven parameter changes is DAW-dependent.

**Session safety rule (implemented):**
- `lvlr.captureToHost01` is forcibly disarmed on state load/preset restore and other safety transitions (e.g., leaving Internal mode / stopping capture) to prevent reopening a session and unintentionally overwriting the host gain automation lane.

### Gain behavior (implemented)
- Rate-limited in **dB/sec** (separate up/down) for internal gain changes.
- Host Gain mode uses a per-block ramp to avoid zipper noise.

### Routing (implemented; reused shared controls)
Leveler v1 reuses existing shared routing policy params:
- `mtdm.mcPolicyChoice`, dialog choices, and LFE overrides
- Supports Linked / Dialog-mask / Unlinked policies consistent with MTDM

### Metering (implemented)
- Leveler gain snapshot:
  - `gainDbCurrent`, `gainDbBlockPeak`, `gainDbHold`
  - sign: positive = boost, negative = cut
- Meter reflects the **actual applied gain** in both Internal and Host Gain modes.

---

## 9. Coding & Patch Workflow (Hard)

- Large files must not be replaced wholesale.
- Add edit anchors:
  - `// [BEGIN <BLOCK-ID>]` ... `// [END <BLOCK-ID>]`
- Prefer patch-style edits: replace only marked blocks or specific functions.

---

## 10. Roadmap Chats (Parallel Work Split)

### Chat A: Modules / DSP
- Refine Leveler + MTDM (quality/performance)
- Add 2D curve displays (UI support)
- Future modules if needed

### Chat B: UI Overhaul Milestone (later)
- IN vs OUT curves + toggles
- module strip (order/enable/bypass)
- overlays/hotspots/policy visualizations

### Chat C: Logic/Policy Proposal Engine
- Presets/specs + proposal outputs (static params + automation curves)
- Operates on analysis snapshots (IN and/or OUT)
- Deterministic standalone/offline behavior