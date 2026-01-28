# LevelScope — Architecture Snapshot (Frozen Contracts)

**Repo/plugin name (current):** LevelScope  
**Intended suite/product name (later):** LevelFlow  
**Vendor (later):** Studio Paradiso  
**Rule:** keep “LevelScope” naming in code/CMake/UI until release.

This document is a **contract** for parallel development chats.  
Anything marked **Implemented (baseline)** must be preserved.  
Anything marked **Planned** must not be assumed to exist until integrated.

---

## 1. Layers

### 1) LevelScopeCore (static lib; no GUI)

**Implemented (baseline):**
- Timeline truth storage (60 Hz) + persistence (chunked gzip state)
- Loudness/features computation (momentary/short-term/integrated/LRA/rolling LRA)
- Running program stats panel support (I running, LRA running, rolling LRA selector/numeric)

**Planned / next milestones:**
- IN snapshot vs OUT live analysis tracks (separate timelines/curves, separate visibility)
- Module DSP chain host (ProcessorCore) + modules
- Future: offline analysis/render hooks, policy/proposal engine hooks

Hard rule: Core must be usable by plugin + future standalone + headless tests.

### 2) Plugin wrapper (JUCE AudioProcessor)

**Implemented (baseline):**
- DAW timebase/transport/discontinuity logic
- transport start ramp guard
- stop-time silence freeze (do not corrupt analysis)
- discontinuity detection (seek/loop)
- smooth zoom/pan/follow playhead

**Planned / next milestones:**
- Controls capture passes for IN snapshot analysis
- Feeds audio blocks through DSP module chain
- Owns APVTS and maps parameters to modules (later)

### 3) UI (JUCE Components)

**Implemented (baseline):**
- Timeline drawing + navigation (performance critical)
- Multi-resolution LOD history pyramid for fast drawing at all zoom levels

**Planned / next milestones:**
- Timeline supports IN vs OUT curve sets with toggles/labels
- Module strip UI (order/enable/bypass)
- Advanced overlays (gates/hotspots/etc.)

---

## 2. Analysis Truth Model: IN snapshot vs OUT live

### IN (Input) snapshot (frozen) — **Planned**
- Captured from **unprocessed input** during a defined capture pass
- Frozen until explicitly re-captured
- Used as stable “timeline truth reference”

### OUT (Output) live/refreshable — **Planned**
- Computed from **post-processing output** (post module chain)
- Can become stale after upstream edits/parameter changes
- Refresh requires capture or offline/standalone render

Practical ordering note (target behavior):
- IN snapshot capture is done from input stream **before** the module chain.
- OUT analysis is computed from stream **after** the module chain.

UI requirement (later): show both curve sets with independent toggles/labels.

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
- In C++17 this is typically implemented via `std::atomic_load/std::atomic_store` free functions on `std::shared_ptr<T>` (standard since C++11), or an equivalent safe atomic pointer mechanism.
- Exact mechanism must preserve: audio thread is read-only + lock-free; non-audio thread builds and swaps.

---

## 4. Multichannel Scope (Hard)

- Support arbitrary `juce::AudioChannelSet` layouts up to **7.1.4** (incl. 7.1.2).
- Channel role mapping uses `juce::AudioChannelSet::getTypeOfChannel (channelIndex)`.
- Loudness analysis follows BS.1770 principles:
  - LFE excluded from loudness by default (configurable later)
- Processing must be channel-role-aware:
  - LFE default excluded from **linked gain** application unless user opts in later.

Development staging rule: **5.1-first** (stereo/mono remain quick test paths).

---

## 5. Module Chain Contracts (Frozen)

### ProcessContext (per block) — **Implemented (scaffolding)**
Contains:
- `juce::AudioBuffer<float>& audio`, optional `juce::MidiBuffer*`
- `sampleRate`, `numSamples`
- `juce::AudioChannelSet channelSet`
- transport flags + discontinuity flag
- absolute sample index and optional 60 Hz frame index

### IAudioModule — **Implemented (scaffolding)**
Required methods:
- `getModuleID()` (stable persistence ID; must not change once released)
- `getDisplayName()`
- `prepare(spec)`, `reset()`
- `process(ctx) noexcept` (**RT-safe; no allocations, no locks**)
- bypass set/get (RT-safe)
- `saveState()/loadState()` (non-audio-thread only)

### ProcessorCore — **Implemented (scaffolding), integration planned**
- Owns active module graph snapshot
- Audio thread reads active graph snapshot and iterates modules
- Non-audio thread builds a new graph and swaps it in (RCU-style)

---

## 6. Persistence Strategy (Additive, Versioned)

**Implemented (baseline):**
- Chunked gzip plugin state container with stable IDs.
- Existing chunk IDs (baseline; do not break):
  - `LSCP` container magic/version
  - `HIST` (energy/history)
  - `LRAG` (LRA gate curve)
  - `TCOF` (user time offset)

**Planned / next milestones:**
- Extend state with additive sections:
  - Timeline truth (IN)
  - Timeline truth (OUT)
  - Module graph (order, enabled/bypass, per-module state)
  - UI preferences (non-audio critical)

Rules:
- Never break old states: new fields are optional with defaults
- Any new binary/state format increments a schema version
- Module states are keyed by `moduleID` and may have per-module state versioning

---

## 7. Coding & Patch Workflow (Hard)

- Large files must not be replaced wholesale.
- Add edit anchors:
  - `// [BEGIN <BLOCK-ID>]` ... `// [END <BLOCK-ID>]`
- Prefer patch-style edits: replace only marked blocks or specific functions.
- Keep audio-thread code visually obvious (comments + `noexcept`).

---

## 8. Roadmap Chats (Parallel Work Split)

### Chat A: Module 1 (Multi-threshold dynamics)
- Implement module DSP + detector/linking policies (linked/dialog-mask/unlinked)
- Add APVTS parameter layout + binding
- Integrate ProcessorCore into PluginProcessor with minimal disruption
- Preserve existing analysis/timeline stability

### Chat B: UI Overhaul Milestone
- Display IN snapshot vs OUT live curves with toggles
- Module strip UI (order/enable/bypass)
- Advanced overlays (gate curves/hotspots later)
- Preserve current timeline LOD performance and navigation guarantees

### Chat C: Logic/Policy Proposal Engine
- Define presets/specs + proposal outputs (static params + automation curves)
- Operates on analysis snapshots (IN and/or OUT)
- Produces recommended settings without directly embedding UI/DSP details
