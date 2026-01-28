#pragma once

// LevelScopeCore umbrella header for the “context pack”.
// NOTE: Keep this umbrella stable and intentionally curated.
// - Do NOT include GUI headers here.
// - Do NOT include huge implementation headers here.
// - Prefer each included header to be self-contained (include what it uses).

// Existing core headers (baseline)
#include "LevelScopeHistoryModel.h"
#include "BS1770KWeighting.h"
#include "RunningLoudnessStats.h"

// New processing/module-chain API
#include "Processing/ProcessContext.h"
#include "Processing/IAudioModule.h"
#include "Processing/ProcessorCore.h"