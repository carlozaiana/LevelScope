#pragma once

// NOTE: Core header intentionally does NOT include juce_audio_processors.
// It only defines stable parameter IDs + numeric defaults/ranges.

namespace levelscope::lvlr
{
    // [BEGIN LVLR-PARAM-IDS]
    namespace ParamIDs
    {
        // Stable IDs (do not change once released)
        static constexpr const char* enabled           = "lvlr.enabled";
        static constexpr const char* targetLufs        = "lvlr.targetLufs";
        static constexpr const char* maxBoostDb        = "lvlr.maxBoostDb";
        static constexpr const char* maxCutDb          = "lvlr.maxCutDb";

        // Measurement selection:
        // 0 = Auto, 1 = Momentary, 2 = Short-term
        static constexpr const char* measChoice        = "lvlr.measChoice";

        // Operating mode:
        // 0 = Adaptive, 1 = Learn-Hold
        static constexpr const char* modeChoice        = "lvlr.modeChoice";

        // Learn toggle (used in Learn-Hold mode)
        static constexpr const char* learn01           = "lvlr.learn01";

        // Rate limiting for applied gain movement (dB/sec)
        static constexpr const char* rateUpDbPerSec    = "lvlr.rateUpDbPerSec";
        static constexpr const char* rateDownDbPerSec  = "lvlr.rateDownDbPerSec";
    }

    namespace Defaults
    {
        static constexpr float enabled01          = 0.0f;   // OFF by default => no audible change
        static constexpr float targetLufs         = -27.0f;
        static constexpr float maxBoostDb         = 12.0f;
        static constexpr float maxCutDb           = 12.0f;  // stored as positive magnitude

        static constexpr int   measChoice         = 0;      // Auto
        static constexpr int   modeChoice         = 0;      // Adaptive
        static constexpr float learn01            = 0.0f;   // false

        static constexpr float rateUpDbPerSec     = 1.0f;
        static constexpr float rateDownDbPerSec   = 3.0f;
    }

    namespace Ranges
    {
        static constexpr float targetMinLufs        = -48.0f;
        static constexpr float targetMaxLufs        = -12.0f;

        static constexpr float maxBoostMinDb        = 0.0f;
        static constexpr float maxBoostMaxDb        = 24.0f;

        static constexpr float maxCutMinDb          = 0.0f;
        static constexpr float maxCutMaxDb          = 24.0f;

        static constexpr float rateUpMinDbPerSec    = 0.1f;
        static constexpr float rateUpMaxDbPerSec    = 24.0f;

        static constexpr float rateDownMinDbPerSec  = 0.1f;
        static constexpr float rateDownMaxDbPerSec  = 24.0f;
    }
    // [END LVLR-PARAM-IDS]
} // namespace levelscope::lvlr