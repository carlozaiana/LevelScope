#pragma once

// NOTE: Core header intentionally does NOT include juce_audio_processors.
// It only defines stable parameter IDs + numeric defaults/ranges.

namespace levelscope::mtdm
{
    // [BEGIN MTDM-PARAM-IDS]
    namespace ParamIDs
    {
        // Stable IDs (do not change once released)
        static constexpr const char* enabled     = "mtdm.enabled";
        static constexpr const char* thresholdDb = "mtdm.thresholdDb";
        static constexpr const char* ratio       = "mtdm.ratio";
    }

    namespace Defaults
    {
        static constexpr float enabled01     = 0.0f;   // 0=disabled, 1=enabled
        static constexpr float thresholdDb   = -24.0f; // placeholder for future DSP
        static constexpr float ratio         = 2.0f;   // placeholder for future DSP
    }

    namespace Ranges
    {
        static constexpr float thresholdMinDb = -80.0f;
        static constexpr float thresholdMaxDb =   0.0f;

        static constexpr float ratioMin = 1.0f;
        static constexpr float ratioMax = 20.0f;
    }
    // [END MTDM-PARAM-IDS]
} // namespace levelscope::mtdm