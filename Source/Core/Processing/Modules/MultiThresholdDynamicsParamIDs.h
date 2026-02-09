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

    // [BEGIN MTDM-PARAM-IDS-STAGE-D1A-ADD]
        namespace ParamIDs
        {
            // --- Stage D1a: T0–T1 upward zone + spectral upward compressor params (additive)
            static constexpr const char* t0Lufs = "mtdm.t0Lufs";
            static constexpr const char* t1Lufs = "mtdm.t1Lufs";

            // Spectral Upward Compressor (SUC)
            static constexpr const char* sucAmount01     = "mtdm.suc.amount01";
            static constexpr const char* sucMaxBoostDb   = "mtdm.suc.maxBoostDb";
            static constexpr const char* sucCurve        = "mtdm.suc.curve";
            static constexpr const char* sucLowKneeDb    = "mtdm.suc.lowKneeDb";
            static constexpr const char* sucHighKneeDb   = "mtdm.suc.highKneeDb";
            static constexpr const char* sucAttackMs     = "mtdm.suc.attackMs";
            static constexpr const char* sucReleaseMs    = "mtdm.suc.releaseMs";

            // “Advanced Quality” / mapping controls
            static constexpr const char* sucFftSizeChoice      = "mtdm.suc.fftSizeChoice";      // choice index
            static constexpr const char* sucBandsPerOctChoice  = "mtdm.suc.bandsPerOctChoice";  // choice index
            static constexpr const char* sucMinFreqHz          = "mtdm.suc.minFreqHz";
            static constexpr const char* sucMaxFreqHz          = "mtdm.suc.maxFreqHz";
        }

        namespace Defaults
        {
            // Zone thresholds shown in UI on LUFS axis (internally translated to spectral domain)
            static constexpr float t0Lufs = -45.0f;
            static constexpr float t1Lufs = -30.0f;

            // SUC defaults (audibly mild; still does something when enabled)
            static constexpr float sucAmount01   = 1.0f;
            static constexpr float sucMaxBoostDb = 8.0f;
            static constexpr float sucCurve      = 0.5f;
            static constexpr float sucLowKneeDb  = 3.0f;
            static constexpr float sucHighKneeDb = 3.0f;
            static constexpr float sucAttackMs   = 5.0f;
            static constexpr float sucReleaseMs  = 50.0f;

            // Advanced
            // fftSizeChoice: 0=1024, 1=2048, 2=4096, 3=8192 (see PluginProcessor layout)
            static constexpr int   sucFftSizeChoice     = 2; // 4096
            // bandsPerOctChoice: 0=1, 1=2, 2=3, 3=6
            static constexpr int   sucBandsPerOctChoice = 2; // 3 bands/oct (third-oct)
            static constexpr float sucMinFreqHz         = 25.0f;
            static constexpr float sucMaxFreqHz         = 20000.0f;
        }

        namespace Ranges
        {
            // Zone
            static constexpr float t0MinLufs = -80.0f;
            static constexpr float t0MaxLufs =   0.0f;

            static constexpr float t1MinLufs = -80.0f;
            static constexpr float t1MaxLufs =   0.0f;

            // SUC
            static constexpr float sucAmountMin01 = 0.0f;
            static constexpr float sucAmountMax01 = 1.0f;

            static constexpr float sucMaxBoostMinDb = 0.0f;
            static constexpr float sucMaxBoostMaxDb = 24.0f;

            static constexpr float sucCurveMin = 0.0f;
            static constexpr float sucCurveMax = 1.0f;

            static constexpr float sucKneeMinDb = 0.0f;
            static constexpr float sucKneeMaxDb = 24.0f;

            static constexpr float sucAttackMinMs  = 1.0f;
            static constexpr float sucAttackMaxMs  = 200.0f;

            static constexpr float sucReleaseMinMs = 5.0f;
            static constexpr float sucReleaseMaxMs = 1000.0f;

            static constexpr float sucMinFreqMinHz = 10.0f;
            static constexpr float sucMinFreqMaxHz = 2000.0f;

            static constexpr float sucMaxFreqMinHz = 1000.0f;
            static constexpr float sucMaxFreqMaxHz = 24000.0f;
        }

    // [BEGIN MTDM-PARAM-IDS-STAGE-D1A-TRIM-AND-CURVETYPE]
        namespace ParamIDs
        {
            // User calibration trim: shifts LUFS->spectral mapping.
            static constexpr const char* sucCalTrimDb        = "mtdm.suc.calTrimDb";

            // Curve type selection
            static constexpr const char* sucCurveTypeChoice  = "mtdm.suc.curveTypeChoice"; // 0=Monotonic, 1=Bell
        }

        namespace Defaults
        {
            static constexpr float sucCalTrimDb = 0.0f;
            static constexpr int   sucCurveTypeChoice = 0; // Monotonic
        }

        namespace Ranges
        {
            static constexpr float sucCalTrimMinDb = -12.0f;
            static constexpr float sucCalTrimMaxDb =  12.0f;
        }
    // [END MTDM-PARAM-IDS-STAGE-D1A-TRIM-AND-CURVETYPE]

    // [BEGIN MTDM-UPWARD-MODE-PARAM]
        namespace ParamIDs
        {
        // Upward processor mode:
        // 0 = Spectral (SUC), 1 = Broadband (time-domain upward)
            static constexpr const char* upwardModeChoice = "mtdm.upwardModeChoice";
        }

        namespace Defaults
        {
            static constexpr int upwardModeChoice = 0; // Spectral default (current behavior)
        }
    // [END MTDM-UPWARD-MODE-PARAM]

    // [BEGIN MTDM-LFE-MASK-PARAMS]
    namespace ParamIDs
    {
        // LFE routing policies (default false):
        // - lfeInDetector: include LFE in level measurement (detector)
        // - lfeInApply:    apply computed gain to LFE
        static constexpr const char* lfeInDetector = "mtdm.lfeInDetector";
        static constexpr const char* lfeInApply    = "mtdm.lfeInApply";
    }

    namespace Defaults
    {
        static constexpr float lfeInDetector01 = 0.0f;
        static constexpr float lfeInApply01    = 0.0f;
    }
    // [END MTDM-LFE-MASK-PARAMS]

    // [BEGIN MTDM-DOWNWARD-PARAMS]
    namespace ParamIDs
    {
        static constexpr const char* t2Lufs = "mtdm.t2Lufs";
        static constexpr const char* t3Lufs = "mtdm.t3Lufs";

        static constexpr const char* downEnabled01 = "mtdm.down.enabled";
        static constexpr const char* downRatio     = "mtdm.down.ratio";
        static constexpr const char* downKneeDb    = "mtdm.down.kneeDb";
        static constexpr const char* downAttackMs  = "mtdm.down.attackMs";
        static constexpr const char* downReleaseMs = "mtdm.down.releaseMs";
        static constexpr const char* downMakeupDb  = "mtdm.down.makeupDb";
    }

    namespace Defaults
    {
        static constexpr float t2Lufs = -12.0f;
        static constexpr float t3Lufs =  -6.0f;

        static constexpr float downEnabled01 = 0.0f; // OFF by default => no audible change
        static constexpr float downRatio     = 2.0f;
        static constexpr float downKneeDb    = 6.0f;
        static constexpr float downAttackMs  = 10.0f;
        static constexpr float downReleaseMs = 100.0f;
        static constexpr float downMakeupDb  = 0.0f;
    }

    namespace Ranges
    {
        static constexpr float t2MinLufs = -80.0f;
        static constexpr float t2MaxLufs =   0.0f;

        static constexpr float t3MinLufs = -80.0f;
        static constexpr float t3MaxLufs =   0.0f;

        static constexpr float downRatioMin = 1.0f;
        static constexpr float downRatioMax = 20.0f;

        static constexpr float downKneeMinDb = 0.0f;
        static constexpr float downKneeMaxDb = 24.0f;

        static constexpr float downAttackMinMs  = 1.0f;
        static constexpr float downAttackMaxMs  = 200.0f;

        static constexpr float downReleaseMinMs = 5.0f;
        static constexpr float downReleaseMaxMs = 1000.0f;

        static constexpr float downMakeupMinDb = -24.0f;
        static constexpr float downMakeupMaxDb =  24.0f;
    }
    // [END MTDM-DOWNWARD-PARAMS]

    // [END MTDM-PARAM-IDS-STAGE-D1A-ADD]
    // [END MTDM-PARAM-IDS]
} // namespace levelscope::mtdm