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

    // [BEGIN MTDM-UPWARD-ENABLED-BYPASS-PARAMS]
    namespace ParamIDs
    {
        // Structural enable (affects latency when Spectral).
        static constexpr const char* upEnabled01 = "mtdm.up.enabled01";

        // Safe bypass (automation-friendly): preserves latency when applicable.
        static constexpr const char* upBypass01  = "mtdm.up.bypass01";
    }

    namespace Defaults
    {
        static constexpr float upEnabled01 = 1.0f; // ON by default (matches existing behavior)
        static constexpr float upBypass01  = 0.0f;
    }
    // [END MTDM-UPWARD-ENABLED-BYPASS-PARAMS]

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

    // [BEGIN MTDM-ZONE-SOLO-MUTE-PARAMS]
    namespace ParamIDs
    {
        static constexpr const char* zoneSoloChoice      = "mtdm.zone.solo";
        static constexpr const char* zoneUpwardMute01    = "mtdm.zone.upward.mute";
        static constexpr const char* zoneDownwardMute01  = "mtdm.zone.downward.mute";
        static constexpr const char* zoneLimiterMute01   = "mtdm.zone.limiter.mute";
    }

    namespace Defaults
    {
        static constexpr int   zoneSoloChoice         = 0;   // 0=None, 1=Upward, 2=Downward, 3=Limiter
        static constexpr float zoneUpwardMute01       = 0.0f;
        static constexpr float zoneDownwardMute01     = 0.0f;
        static constexpr float zoneLimiterMute01      = 0.0f;
    }

    // [BEGIN MTDM-ZONE-UNTOUCHED-AUDITION-PARAM]
    namespace ParamIDs
    {
        static constexpr const char* zoneUntouchedMute01 = "mtdm.zone.untouched.mute";
    }

    namespace Defaults
    {
        static constexpr float zoneUntouchedMute01 = 0.0f;
    }
    // [END MTDM-ZONE-UNTOUCHED-AUDITION-PARAM]
    // [END MTDM-ZONE-SOLO-MUTE-PARAMS]

    // [BEGIN MTDM-ZONE-AUDITION-PARAMS]
    namespace ParamIDs
    {
        // Zone audition (combinable "solo zones"):
        // If any of these are ON => audition is active and ONLY the selected zones pass.
        static constexpr const char* zoneAudBelowT0_01 = "mtdm.zoneAud.belowT0";
        static constexpr const char* zoneAudT0T1_01    = "mtdm.zoneAud.t0t1";
        static constexpr const char* zoneAudT1T2_01    = "mtdm.zoneAud.t1t2";
        static constexpr const char* zoneAudT2T3_01    = "mtdm.zoneAud.t2t3";
        static constexpr const char* zoneAudAboveT3_01 = "mtdm.zoneAud.aboveT3";
    }

    namespace Defaults
    {
        static constexpr float zoneAudBelowT0_01 = 0.0f;
        static constexpr float zoneAudT0T1_01    = 0.0f;
        static constexpr float zoneAudT1T2_01    = 0.0f;
        static constexpr float zoneAudT2T3_01    = 0.0f;
        static constexpr float zoneAudAboveT3_01 = 0.0f;
    }
    // [END MTDM-ZONE-AUDITION-PARAMS]

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

    // [BEGIN MTDM-DOWNWARD-BYPASS-PARAM]
    namespace ParamIDs
    {
        static constexpr const char* downBypass01 = "mtdm.down.bypass01";
    }
    namespace Defaults
    {
        static constexpr float downBypass01 = 0.0f;
    }
    // [END MTDM-DOWNWARD-BYPASS-PARAM]

    // [BEGIN MTDM-LIMITER-PARAMS]
    namespace ParamIDs
    {
        static constexpr const char* limEnabled01    = "mtdm.lim.enabled";
        static constexpr const char* limCeilingDb    = "mtdm.lim.ceilingDb";
        static constexpr const char* limLookaheadMs  = "mtdm.lim.lookaheadMs";
        static constexpr const char* limReleaseMs    = "mtdm.lim.releaseMs";
    }

    namespace Defaults
    {
        static constexpr float limEnabled01   = 0.0f;   // OFF by default => no audible change
        static constexpr float limCeilingDb   = -1.0f;  // dBFS
        static constexpr float limLookaheadMs = 5.0f;   // ms
        static constexpr float limReleaseMs   = 100.0f; // ms
    }

    namespace Ranges
    {
        static constexpr float limCeilingMinDb = -20.0f;
        static constexpr float limCeilingMaxDb =   0.0f;

        static constexpr float limLookaheadMinMs = 0.0f;
        static constexpr float limLookaheadMaxMs = 50.0f;

        static constexpr float limReleaseMinMs = 5.0f;
        static constexpr float limReleaseMaxMs = 2000.0f;
    }

    // [BEGIN MTDM-LIMITER-PARAMS-TP-ATTACK-DRIVE]
    namespace ParamIDs
    {
        static constexpr const char* limAttackMs            = "mtdm.lim.attackMs";
        static constexpr const char* limDriveDb             = "mtdm.lim.driveDb";
        static constexpr const char* limOversamplingChoice  = "mtdm.lim.oversamplingChoice"; // 0=Off, 1=2x, 2=4x
    }

    namespace Defaults
    {
        static constexpr float limAttackMs = 0.5f;
        static constexpr float limDriveDb  = 0.0f;
        static constexpr int   limOversamplingChoice = 0; // Off
    }

    namespace Ranges
    {
        static constexpr float limAttackMinMs = 0.0f;
        static constexpr float limAttackMaxMs = 5.0f;

        static constexpr float limDriveMinDb = -24.0f;
        static constexpr float limDriveMaxDb =  24.0f;
    }
    // [END MTDM-LIMITER-PARAMS-TP-ATTACK-DRIVE]
    // [END MTDM-LIMITER-PARAMS]

    // [BEGIN MTDM-LIMITER-BYPASS-PARAM]
    namespace ParamIDs
    {
        static constexpr const char* limBypass01 = "mtdm.lim.bypass01";
    }
    namespace Defaults
    {
        static constexpr float limBypass01 = 0.0f;
    }
    // [END MTDM-LIMITER-BYPASS-PARAM]

    // [BEGIN MTDM-MC-POLICY-PARAMS]
    namespace ParamIDs
    {
        // Multichannel policy:
        // 0 = Linked (default)
        // 1 = Dialog-mask (detector from C or LCR; apply to C or LCR)
        // 2 = Unlinked (per-channel detectors)
        static constexpr const char* mcPolicyChoice = "mtdm.mcPolicyChoice";

        // Dialog-mask controls
        // 0 = C, 1 = LCR
        static constexpr const char* dialogDetectorChoice = "mtdm.dialogDetectorChoice";
        static constexpr const char* dialogApplyChoice    = "mtdm.dialogApplyChoice";
    }

    namespace Defaults
    {
        static constexpr int mcPolicyChoice        = 0; // Linked
        static constexpr int dialogDetectorChoice  = 0; // C
        static constexpr int dialogApplyChoice     = 0; // C
    }
    // [END MTDM-MC-POLICY-PARAMS]
    
    // [END MTDM-PARAM-IDS-STAGE-D1A-ADD]
    // [END MTDM-PARAM-IDS]
} // namespace levelscope::mtdm