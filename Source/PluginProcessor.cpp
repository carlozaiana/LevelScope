#include "PluginProcessor.h"
#include "PluginEditor.h"

// [BEGIN LVLR-PLUGIN-INCLUDE]
#include "Core/Processing/Modules/LevelerModule.h"
// [END LVLR-PLUGIN-INCLUDE]

// [BEGIN MTDM-PLUGIN-INCLUDE]
#include "Core/Processing/Modules/MultiThresholdDynamicsModule.h"
// [END MTDM-PLUGIN-INCLUDE]

// [BEGIN LVLR-PARAMIDS-INCLUDE]
#include "Core/Processing/Modules/LevelerParamIDs.h"
// [END LVLR-PARAMIDS-INCLUDE]

// [BEGIN MTDM-PARAMIDS-INCLUDE]
#include "Core/Processing/Modules/MultiThresholdDynamicsParamIDs.h"
// [END MTDM-PARAMIDS-INCLUDE]

#include <cmath>
#include <algorithm>
#include <type_traits>
#include <cstring>

namespace
{
    // Treat arithmetic types (int, int64, bool, double) as always "present".
    template <typename T, std::enable_if_t<std::is_arithmetic_v<std::decay_t<T>>, int> = 0>
    bool optHasValue (T) noexcept { return true; }

    template <typename T, std::enable_if_t<std::is_arithmetic_v<std::decay_t<T>>, int> = 0>
    T optValue (T v) noexcept { return v; }

    // Treat optional-like types (std::optional, juce::Optional, etc.) as "present" if bool(o) is true.
    template <typename Opt, std::enable_if_t<!std::is_arithmetic_v<std::decay_t<Opt>>, int> = 0>
    bool optHasValue (const Opt& o) noexcept { return (bool) o; }

    template <typename Opt, std::enable_if_t<!std::is_arithmetic_v<std::decay_t<Opt>>, int> = 0>
    auto optValue (const Opt& o) noexcept -> decltype(*o) { return *o; }

    // [BEGIN LS-IO-METERING-HELPER]
    struct BlockPeakRms
    {
        float peakDb = -200.0f;
        float rmsDb  = -200.0f;
    };

    static BlockPeakRms computeBlockPeakRmsMaxAcrossChannels (const juce::AudioBuffer<float>& b,
                                                              int numCh,
                                                              int numSamples) noexcept
    {
        BlockPeakRms out;
        if (numCh <= 0 || numSamples <= 0)
            return out;

        numCh = juce::jlimit (1, 12, numCh);

        float peakMax = 0.0f;
        float rmsMax  = 0.0f;

        for (int ch = 0; ch < numCh; ++ch)
        {
            const float* d = b.getReadPointer (ch);
            float p = 0.0f;
            double ss = 0.0;

            for (int i = 0; i < numSamples; ++i)
            {
                const float x = d[i];
                const float ax = std::abs (x);
                if (ax > p) p = ax;
                ss += (double) x * (double) x;
            }

            const float rms = (float) std::sqrt (ss / (double) numSamples);

            if (p > peakMax) peakMax = p;
            if (rms > rmsMax) rmsMax = rms;
        }

        out.peakDb = juce::Decibels::gainToDecibels (peakMax, -200.0f);
        out.rmsDb  = juce::Decibels::gainToDecibels (rmsMax,  -200.0f);
        return out;
    }

    static float decayHoldDb (float holdDb, float currentDb, float decayDb) noexcept
    {
        // Linear decay in dB for simplicity (stable and cheap).
        holdDb = juce::jmax (holdDb - decayDb, -200.0f);
        return juce::jmax (holdDb, currentDb);
    }
    // [END LS-IO-METERING-HELPER]
}

// [BEGIN MTDM-APVTS-PARAM-LAYOUT]
juce::AudioProcessorValueTreeState::ParameterLayout LevelScopeAudioProcessor::createParameterLayout()
{
    juce::AudioProcessorValueTreeState::ParameterLayout layout;

    using namespace levelscope::mtdm;

    // [BEGIN LVLR-STRUCTURAL-PARAM-SOURCE-OF-TRUTH]
    // Leveler parameter policy/source-of-truth block.
    // Keep this block in sync with:
    //   1) APVTS automatable/non-automatable flags
    //   2) UI enable/disable rules
    //   3) module binding order
    //   4) host-gain capture logic in PluginProcessor timerCallback()
    //
    // Leveler-specific params
    // - levelscope::lvlr::ParamIDs::enabled
    // - levelscope::lvlr::ParamIDs::targetLufs
    // - levelscope::lvlr::ParamIDs::maxBoostDb
    // - levelscope::lvlr::ParamIDs::maxCutDb
    // - levelscope::lvlr::ParamIDs::measChoice          // non-automatable
    // - levelscope::lvlr::ParamIDs::modeChoice          // non-automatable
    // - levelscope::lvlr::ParamIDs::learn01             // non-automatable
    // - levelscope::lvlr::ParamIDs::rateUpDbPerSec
    // - levelscope::lvlr::ParamIDs::rateDownDbPerSec
    // - levelscope::lvlr::ParamIDs::controlModeChoice   // non-automatable
    // - levelscope::lvlr::ParamIDs::hostGainDb          // automatable host-editable lane
    // - levelscope::lvlr::ParamIDs::captureToHost01     // non-automatable; forcibly disarmed on state load
    //
    // Reused shared routing-policy params (currently under mtdm.* namespace)
    // - ParamIDs::mcPolicyChoice
    // - ParamIDs::dialogDetectorChoice
    // - ParamIDs::dialogApplyChoice
    // - ParamIDs::lfeInDetector
    // - ParamIDs::lfeInApply
    // [END LVLR-STRUCTURAL-PARAM-SOURCE-OF-TRUTH]

    // [BEGIN MTDM-STRUCTURAL-PARAM-SOURCE-OF-TRUTH]
    // MTDM "structural / stop-playback / quality" parameter list.
    // Keep this comment block as the single source of truth so:
    //   1) APVTS non-automatable flags,
    //   2) UI playback-lock rules,
    //   3) latency-update bookkeeping
    // do not drift apart over time.
    //
    // Maintenance rule:
    // - If a param here should not be host-automatable, set .withAutomatable (false) on its APVTS parameter.
    // - If a param here should be blocked during playback, include it in the UI effective-playback lock list.
    // - If a param here changes reported latency, also update:
    //     cacheLatencyParamPointers()
    //     registerLatencyParamListeners()
    //     unregisterLatencyParamListeners()
    //     computeTotalLatencySamplesFromCachedParams()
    //
    // Structural enable params
    // (non-automatable + playback-locked; these may change processing topology/latency):
    // - ParamIDs::enabled              // mtdm.enabled
    // - ParamIDs::upEnabled01          // mtdm.up.enabled01
    // - ParamIDs::limEnabled01         // mtdm.lim.enabled
    //
    // Structural / quality params
    // (currently non-automatable + playback-locked):
    // - ParamIDs::upwardModeChoice     // mtdm.upwardModeChoice
    // - ParamIDs::sucFftSizeChoice     // mtdm.suc.fftSizeChoice
    // - ParamIDs::sucBandsPerOctChoice // mtdm.suc.bandsPerOctChoice
    // - ParamIDs::sucMinFreqHz         // mtdm.suc.minFreqHz
    // - ParamIDs::sucMaxFreqHz         // mtdm.suc.maxFreqHz
    // - ParamIDs::limLookaheadMs       // mtdm.lim.lookaheadMs
    // - ParamIDs::limOversamplingChoice// mtdm.lim.oversamplingChoice
    // [END MTDM-STRUCTURAL-PARAM-SOURCE-OF-TRUTH]

    // [BEGIN LVLR-APVTS-PARAM-LAYOUT]
    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { levelscope::lvlr::ParamIDs::enabled, 1 },
        "Leveler Enabled",
        (levelscope::lvlr::Defaults::enabled01 >= 0.5f)));

    layout.add (std::make_unique<juce::AudioParameterChoice> (
        juce::ParameterID { levelscope::lvlr::ParamIDs::controlModeChoice, 1 },
        "Leveler Control Mode",
        juce::StringArray { "Internal", "Host Gain" },
        levelscope::lvlr::Defaults::controlModeChoice,
        juce::AudioParameterChoiceAttributes().withAutomatable (false)));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { levelscope::lvlr::ParamIDs::hostGainDb, 1 },
        "Leveler Host Gain (dB)",
        juce::NormalisableRange<float> (levelscope::lvlr::Ranges::hostGainMinDb,
                                        levelscope::lvlr::Ranges::hostGainMaxDb,
                                        0.1f),
        levelscope::lvlr::Defaults::hostGainDb));

    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { levelscope::lvlr::ParamIDs::captureToHost01, 1 },
        "Leveler Capture To Host Gain",
        (levelscope::lvlr::Defaults::captureToHost01 >= 0.5f),
        juce::AudioParameterBoolAttributes().withAutomatable (false)));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { levelscope::lvlr::ParamIDs::targetLufs, 1 },
        "Leveler Target (LUFS)",
        juce::NormalisableRange<float> (levelscope::lvlr::Ranges::targetMinLufs,
                                        levelscope::lvlr::Ranges::targetMaxLufs,
                                        0.1f),
        levelscope::lvlr::Defaults::targetLufs));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { levelscope::lvlr::ParamIDs::maxBoostDb, 1 },
        "Leveler Max Boost (dB)",
        juce::NormalisableRange<float> (levelscope::lvlr::Ranges::maxBoostMinDb,
                                        levelscope::lvlr::Ranges::maxBoostMaxDb,
                                        0.1f),
        levelscope::lvlr::Defaults::maxBoostDb));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { levelscope::lvlr::ParamIDs::maxCutDb, 1 },
        "Leveler Max Cut (dB)",
        juce::NormalisableRange<float> (levelscope::lvlr::Ranges::maxCutMinDb,
                                        levelscope::lvlr::Ranges::maxCutMaxDb,
                                        0.1f),
        levelscope::lvlr::Defaults::maxCutDb));

    layout.add (std::make_unique<juce::AudioParameterChoice> (
        juce::ParameterID { levelscope::lvlr::ParamIDs::measChoice, 1 },
        "Leveler Measurement",
        juce::StringArray { "Auto", "Momentary", "Short-term" },
        levelscope::lvlr::Defaults::measChoice,
        juce::AudioParameterChoiceAttributes().withAutomatable (false)));

    layout.add (std::make_unique<juce::AudioParameterChoice> (
        juce::ParameterID { levelscope::lvlr::ParamIDs::modeChoice, 1 },
        "Leveler Algorithm Mode",
        juce::StringArray { "Continuous", "Learn-Hold" },
        levelscope::lvlr::Defaults::modeChoice,
        juce::AudioParameterChoiceAttributes().withAutomatable (false)));

    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { levelscope::lvlr::ParamIDs::learn01, 1 },
        "Leveler Learn",
        (levelscope::lvlr::Defaults::learn01 >= 0.5f),
        juce::AudioParameterBoolAttributes().withAutomatable (false)));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { levelscope::lvlr::ParamIDs::rateUpDbPerSec, 1 },
        "Leveler Rate Up (dB/s)",
        juce::NormalisableRange<float> (levelscope::lvlr::Ranges::rateUpMinDbPerSec,
                                        levelscope::lvlr::Ranges::rateUpMaxDbPerSec,
                                        0.1f),
        levelscope::lvlr::Defaults::rateUpDbPerSec));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { levelscope::lvlr::ParamIDs::rateDownDbPerSec, 1 },
        "Leveler Rate Down (dB/s)",
        juce::NormalisableRange<float> (levelscope::lvlr::Ranges::rateDownMinDbPerSec,
                                        levelscope::lvlr::Ranges::rateDownMaxDbPerSec,
                                        0.1f),
        levelscope::lvlr::Defaults::rateDownDbPerSec));
    // [END LVLR-APVTS-PARAM-LAYOUT]

    // [BEGIN MTDM-NONAUTOMATABLE-STRUCTURAL-ENABLES]
    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::enabled, 1 },
        "MTDM Enabled",
        (Defaults::enabled01 >= 0.5f),
        juce::AudioParameterBoolAttributes().withAutomatable (false)));
    // [END MTDM-NONAUTOMATABLE-STRUCTURAL-ENABLES]

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::thresholdDb, 1 },
        "MTDM Threshold (dB)",
        juce::NormalisableRange<float> (Ranges::thresholdMinDb, Ranges::thresholdMaxDb, 0.1f),
        Defaults::thresholdDb));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::ratio, 1 },
        "MTDM Ratio",
        juce::NormalisableRange<float> (Ranges::ratioMin, Ranges::ratioMax, 0.01f),
        Defaults::ratio));

    // [BEGIN MTDM-APVTS-PARAM-LAYOUT-STAGE-D1A-ADD]
    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::t0Lufs, 1 },
        "MTDM T0 (LUFS)",
        juce::NormalisableRange<float> (Ranges::t0MinLufs, Ranges::t0MaxLufs, 0.1f),
        Defaults::t0Lufs));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::t1Lufs, 1 },
        "MTDM T1 (LUFS)",
        juce::NormalisableRange<float> (Ranges::t1MinLufs, Ranges::t1MaxLufs, 0.1f),
        Defaults::t1Lufs));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::sucAmount01, 1 },
        "SUC Amount",
        juce::NormalisableRange<float> (Ranges::sucAmountMin01, Ranges::sucAmountMax01, 0.001f),
        Defaults::sucAmount01));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::sucMaxBoostDb, 1 },
        "SUC Max Boost (dB)",
        juce::NormalisableRange<float> (Ranges::sucMaxBoostMinDb, Ranges::sucMaxBoostMaxDb, 0.1f),
        Defaults::sucMaxBoostDb));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::sucCurve, 1 },
        "SUC Curve",
        juce::NormalisableRange<float> (Ranges::sucCurveMin, Ranges::sucCurveMax, 0.001f),
        Defaults::sucCurve));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::sucLowKneeDb, 1 },
        "SUC Low Knee (dB)",
        juce::NormalisableRange<float> (Ranges::sucKneeMinDb, Ranges::sucKneeMaxDb, 0.1f),
        Defaults::sucLowKneeDb));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::sucHighKneeDb, 1 },
        "SUC High Knee (dB)",
        juce::NormalisableRange<float> (Ranges::sucKneeMinDb, Ranges::sucKneeMaxDb, 0.1f),
        Defaults::sucHighKneeDb));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::sucAttackMs, 1 },
        "SUC Attack (ms)",
        juce::NormalisableRange<float> (Ranges::sucAttackMinMs, Ranges::sucAttackMaxMs, 0.1f),
        Defaults::sucAttackMs));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::sucReleaseMs, 1 },
        "SUC Release (ms)",
        juce::NormalisableRange<float> (Ranges::sucReleaseMinMs, Ranges::sucReleaseMaxMs, 0.1f),
        Defaults::sucReleaseMs));

    layout.add (std::make_unique<juce::AudioParameterChoice> (
        juce::ParameterID { ParamIDs::sucFftSizeChoice, 1 },
        "SUC FFT Size",
        juce::StringArray { "1024", "2048", "4096", "8192" },
        Defaults::sucFftSizeChoice,
        juce::AudioParameterChoiceAttributes().withAutomatable (false)));

    layout.add (std::make_unique<juce::AudioParameterChoice> (
        juce::ParameterID { ParamIDs::sucBandsPerOctChoice, 1 },
        "SUC Bands/Oct",
        juce::StringArray { "1", "2", "3", "6" },
        Defaults::sucBandsPerOctChoice,
        juce::AudioParameterChoiceAttributes().withAutomatable (false)));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::sucMinFreqHz, 1 },
        "SUC Min Freq (Hz)",
        juce::NormalisableRange<float> (Ranges::sucMinFreqMinHz, Ranges::sucMinFreqMaxHz, 1.0f),
        Defaults::sucMinFreqHz,
        juce::AudioParameterFloatAttributes().withAutomatable (false)));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::sucMaxFreqHz, 1 },
        "SUC Max Freq (Hz)",
        juce::NormalisableRange<float> (Ranges::sucMaxFreqMinHz, Ranges::sucMaxFreqMaxHz, 1.0f),
        Defaults::sucMaxFreqHz,
        juce::AudioParameterFloatAttributes().withAutomatable (false)));
    // [END MTDM-APVTS-PARAM-LAYOUT-STAGE-D1A-ADD]

    // [BEGIN MTDM-APVTS-PARAM-LAYOUT-TRIM-AND-CURVETYPE]
    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::sucCalTrimDb, 1 },
        "SUC Cal Trim (dB)",
        juce::NormalisableRange<float> (Ranges::sucCalTrimMinDb, Ranges::sucCalTrimMaxDb, 0.1f),
        Defaults::sucCalTrimDb));

    layout.add (std::make_unique<juce::AudioParameterChoice> (
        juce::ParameterID { ParamIDs::sucCurveTypeChoice, 1 },
        "SUC Curve Type",
        juce::StringArray { "Monotonic", "Bell" },
        Defaults::sucCurveTypeChoice));
    // [END MTDM-APVTS-PARAM-LAYOUT-TRIM-AND-CURVETYPE]

    // [BEGIN MTDM-APVTS-PARAM-LAYOUT-UPWARD-MODE]
    layout.add (std::make_unique<juce::AudioParameterChoice> (
        juce::ParameterID { ParamIDs::upwardModeChoice, 1 },
        "Upward Mode",
        juce::StringArray { "Spectral", "Broadband" },
        Defaults::upwardModeChoice,
        juce::AudioParameterChoiceAttributes().withAutomatable (false)));
    // [END MTDM-APVTS-PARAM-LAYOUT-UPWARD-MODE]

    // [BEGIN MTDM-APVTS-PARAM-LAYOUT-UPWARD-ENABLED-BYPASS]
    // [BEGIN MTDM-NONAUTOMATABLE-UPWARD-STRUCTURAL-ENABLE]
    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::upEnabled01, 1 },
        "Upward Enabled",
        (levelscope::mtdm::Defaults::upEnabled01 >= 0.5f),
        juce::AudioParameterBoolAttributes().withAutomatable (false)));
    // [END MTDM-NONAUTOMATABLE-UPWARD-STRUCTURAL-ENABLE]

    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::upBypass01, 1 },
        "Upward Bypass (safe)",
        (levelscope::mtdm::Defaults::upBypass01 >= 0.5f)));
    // [END MTDM-APVTS-PARAM-LAYOUT-UPWARD-ENABLED-BYPASS]

    // [BEGIN MTDM-APVTS-PARAM-LAYOUT-LFE-MASK]
    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::lfeInDetector, 1 },
        "LFE In Detector",
        (levelscope::mtdm::Defaults::lfeInDetector01 >= 0.5f)));

    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::lfeInApply, 1 },
        "LFE In Apply",
        (levelscope::mtdm::Defaults::lfeInApply01 >= 0.5f)));
    // [END MTDM-APVTS-PARAM-LAYOUT-LFE-MASK]

    // [BEGIN MTDM-APVTS-PARAM-LAYOUT-MC-POLICY]
    layout.add (std::make_unique<juce::AudioParameterChoice> (
        juce::ParameterID { ParamIDs::mcPolicyChoice, 1 },
        "Multichannel Policy",
        juce::StringArray { "Linked", "Dialog-mask", "Unlinked" },
        Defaults::mcPolicyChoice));

    layout.add (std::make_unique<juce::AudioParameterChoice> (
        juce::ParameterID { ParamIDs::dialogDetectorChoice, 1 },
        "Dialog Detector",
        juce::StringArray { "C", "LCR" },
        Defaults::dialogDetectorChoice));

    layout.add (std::make_unique<juce::AudioParameterChoice> (
        juce::ParameterID { ParamIDs::dialogApplyChoice, 1 },
        "Dialog Apply",
        juce::StringArray { "C", "LCR" },
        Defaults::dialogApplyChoice));
    // [END MTDM-APVTS-PARAM-LAYOUT-MC-POLICY]

    // [BEGIN MTDM-APVTS-PARAM-LAYOUT-DOWNWARD]
    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::t2Lufs, 1 },
        "MTDM T2 (LUFS)",
        juce::NormalisableRange<float> (Ranges::t2MinLufs, Ranges::t2MaxLufs, 0.1f),
        Defaults::t2Lufs));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::t3Lufs, 1 },
        "MTDM T3 (LUFS)",
        juce::NormalisableRange<float> (Ranges::t3MinLufs, Ranges::t3MaxLufs, 0.1f),
        Defaults::t3Lufs));

    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::downEnabled01, 1 },
        "Downward Enabled",
        (Defaults::downEnabled01 >= 0.5f)));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::downRatio, 1 },
        "Downward Ratio",
        juce::NormalisableRange<float> (Ranges::downRatioMin, Ranges::downRatioMax, 0.01f),
        Defaults::downRatio));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::downKneeDb, 1 },
        "Downward Knee (dB)",
        juce::NormalisableRange<float> (Ranges::downKneeMinDb, Ranges::downKneeMaxDb, 0.1f),
        Defaults::downKneeDb));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::downAttackMs, 1 },
        "Downward Attack (ms)",
        juce::NormalisableRange<float> (Ranges::downAttackMinMs, Ranges::downAttackMaxMs, 0.1f),
        Defaults::downAttackMs));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::downReleaseMs, 1 },
        "Downward Release (ms)",
        juce::NormalisableRange<float> (Ranges::downReleaseMinMs, Ranges::downReleaseMaxMs, 0.1f),
        Defaults::downReleaseMs));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::downMakeupDb, 1 },
        "Downward Makeup (dB)",
        juce::NormalisableRange<float> (Ranges::downMakeupMinDb, Ranges::downMakeupMaxDb, 0.1f),
        Defaults::downMakeupDb));
    // [END MTDM-APVTS-PARAM-LAYOUT-DOWNWARD]

    // [BEGIN MTDM-APVTS-PARAM-LAYOUT-DOWNWARD-BYPASS]
    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::downBypass01, 1 },
        "Downward Bypass (safe)",
        (levelscope::mtdm::Defaults::downBypass01 >= 0.5f)));
    // [END MTDM-APVTS-PARAM-LAYOUT-DOWNWARD-BYPASS]

    // [BEGIN MTDM-APVTS-PARAM-LAYOUT-LIMITER]
    // [BEGIN MTDM-NONAUTOMATABLE-LIMITER-STRUCTURAL-ENABLE]
    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::limEnabled01, 1 },
        "Limiter Enabled",
        (Defaults::limEnabled01 >= 0.5f),
        juce::AudioParameterBoolAttributes().withAutomatable (false)));
    // [END MTDM-NONAUTOMATABLE-LIMITER-STRUCTURAL-ENABLE]

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::limCeilingDb, 1 },
        "Limiter Ceiling (dBFS)",
        juce::NormalisableRange<float> (Ranges::limCeilingMinDb, Ranges::limCeilingMaxDb, 0.1f),
        Defaults::limCeilingDb));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::limLookaheadMs, 1 },
        "Limiter Lookahead (ms)",
        juce::NormalisableRange<float> (Ranges::limLookaheadMinMs, Ranges::limLookaheadMaxMs, 0.1f),
        Defaults::limLookaheadMs,
        juce::AudioParameterFloatAttributes().withAutomatable (false)));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::limReleaseMs, 1 },
        "Limiter Release (ms)",
        juce::NormalisableRange<float> (Ranges::limReleaseMinMs, Ranges::limReleaseMaxMs, 0.1f),
        Defaults::limReleaseMs));

    // [BEGIN MTDM-APVTS-PARAM-LAYOUT-LIMITER-TP-ATTACK-DRIVE]
    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::limAttackMs, 1 },
        "Limiter Attack (ms)",
        juce::NormalisableRange<float> (Ranges::limAttackMinMs, Ranges::limAttackMaxMs, 0.01f),
        Defaults::limAttackMs));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::limDriveDb, 1 },
        "Limiter Drive (dB)",
        juce::NormalisableRange<float> (Ranges::limDriveMinDb, Ranges::limDriveMaxDb, 0.1f),
        Defaults::limDriveDb));

    layout.add (std::make_unique<juce::AudioParameterChoice> (
        juce::ParameterID { ParamIDs::limOversamplingChoice, 1 },
        "Limiter Oversampling",
        juce::StringArray { "Off", "2x", "4x" },
        Defaults::limOversamplingChoice,
        juce::AudioParameterChoiceAttributes().withAutomatable (false)));
    // [END MTDM-APVTS-PARAM-LAYOUT-LIMITER-TP-ATTACK-DRIVE]
    // [END MTDM-APVTS-PARAM-LAYOUT-LIMITER]

    // [BEGIN MTDM-APVTS-PARAM-LAYOUT-LIMITER-BYPASS]
    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::limBypass01, 1 },
        "Limiter Bypass (safe)",
        (levelscope::mtdm::Defaults::limBypass01 >= 0.5f)));
    // [END MTDM-APVTS-PARAM-LAYOUT-LIMITER-BYPASS]

    // [BEGIN MTDM-APVTS-PARAM-LAYOUT-ZONE-SOLO-MUTE]
    layout.add (std::make_unique<juce::AudioParameterChoice> (
        juce::ParameterID { ParamIDs::zoneSoloChoice, 1 },
        "Zone Solo",
        juce::StringArray { "None", "Upward", "Downward", "Limiter", "Untouched" },
        Defaults::zoneSoloChoice));

    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::zoneUpwardMute01, 1 },
        "Upward Mute",
        (Defaults::zoneUpwardMute01 >= 0.5f)));

    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::zoneDownwardMute01, 1 },
        "Downward Mute",
        (Defaults::zoneDownwardMute01 >= 0.5f)));

    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::zoneLimiterMute01, 1 },
        "Limiter Mute",
        (Defaults::zoneLimiterMute01 >= 0.5f)));

    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::zoneUntouchedMute01, 1 },
        "Untouched Mute (audition)",
        (Defaults::zoneUntouchedMute01 >= 0.5f)));
    // [END MTDM-APVTS-PARAM-LAYOUT-ZONE-SOLO-MUTE]

    // [BEGIN MTDM-APVTS-PARAM-LAYOUT-ZONE-AUDITION]
    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::zoneAudBelowT0_01, 1 },
        "Zone Audition: Below T0",
        (Defaults::zoneAudBelowT0_01 >= 0.5f)));

    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::zoneAudT0T1_01, 1 },
        "Zone Audition: T0..T1",
        (Defaults::zoneAudT0T1_01 >= 0.5f)));

    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::zoneAudT1T2_01, 1 },
        "Zone Audition: T1..T2",
        (Defaults::zoneAudT1T2_01 >= 0.5f)));

    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::zoneAudT2T3_01, 1 },
        "Zone Audition: T2..T3",
        (Defaults::zoneAudT2T3_01 >= 0.5f)));

    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::zoneAudAboveT3_01, 1 },
        "Zone Audition: Above T3",
        (Defaults::zoneAudAboveT3_01 >= 0.5f)));
    // [END MTDM-APVTS-PARAM-LAYOUT-ZONE-AUDITION]

    return layout;
}
// [END MTDM-APVTS-PARAM-LAYOUT]

juce::int64 LevelScopeAudioProcessor::floorDivInt64 (juce::int64 a, juce::int64 b) noexcept
{
    if (b <= 0) return 0;
    if (a >= 0) return a / b;
    return - ( ( -a + b - 1 ) / b );
}

//==============================================================================

// [BEGIN MTDM-APVTS-CTOR-INIT]
LevelScopeAudioProcessor::LevelScopeAudioProcessor()
    : AudioProcessor (BusesProperties()
        .withInput  ("Input",  juce::AudioChannelSet::stereo(), true)
        .withOutput ("Output", juce::AudioChannelSet::stereo(), true)),
      apvts (*this, nullptr, "PARAMS", createParameterLayout())
// [END MTDM-APVTS-CTOR-INIT]
{
    // Phase 2A: core publishes LUFS (UI will be updated in Phase 2B)
    historyModel.setOutputMode (LevelScopeHistoryModel::OutputMode::lufs);

    // [BEGIN LS-PROCESSORCORE-CONSTRUCTOR-GRAPH-WITH-MTDM]
        // Stage C1/C2: default graph with MTDM + RT-safe APVTS binding.
        // Default: MTDM disabled => pass-through (no audible change).
    rebuildModuleGraphFromState (nullptr);
    // [END LS-PROCESSORCORE-CONSTRUCTOR-GRAPH-WITH-MTDM]

    // [BEGIN LS-LATENCY-CACHE-PARAMS-CTOR]
    cacheLatencyParamPointers();
    // [END LS-LATENCY-CACHE-PARAMS-CTOR]

    // [BEGIN LS-LATENCY-LISTENER-CTOR-START]
    registerLatencyParamListeners();
    // [END LS-LATENCY-LISTENER-CTOR-START]

    // [BEGIN LS-LATENCY-CTOR-UPDATE]
    updateLatencyFromAPVTS_NonRT();
    // [END LS-LATENCY-CTOR-UPDATE]
}

// [BEGIN LS-C2-MODGRAPH-IMPL]
void LevelScopeAudioProcessor::rebuildModuleGraphFromState (const juce::MemoryBlock* modgChunkData)
{
    // Non-audio-thread only. Builds a new graph snapshot and swaps it in.
    // [BEGIN LVLR-MODGRAPH-MODULE-IDS]
    const juce::String levelerId = juce::String (levelscope::LevelerModule().getModuleID());
    const juce::String mtdmId    = juce::String (levelscope::MultiThresholdDynamicsModule().getModuleID());
    // [END LVLR-MODGRAPH-MODULE-IDS]

    std::vector<juce::String> orderedModuleIds;
    orderedModuleIds.reserve (4);

    // [BEGIN LVLR-MODGRAPH-DEFAULT-ORDER]
    // Default graph if no MODG chunk (or invalid):
    // Leveler first, then MTDM.
    orderedModuleIds.push_back (levelerId);
    orderedModuleIds.push_back (mtdmId);
    // [END LVLR-MODGRAPH-DEFAULT-ORDER]

    if (modgChunkData != nullptr && modgChunkData->getSize() > 0)
    {
        juce::MemoryInputStream in (modgChunkData->getData(), modgChunkData->getSize(), false);

        const int schema = in.readInt();
        if (schema == 1)
        {
            const int numModules = in.readInt();
            if (numModules >= 0 && numModules < 1024)
            {
                orderedModuleIds.clear();
                orderedModuleIds.reserve ((size_t) numModules);

                for (int i = 0; i < numModules; ++i)
                {
                    if (in.getNumBytesRemaining() <= 0)
                        break;

                    const auto id = in.readString();
                    (void) in.readByte(); // bypass flag (stored, but not applied separately yet)
                    if (id.isNotEmpty())
                        orderedModuleIds.push_back (id);
                }

                if (orderedModuleIds.empty())
                {
                    orderedModuleIds.push_back (levelerId);
                    orderedModuleIds.push_back (mtdmId);
                }
            }
        }
    }

    // [BEGIN LVLR-MODGRAPH-ENSURE-ADDITIVE-INSERT]
    // Backward-compatible additive graph upgrade:
    // if an older session stored only MTDM in MODG, insert Leveler before it.
    {
        const auto levelerIt = std::find (orderedModuleIds.begin(), orderedModuleIds.end(), levelerId);
        const auto mtdmIt    = std::find (orderedModuleIds.begin(), orderedModuleIds.end(), mtdmId);

        if (levelerIt == orderedModuleIds.end())
        {
            if (mtdmIt != orderedModuleIds.end())
                orderedModuleIds.insert (mtdmIt, levelerId);
            else
                orderedModuleIds.insert (orderedModuleIds.begin(), levelerId);
        }

        if (std::find (orderedModuleIds.begin(), orderedModuleIds.end(), mtdmId) == orderedModuleIds.end())
            orderedModuleIds.push_back (mtdmId);
    }
    // [END LVLR-MODGRAPH-ENSURE-ADDITIVE-INSERT]

    auto graph = std::make_shared<levelscope::ModuleGraph>();
    graph->revision = 3;

    for (const auto& id : orderedModuleIds)
    {
        // [BEGIN LVLR-BINDPARAMS-CALL]
        if (id == levelerId)
        {
            auto leveler = std::make_shared<levelscope::LevelerModule>();

            leveler->bindParameters (
                apvts.getRawParameterValue (levelscope::lvlr::ParamIDs::enabled),
                apvts.getRawParameterValue (levelscope::lvlr::ParamIDs::targetLufs),
                apvts.getRawParameterValue (levelscope::lvlr::ParamIDs::maxBoostDb),
                apvts.getRawParameterValue (levelscope::lvlr::ParamIDs::maxCutDb),
                apvts.getRawParameterValue (levelscope::lvlr::ParamIDs::measChoice),
                apvts.getRawParameterValue (levelscope::lvlr::ParamIDs::modeChoice),
                apvts.getRawParameterValue (levelscope::lvlr::ParamIDs::learn01),
                apvts.getRawParameterValue (levelscope::lvlr::ParamIDs::rateUpDbPerSec),
                apvts.getRawParameterValue (levelscope::lvlr::ParamIDs::rateDownDbPerSec),
                apvts.getRawParameterValue (levelscope::lvlr::ParamIDs::controlModeChoice),
                apvts.getRawParameterValue (levelscope::lvlr::ParamIDs::hostGainDb),

                apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::mcPolicyChoice),
                apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::dialogDetectorChoice),
                apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::dialogApplyChoice),
                apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::lfeInDetector),
                apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::lfeInApply));

            // [BEGIN LS-LVLR-UI-HANDLE-STORE]
            std::atomic_store_explicit (&levelerForUI, leveler, std::memory_order_release);
            // [END LS-LVLR-UI-HANDLE-STORE]

            graph->modules.push_back (leveler);
        }
        // [END LVLR-BINDPARAMS-CALL]
        // [BEGIN MTDM-BINDPARAMS-CALL-FULL-D1C]
        if (id == mtdmId)
        {
            auto mtdm = std::make_shared<levelscope::MultiThresholdDynamicsModule>();

            mtdm->bindParameters (apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::enabled),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::thresholdDb),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::ratio),

                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::t0Lufs),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::t1Lufs),

                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::sucAmount01),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::sucMaxBoostDb),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::sucCurve),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::sucLowKneeDb),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::sucHighKneeDb),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::sucAttackMs),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::sucReleaseMs),

                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::sucFftSizeChoice),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::sucBandsPerOctChoice),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::sucMinFreqHz),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::sucMaxFreqHz),

                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::sucCalTrimDb),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::sucCurveTypeChoice),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::upwardModeChoice),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::lfeInDetector),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::lfeInApply),

                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::t2Lufs),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::t3Lufs),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::downEnabled01),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::downRatio),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::downKneeDb),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::downAttackMs),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::downReleaseMs),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::downMakeupDb),

                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::limEnabled01),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::limCeilingDb),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::limLookaheadMs),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::limReleaseMs),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::limAttackMs),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::limDriveDb),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::limOversamplingChoice),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::zoneSoloChoice),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::zoneUpwardMute01),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::zoneDownwardMute01),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::zoneLimiterMute01),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::zoneUntouchedMute01),
                                  // [BEGIN MTDM-BINDPARAMS-STAGE-E-MC]
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::mcPolicyChoice),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::dialogDetectorChoice),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::dialogApplyChoice),

                                  // [BEGIN MTDM-BINDPARAMS-ZONE-AUDITION]
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::zoneAudBelowT0_01),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::zoneAudT0T1_01),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::zoneAudT1T2_01),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::zoneAudT2T3_01),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::zoneAudAboveT3_01),

                                  // [BEGIN MTDM-BINDPARAMS-STAGE-BYPASS]
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::upEnabled01),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::upBypass01),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::downBypass01),
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::limBypass01));
                                  // [END MTDM-BINDPARAMS-STAGE-BYPASS]
                                  // [END MTDM-BINDPARAMS-ZONE-AUDITION]
                                  // [END MTDM-BINDPARAMS-STAGE-E-MC]

            // [BEGIN LS-MTDM-UI-HANDLE-STORE]
            std::atomic_store_explicit (&mtdmForUI, mtdm, std::memory_order_release);
            // [END LS-MTDM-UI-HANDLE-STORE]                  

            graph->modules.push_back (mtdm);
        }
        // [END MTDM-BINDPARAMS-CALL-FULL-D1C]
        // Unknown modules are ignored (forward-compat)
    }

    // [BEGIN LS-LVLR-UI-HANDLE-CLEAR-IF-NOT-SET]
    if (std::find (orderedModuleIds.begin(), orderedModuleIds.end(), levelerId) == orderedModuleIds.end())
        std::atomic_store_explicit (&levelerForUI, std::shared_ptr<levelscope::LevelerModule>{}, std::memory_order_release);
    // [END LS-LVLR-UI-HANDLE-CLEAR-IF-NOT-SET]

    // [BEGIN LS-MTDM-UI-HANDLE-CLEAR-IF-NOT-SET]
    if (std::find (orderedModuleIds.begin(), orderedModuleIds.end(), mtdmId) == orderedModuleIds.end())
        std::atomic_store_explicit (&mtdmForUI, std::shared_ptr<levelscope::MultiThresholdDynamicsModule>{}, std::memory_order_release);
    // [END LS-MTDM-UI-HANDLE-CLEAR-IF-NOT-SET]

    processorCore.setActiveGraph (std::move (graph));

    // If host already called prepareToPlay, immediately prepare the newly-swapped graph too.
    if (processorCorePrepared)
    {
        levelscope::ModulePrepareSpec spec;
        spec.sampleRate   = currentSampleRate;
        spec.maxBlockSize = lastMaxBlockSizeForSpec;
        spec.channelSet   = getBusesLayout().getMainInputChannelSet();

        processorCore.prepare (spec);
        processorCore.reset();
    }
}
// [END LS-C2-MODGRAPH-IMPL]

// [BEGIN LS-LATENCY-CACHED-PARAMS-IMPL]
void LevelScopeAudioProcessor::cacheLatencyParamPointers()
{
    using namespace levelscope::mtdm::ParamIDs;

    pLatMtdmEnabled01         = apvts.getRawParameterValue (enabled);
    pLatUpwardModeChoice      = apvts.getRawParameterValue (upwardModeChoice);
    pLatUpEnabled01           = apvts.getRawParameterValue (upEnabled01);
    pLatSucFftSizeChoice      = apvts.getRawParameterValue (sucFftSizeChoice);

    pLatLimEnabled01          = apvts.getRawParameterValue (limEnabled01);
    pLatLimLookaheadMs        = apvts.getRawParameterValue (limLookaheadMs);
    pLatLimOversamplingChoice = apvts.getRawParameterValue (limOversamplingChoice);
}

int LevelScopeAudioProcessor::computeTotalLatencySamplesFromCachedParams() const noexcept
{
    const bool mtdmEnabled =
        (pLatMtdmEnabled01 != nullptr && pLatMtdmEnabled01->load (std::memory_order_relaxed) >= 0.5f);

    if (! mtdmEnabled)
        return 0;

    // [BEGIN LS-LATENCY-UPWARD-RESPECT-UPENABLED]
    // --- Upward latency
    int upwardLatency = 0;

    const bool upEnabled =
        (pLatUpEnabled01 != nullptr && pLatUpEnabled01->load (std::memory_order_relaxed) >= 0.5f);

    const int upwardMode =
        (pLatUpwardModeChoice != nullptr ? (int) std::lround (pLatUpwardModeChoice->load (std::memory_order_relaxed))
                                         : levelscope::mtdm::Defaults::upwardModeChoice);

    if (upEnabled && upwardMode == 0) // Spectral
    {
        const int choice =
            (pLatSucFftSizeChoice != nullptr ? (int) std::lround (pLatSucFftSizeChoice->load (std::memory_order_relaxed))
                                             : levelscope::mtdm::Defaults::sucFftSizeChoice);

        const int fftSizes[] = { 1024, 2048, 4096, 8192 };
        upwardLatency = fftSizes[juce::jlimit (0, 3, choice)];
    }
    // [END LS-LATENCY-UPWARD-RESPECT-UPENABLED]

    // --- Limiter latency (lookahead + FIR detector latency)
    int limiterLatency = 0;

    const bool limEnabled =
        (pLatLimEnabled01 != nullptr && pLatLimEnabled01->load (std::memory_order_relaxed) >= 0.5f);

    if (limEnabled)
    {
        const float lookMs =
            (pLatLimLookaheadMs != nullptr ? pLatLimLookaheadMs->load (std::memory_order_relaxed)
                                           : levelscope::mtdm::Defaults::limLookaheadMs);

        const double sr = (currentSampleRate > 0.0 ? currentSampleRate : 48000.0);
        const double ms = juce::jlimit (0.0, 50.0, (double) lookMs);
        const int lookSamples = (int) std::lround (sr * ms * 0.001);

        int detDelay = 0;
        const int osChoice =
            (pLatLimOversamplingChoice != nullptr ? (int) std::lround (pLatLimOversamplingChoice->load (std::memory_order_relaxed))
                                                  : levelscope::mtdm::Defaults::limOversamplingChoice);

        const int clamped = juce::jlimit (0, 2, osChoice);
        if (clamped > 0)
        {
            const int stages = (clamped == 1 ? 1 : 2);
            // Using initProcessing(maxBlockSize) gives correct latency for the configured filter.
            juce::dsp::Oversampling<float> os (1, stages,
                juce::dsp::Oversampling<float>::filterHalfBandFIREquiripple, true);

            const int bs = juce::jmax (1, lastMaxBlockSizeForSpec);
            os.initProcessing ((size_t) bs);
            detDelay = os.getLatencyInSamples();
        }

        limiterLatency = lookSamples + detDelay;
    }

    return upwardLatency + limiterLatency;
}
// [END LS-LATENCY-CACHED-PARAMS-IMPL]

// [BEGIN LS-LATENCY-HELPER-IMPL]
void LevelScopeAudioProcessor::updateLatencyFromAPVTS_NonRT()
{
    const int totalLatency = computeTotalLatencySamplesFromCachedParams();
    setLatencySamples (totalLatency);
    updateHostDisplay();
}
// [END LS-LATENCY-HELPER-IMPL]

// [BEGIN LS-LATENCY-LISTENER-IMPL]

LevelScopeAudioProcessor::~LevelScopeAudioProcessor()
{
    stopTimer();
    endLevelerHostGainCaptureGesture_NonRT();
    unregisterLatencyParamListeners();
}

void LevelScopeAudioProcessor::registerLatencyParamListeners()
{
    using namespace levelscope::mtdm::ParamIDs;

    // Only params that can change total latency.
    apvts.addParameterListener (enabled, this);
    apvts.addParameterListener (upwardModeChoice, this);
    apvts.addParameterListener (upEnabled01, this);
    apvts.addParameterListener (sucFftSizeChoice, this);

    apvts.addParameterListener (limEnabled01, this);
    apvts.addParameterListener (limLookaheadMs, this);
    apvts.addParameterListener (limOversamplingChoice, this);

    // [BEGIN LVLR-HOSTGAIN-CAPTURE-TIMER-RATE]
    // 30 Hz message-thread timer:
    // - latency updates remain cheap
    // - transport freshness watchdog still works
    // - Leveler host-gain capture uses this cadence
    startTimerHz (30);
    // [END LVLR-HOSTGAIN-CAPTURE-TIMER-RATE]
}

void LevelScopeAudioProcessor::unregisterLatencyParamListeners()
{
    using namespace levelscope::mtdm::ParamIDs;

    apvts.removeParameterListener (enabled, this);
    apvts.removeParameterListener (upwardModeChoice, this);
    apvts.removeParameterListener (upEnabled01, this);
    apvts.removeParameterListener (sucFftSizeChoice, this);

    apvts.removeParameterListener (limEnabled01, this);
    apvts.removeParameterListener (limLookaheadMs, this);
    apvts.removeParameterListener (limOversamplingChoice, this);
}

// [BEGIN LS-LATENCY-LISTENER-PARAMCHANGED-IMMEDIATE]
void LevelScopeAudioProcessor::parameterChanged (const juce::String& parameterID, float newValue)
{
    juce::ignoreUnused (parameterID, newValue);

    // If this callback arrives on the message thread, update latency immediately
    // (avoids waiting for the 10 Hz timer, which can matter before an export).
    if (auto* mm = juce::MessageManager::getInstanceWithoutCreating())
    {
        if (mm->isThisTheMessageThread())
        {
            updateLatencyFromAPVTS_NonRT();
            updateHostDisplay();
            return;
        }
    }

    // RT-safe fallback (audio thread or unknown thread): mark dirty for timerCallback().
    latencyDirty.store (true, std::memory_order_release);
}
// [END LS-LATENCY-LISTENER-PARAMCHANGED-IMMEDIATE]

void LevelScopeAudioProcessor::timerCallback()
{
    // Transport freshness watchdog.
    {
        const auto seq = processBlockSequence.load (std::memory_order_relaxed);

        if (seq != lastObservedProcessBlockSequence)
        {
            lastObservedProcessBlockSequence = seq;
            processBlockQuietTimerTicks = 0;
            transportStateFresh.store (true, std::memory_order_relaxed);
        }
        else
        {
            ++processBlockQuietTimerTicks;

            if (processBlockQuietTimerTicks >= processBlockFreshnessTimeoutTicks)
                transportStateFresh.store (false, std::memory_order_relaxed);
        }
    }

    if (latencyDirty.exchange (false, std::memory_order_acq_rel))
    {
        updateLatencyFromAPVTS_NonRT();

        // Helps some hosts refresh; safe on message thread.
        updateHostDisplay();
    }

    // [BEGIN LVLR-HOSTGAIN-CAPTURE-TIMER]
    updateLevelerHostGainCapture_NonRT();
    // [END LVLR-HOSTGAIN-CAPTURE-TIMER]
}

// [BEGIN LVLR-HOSTGAIN-CAPTURE-IMPL]
void LevelScopeAudioProcessor::endLevelerHostGainCaptureGesture_NonRT()
{
    if (levelerHostGainGestureActive)
    {
        if (auto* hostGainParam = apvts.getParameter (levelscope::lvlr::ParamIDs::hostGainDb))
            hostGainParam->endChangeGesture();
    }

    levelerHostGainGestureActive = false;
    haveLastSentLevelerHostGainDb = false;
}

void LevelScopeAudioProcessor::forceDisarmLevelerHostGainCapture_NonRT()
{
    endLevelerHostGainCaptureGesture_NonRT();

    if (auto* captureParam = apvts.getParameter (levelscope::lvlr::ParamIDs::captureToHost01))
        captureParam->setValueNotifyingHost (0.0f);
}

void LevelScopeAudioProcessor::updateLevelerHostGainCapture_NonRT()
{
    const auto* controlModeRaw = apvts.getRawParameterValue (levelscope::lvlr::ParamIDs::controlModeChoice);
    const auto* captureRaw     = apvts.getRawParameterValue (levelscope::lvlr::ParamIDs::captureToHost01);

    const int controlMode =
        (controlModeRaw != nullptr ? (int) std::lround (controlModeRaw->load (std::memory_order_relaxed))
                                   : levelscope::lvlr::Defaults::controlModeChoice);

    const bool captureArmed =
        (captureRaw != nullptr && captureRaw->load (std::memory_order_relaxed) >= 0.5f);

    const bool shouldCapture =
        (controlMode == levelscope::LevelerModule::controlInternal) && captureArmed;

    auto* hostGainParam = apvts.getParameter (levelscope::lvlr::ParamIDs::hostGainDb);

    if (! shouldCapture || hostGainParam == nullptr)
    {
        endLevelerHostGainCaptureGesture_NonRT();
        return;
    }

    if (! levelerHostGainGestureActive)
    {
        hostGainParam->beginChangeGesture();
        levelerHostGainGestureActive = true;
        haveLastSentLevelerHostGainDb = false;
    }

    const auto snap = getLevelerMeteringSnapshot();
    const float capturedDb = juce::jlimit (levelscope::lvlr::Ranges::hostGainMinDb,
                                           levelscope::lvlr::Ranges::hostGainMaxDb,
                                           snap.gainDbCurrent);

    if (! haveLastSentLevelerHostGainDb
        || std::abs (capturedDb - lastSentLevelerHostGainDb) >= 0.1f)
    {
        hostGainParam->setValueNotifyingHost (hostGainParam->convertTo0to1 (capturedDb));
        lastSentLevelerHostGainDb = capturedDb;
        haveLastSentLevelerHostGainDb = true;
    }
}
// [END LVLR-HOSTGAIN-CAPTURE-IMPL]
// [END LS-LATENCY-LISTENER-IMPL]

// [BEGIN LS-UPWARD-METERING-SNAPSHOT-IMPL]
LevelScopeAudioProcessor::UpwardMeteringSnapshot LevelScopeAudioProcessor::getUpwardMeteringSnapshot() const noexcept
{
    UpwardMeteringSnapshot s;

    auto m = std::atomic_load_explicit (&mtdmForUI, std::memory_order_acquire);
    if (! m)
        return s;

    const auto& met = m->getUpwardMetering();
    s.boostDbCurrent   = met.boostDbCurrent.load   (std::memory_order_relaxed);
    s.boostDbBlockPeak = met.boostDbBlockPeak.load (std::memory_order_relaxed);
    s.boostDbHold      = met.boostDbHold.load      (std::memory_order_relaxed);
    return s;
}
// [END LS-UPWARD-METERING-SNAPSHOT-IMPL]

// [BEGIN LS-LVLR-METERING-SNAPSHOT-IMPL]
LevelScopeAudioProcessor::LevelerMeteringSnapshot LevelScopeAudioProcessor::getLevelerMeteringSnapshot() const noexcept
{
    LevelerMeteringSnapshot s;

    auto m = std::atomic_load_explicit (&levelerForUI, std::memory_order_acquire);
    if (! m)
        return s;

    const auto& met = m->getMetering();
    s.gainDbCurrent   = met.gainDbCurrent.load   (std::memory_order_relaxed);
    s.gainDbBlockPeak = met.gainDbBlockPeak.load (std::memory_order_relaxed);
    s.gainDbHold      = met.gainDbHold.load      (std::memory_order_relaxed);
    return s;
}
// [END LS-LVLR-METERING-SNAPSHOT-IMPL]

// [BEGIN LS-LIMITER-METERING-SNAPSHOT-IMPL]
LevelScopeAudioProcessor::LimiterMeteringSnapshot LevelScopeAudioProcessor::getLimiterMeteringSnapshot() const noexcept
{
    LimiterMeteringSnapshot s;

    auto m = std::atomic_load_explicit (&mtdmForUI, std::memory_order_acquire);
    if (! m)
        return s;

    const auto& met = m->getLimiterMetering();
    s.grDbCurrent   = met.grDbCurrent.load   (std::memory_order_relaxed);
    s.grDbBlockPeak = met.grDbBlockPeak.load (std::memory_order_relaxed);
    s.grDbHold      = met.grDbHold.load      (std::memory_order_relaxed);
    return s;
}
// [END LS-LIMITER-METERING-SNAPSHOT-IMPL]

// [BEGIN LS-DOWNWARD-METERING-SNAPSHOT-IMPL]
LevelScopeAudioProcessor::DownwardMeteringSnapshot LevelScopeAudioProcessor::getDownwardMeteringSnapshot() const noexcept
{
    DownwardMeteringSnapshot s;

    auto m = std::atomic_load_explicit (&mtdmForUI, std::memory_order_acquire);
    if (! m)
        return s;

    const auto& met = m->getDownwardMetering();
    s.grDbCurrent   = met.grDbCurrent.load   (std::memory_order_relaxed);
    s.grDbBlockPeak = met.grDbBlockPeak.load (std::memory_order_relaxed);
    s.grDbHold      = met.grDbHold.load      (std::memory_order_relaxed);
    return s;
}
// [END LS-DOWNWARD-METERING-SNAPSHOT-IMPL]

// [BEGIN LS-IO-METERING-SNAPSHOT-IMPL]
LevelScopeAudioProcessor::IOMeteringSnapshot LevelScopeAudioProcessor::getIOMeteringSnapshot() const noexcept
{
    IOMeteringSnapshot s;
    s.inPeakDbCurrent  = ioInPeakDbCurrent.load  (std::memory_order_relaxed);
    s.inPeakDbHold     = ioInPeakDbHold.load     (std::memory_order_relaxed);
    s.inRmsDbCurrent   = ioInRmsDbCurrent.load   (std::memory_order_relaxed);

    s.outPeakDbCurrent = ioOutPeakDbCurrent.load (std::memory_order_relaxed);
    s.outPeakDbHold    = ioOutPeakDbHold.load    (std::memory_order_relaxed);
    s.outRmsDbCurrent  = ioOutRmsDbCurrent.load  (std::memory_order_relaxed);
    return s;
}
// [END LS-IO-METERING-SNAPSHOT-IMPL]

//==============================================================================

// [BEGIN LS-UIST-IMPL]
LevelScopeAudioProcessor::UIStateSnapshot LevelScopeAudioProcessor::getPersistedUIStateSnapshot() const noexcept
{
    UIStateSnapshot s;
    s.bottomPanelFraction01 = uiBottomPanelFraction01.load (std::memory_order_relaxed);

    s.rightStripWidthPxUser = uiRightStripWidthPxUser.load (std::memory_order_relaxed);
    s.rollingLaneHeightPx   = uiRollingLaneHeightPx.load   (std::memory_order_relaxed);

    s.cardHeights.levelling = uiCardLevellingPx.load (std::memory_order_relaxed);
    s.cardHeights.zones     = uiCardZonesPx.load     (std::memory_order_relaxed);
    s.cardHeights.audition  = uiCardAuditionPx.load  (std::memory_order_relaxed);
    s.cardHeights.upward    = uiCardUpwardPx.load    (std::memory_order_relaxed);
    s.cardHeights.downward  = uiCardDownwardPx.load  (std::memory_order_relaxed);
    s.cardHeights.limiter   = uiCardLimiterPx.load   (std::memory_order_relaxed);

    s.upwardAdvancedOpen   = (uiUpwardAdvancedOpen01.load   (std::memory_order_relaxed) != 0);
    s.downwardAdvancedOpen = (uiDownwardAdvancedOpen01.load (std::memory_order_relaxed) != 0);
    s.limiterAdvancedOpen  = (uiLimiterAdvancedOpen01.load  (std::memory_order_relaxed) != 0);

    s.showMomentaryCurve = (uiShowMomentaryCurve01.load (std::memory_order_relaxed) != 0);
    s.showShortTermCurve = (uiShowShortTermCurve01.load (std::memory_order_relaxed) != 0);
    s.showGateCurve      = (uiShowGateCurve01.load      (std::memory_order_relaxed) != 0);
    s.showRollingLra     = (uiShowRollingLra01.load     (std::memory_order_relaxed) != 0);
    s.followRightEdge    = (uiFollowRightEdge01.load    (std::memory_order_relaxed) != 0);

    s.rollingWindowSeconds = uiRollingWindowSeconds.load (std::memory_order_relaxed);

    s.historyViewStateValid = (uiHistoryViewStateValid01.load (std::memory_order_relaxed) != 0);
    s.historyHasCustomZoomX = (uiHistoryHasCustomZoomX01.load (std::memory_order_relaxed) != 0);
    s.historyZoomX          = uiHistoryZoomX.load          (std::memory_order_relaxed);
    s.historyZoomY          = uiHistoryZoomY.load          (std::memory_order_relaxed);
    s.historyViewRightFrame = uiHistoryViewRightFrame.load (std::memory_order_relaxed);
    s.historyViewTopDb      = uiHistoryViewTopDb.load      (std::memory_order_relaxed);

    return s;
}

void LevelScopeAudioProcessor::setUiBottomPanelFraction01 (float v) noexcept
{
    v = juce::jlimit (0.10f, 0.60f, v);
    uiBottomPanelFraction01.store (v, std::memory_order_relaxed);
}

void LevelScopeAudioProcessor::setUiRightStripWidthPx (int v) noexcept
{
    v = juce::jlimit (112, 300, v);
    uiRightStripWidthPxUser.store (v, std::memory_order_relaxed);
}

void LevelScopeAudioProcessor::setUiRollingLaneHeightPx (int v) noexcept
{
    v = juce::jlimit (28, 240, v);
    uiRollingLaneHeightPx.store (v, std::memory_order_relaxed);
}

void LevelScopeAudioProcessor::setUiCardHeights (const UICardHeightsState& h) noexcept
{
    uiCardLevellingPx.store (juce::jmax (50, h.levelling), std::memory_order_relaxed);
    uiCardZonesPx.store     (juce::jmax (34, h.zones),     std::memory_order_relaxed);
    uiCardAuditionPx.store  (juce::jmax (34, h.audition),  std::memory_order_relaxed);
    uiCardUpwardPx.store    (juce::jmax (34, h.upward),    std::memory_order_relaxed);
    uiCardDownwardPx.store  (juce::jmax (26, h.downward),  std::memory_order_relaxed);
    uiCardLimiterPx.store   (juce::jmax (26, h.limiter),   std::memory_order_relaxed);
}

void LevelScopeAudioProcessor::setUiUpwardAdvancedOpen (bool b) noexcept
{
    uiUpwardAdvancedOpen01.store (b ? 1 : 0, std::memory_order_relaxed);
}

void LevelScopeAudioProcessor::setUiDownwardAdvancedOpen (bool b) noexcept
{
    uiDownwardAdvancedOpen01.store (b ? 1 : 0, std::memory_order_relaxed);
}

void LevelScopeAudioProcessor::setUiLimiterAdvancedOpen (bool b) noexcept
{
    uiLimiterAdvancedOpen01.store (b ? 1 : 0, std::memory_order_relaxed);
}

void LevelScopeAudioProcessor::setUiHistoryToggleState (bool showMomentary,
                                                        bool showShortTerm,
                                                        bool showGate,
                                                        bool showRollingLra,
                                                        bool followRightEdge) noexcept
{
    uiShowMomentaryCurve01.store (showMomentary ? 1 : 0, std::memory_order_relaxed);
    uiShowShortTermCurve01.store (showShortTerm ? 1 : 0, std::memory_order_relaxed);
    uiShowGateCurve01.store      (showGate ? 1 : 0, std::memory_order_relaxed);
    uiShowRollingLra01.store     (showRollingLra ? 1 : 0, std::memory_order_relaxed);
    uiFollowRightEdge01.store    (followRightEdge ? 1 : 0, std::memory_order_relaxed);
}

void LevelScopeAudioProcessor::setUiHistoryViewState (double zoomX,
                                                      double zoomY,
                                                      double viewRightFrame,
                                                      double viewTopDb,
                                                      bool hasCustomZoomX,
                                                      bool valid) noexcept
{
    uiHistoryZoomX.store          (zoomX,          std::memory_order_relaxed);
    uiHistoryZoomY.store          (zoomY,          std::memory_order_relaxed);
    uiHistoryViewRightFrame.store (viewRightFrame, std::memory_order_relaxed);
    uiHistoryViewTopDb.store      (viewTopDb,      std::memory_order_relaxed);

    uiHistoryHasCustomZoomX01.store (hasCustomZoomX ? 1 : 0, std::memory_order_relaxed);
    uiHistoryViewStateValid01.store (valid ? 1 : 0, std::memory_order_relaxed);
}

void LevelScopeAudioProcessor::buildUistStateChunk (juce::MemoryBlock& destData) const
{
    destData.setSize (0);

    const auto s = getPersistedUIStateSnapshot();

    juce::MemoryOutputStream os (destData, true);
    os.writeInt (2); // UIST schema v2

    // v1 fields
    os.writeFloat (juce::jlimit (0.10f, 0.60f, s.bottomPanelFraction01));
    os.writeInt   (juce::jlimit (112, 300, s.rightStripWidthPxUser));
    os.writeInt   (juce::jlimit (28, 240, s.rollingLaneHeightPx));

    os.writeInt (juce::jmax (50, s.cardHeights.levelling));
    os.writeInt (juce::jmax (34, s.cardHeights.zones));
    os.writeInt (juce::jmax (34, s.cardHeights.audition));
    os.writeInt (juce::jmax (34, s.cardHeights.upward));
    os.writeInt (juce::jmax (26, s.cardHeights.downward));
    os.writeInt (juce::jmax (26, s.cardHeights.limiter));

    os.writeByte ((char) (s.upwardAdvancedOpen   ? 1 : 0));
    os.writeByte ((char) (s.downwardAdvancedOpen ? 1 : 0));
    os.writeByte ((char) (s.limiterAdvancedOpen  ? 1 : 0));
    os.writeByte ((char) 0); // reserved/padding

    // v2 fields
    os.writeByte ((char) (s.showMomentaryCurve ? 1 : 0));
    os.writeByte ((char) (s.showShortTermCurve ? 1 : 0));
    os.writeByte ((char) (s.showGateCurve      ? 1 : 0));
    os.writeByte ((char) (s.showRollingLra     ? 1 : 0));
    os.writeByte ((char) (s.followRightEdge    ? 1 : 0));
    os.writeByte ((char) (s.historyViewStateValid ? 1 : 0));
    os.writeByte ((char) (s.historyHasCustomZoomX ? 1 : 0));
    os.writeByte ((char) 0); // reserved

    const int rollingWindow = (s.rollingWindowSeconds == 30 || s.rollingWindowSeconds == 60 || s.rollingWindowSeconds == 120
                                ? s.rollingWindowSeconds : 60);
    os.writeInt (rollingWindow);

    os.writeDouble (s.historyZoomX);
    os.writeDouble (s.historyZoomY);
    os.writeDouble (s.historyViewRightFrame);
    os.writeDouble (s.historyViewTopDb);
}

bool LevelScopeAudioProcessor::applyUistStateChunk (const juce::MemoryBlock& uistChunk)
{
    constexpr size_t minBytesV1 = 40;
    constexpr size_t minBytesV2 = 88;

    if (uistChunk.getSize() < minBytesV1)
        return false;

    juce::MemoryInputStream in (uistChunk.getData(), uistChunk.getSize(), false);
    const int schema = in.readInt();

    if (schema != 1 && schema != 2)
        return false;

    setUiBottomPanelFraction01 (in.readFloat());
    setUiRightStripWidthPx     (in.readInt());
    setUiRollingLaneHeightPx   (in.readInt());

    UICardHeightsState h;
    h.levelling = in.readInt();
    h.zones     = in.readInt();
    h.audition  = in.readInt();
    h.upward    = in.readInt();
    h.downward  = in.readInt();
    h.limiter   = in.readInt();
    setUiCardHeights (h);

    setUiUpwardAdvancedOpen   (in.readByte() != 0);
    setUiDownwardAdvancedOpen (in.readByte() != 0);
    setUiLimiterAdvancedOpen  (in.readByte() != 0);
    (void) in.readByte(); // reserved

    if (schema == 1)
        return true;

    if (uistChunk.getSize() < minBytesV2)
        return false;

    const bool showMomentary = (in.readByte() != 0);
    const bool showShortTerm = (in.readByte() != 0);
    const bool showGate      = (in.readByte() != 0);
    const bool showRolling   = (in.readByte() != 0);
    const bool follow        = (in.readByte() != 0);
    const bool viewValid     = (in.readByte() != 0);
    const bool customZoomX   = (in.readByte() != 0);
    (void) in.readByte(); // reserved

    setUiHistoryToggleState (showMomentary, showShortTerm, showGate, showRolling, follow);
    setRollingLraWindowSeconds (in.readInt());

    setUiHistoryViewState (in.readDouble(),
                           in.readDouble(),
                           in.readDouble(),
                           in.readDouble(),
                           customZoomX,
                           viewValid);

    return true;
}
// [END LS-UIST-IMPL]

// [BEGIN LS-PROCESSORCORE-PREPARETOPLAY]
void LevelScopeAudioProcessor::prepareToPlay (double sampleRate,
                                             int samplesPerBlockExpected)
{
    currentSampleRate = (sampleRate > 0.0 ? sampleRate : 44100.0);

    frameSamples = juce::jmax (1,
        (int) std::round (currentSampleRate / loudnessFrameRate));

    historyModel.setFrameSamplesForMetadata (frameSamples);

    // [BS1770] prepare K-weighting for current channel count
    const int numCh = juce::jmax (1, getTotalNumInputChannels());
    kWeight.prepare (currentSampleRate, numCh);

    // [BS1770] channel weights (LFE = 0, everything else = 1)
    bs1770ChannelWeights.assign ((size_t) numCh, 1.0f);
    const auto layout = getBusesLayout().getMainInputChannelSet();
    for (int ch = 0; ch < numCh; ++ch)
        if (layout.getTypeOfChannel (ch) == juce::AudioChannelSet::LFE)
            bs1770ChannelWeights[(size_t) ch] = 0.0f;

    // Stage A: prepare ProcessorCore (empty graph => no-op).
    // We avoid juce::dsp::ProcessSpec (juce_dsp not enabled).
    // Also, we don't assume ModulePrepareSpec field names here.
    // [BEGIN LS-PROCESSORCORE-PREPARE-SPEC]
    levelscope::ModulePrepareSpec spec;
    spec.sampleRate  = currentSampleRate;
    spec.maxBlockSize = juce::jmax (0, samplesPerBlockExpected);
    spec.channelSet  = getBusesLayout().getMainInputChannelSet();

    processorCore.prepare (spec);
    // [END LS-PROCESSORCORE-PREPARE-SPEC]

    // [BEGIN LS-C2-PREPARE-FLAGS]
        processorCorePrepared = true;
        lastMaxBlockSizeForSpec = juce::jmax (0, samplesPerBlockExpected);
    // [END LS-C2-PREPARE-FLAGS]

    // [BEGIN LS-LATENCY-PREPARE-UPDATE]
        updateLatencyFromAPVTS_NonRT();
    // [END LS-LATENCY-PREPARE-UPDATE]

    resetLoudnessState();
}
// [END LS-PROCESSORCORE-PREPARETOPLAY]

void LevelScopeAudioProcessor::releaseResources()
{
    resetLoudnessState();

    // [BEGIN LS-C2-RELEASE-FLAGS]
    processorCorePrepared = false;
    // [END LS-C2-RELEASE-FLAGS]
}

// [BEGIN LS-PROCESSORCORE-RESETLOUDNESSSTATE]
void LevelScopeAudioProcessor::resetLoudnessState() noexcept
{
    samplesUntilNextFrame = frameSamples;
    frameEnergyAccum = 0.0;
    historyModel.resetRealtimeFifo();
    kWeight.reset();
    runningStats.reset();
    transportIsPlaying.store (false, std::memory_order_relaxed);
    transportStateFresh.store (false, std::memory_order_relaxed);
    transportIsPlaying.store (false, std::memory_order_relaxed);

    // RT-safe reset; Stage A graph is empty => no audible change
    processorCore.reset();
}
// [END LS-PROCESSORCORE-RESETLOUDNESSSTATE]

//==============================================================================

// [BEGIN LS-BUSES-LAYOUT-SUPPORTED-STAGE-E]
bool LevelScopeAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
    const auto& mainIn  = layouts.getMainInputChannelSet();
    const auto& mainOut = layouts.getMainOutputChannelSet();

    if (mainIn.isDisabled() || mainOut.isDisabled())
        return false;

    // Must be symmetrical I/O for this processor (same channel set).
    if (mainIn != mainOut)
        return false;

    // Disallow ambisonics for now (channel-role policies assume speaker roles).
    // [BEGIN LS-BUSES-AMBISONIC-CHECK-COMPAT]
    // JUCE compatibility: some JUCE versions don't provide AudioChannelSet::isAmbisonic().
    // Heuristic is fine here (non-RT path): ambisonic layouts are labeled "Ambisonic" in description.
    if (mainIn.getDescription().containsIgnoreCase ("ambisonic"))
        return false;
    // [END LS-BUSES-AMBISONIC-CHECK-COMPAT]

    const int ch = mainIn.size();

    // Allow up to 7.1.4 (12ch). (We keep headroom for future by changing this constant later.)
    if (ch <= 0 || ch > 12)
        return false;

    return true;
}
// [END LS-BUSES-LAYOUT-SUPPORTED-STAGE-E]

//==============================================================================

void LevelScopeAudioProcessor::processSampleForLoudness (const float* const* channelData,
                                                        int numChannels,
                                                        int sampleIndex) noexcept
{
    // [BS1770] K-weighted mean-square energy (sum over channels; LFE weight=0)
    double e = 0.0;

    for (int ch = 0; ch < numChannels; ++ch)
    {
        const float w = (ch < (int) bs1770ChannelWeights.size() ? bs1770ChannelWeights[(size_t) ch] : 1.0f);
        if (w <= 0.0f)
            continue;

        const float x = channelData[ch][sampleIndex];
        const float y = kWeight.processSample (ch, x);

        e += (double) w * (double) y * (double) y;
    }

    frameEnergyAccum += e;

    // 60 Hz frame scheduling
    if (--samplesUntilNextFrame <= 0)
    {
        samplesUntilNextFrame += frameSamples;

        const double meanSquare = frameEnergyAccum / (double) juce::jmax (1, frameSamples);
        frameEnergyAccum = 0.0;

        if (skipNextPartialFrameWrite)
        {
            skipNextPartialFrameWrite = false;
            return;
        }

        const juce::int64 frameProjectSample =
            currentBlockStartProjectSample + (juce::int64) currentBlockSampleIndex;

        const juce::int64 frameIndex =
            floorDivInt64 (frameProjectSample, (juce::int64) frameSamples);

        // [FIX-START-RAMP] do not overwrite existing timeline truth during guard period
        if (transportStartOverwriteGuardFrames > 0)
        {
            const bool alreadyExists = historyModel.frameExists (frameIndex);
            --transportStartOverwriteGuardFrames;

            if (alreadyExists)
                return;
        }

        const float energyMS = (float) juce::jmax (0.0, meanSquare);

        const int momentaryFrames = (int) std::round (momentaryWindowSeconds * loudnessFrameRate); // 24
        const int shortTermFrames = (int) std::round (shortTermWindowSeconds * loudnessFrameRate); // 180

        // Update running stats first (integrated + LRA, etc.)
        runningStats.pushFrameEnergy (energyMS);

        // LRA "relative gate" curve value (integratedRunning - 20 LU)
        float gateLufs = runningStats.getLraGateLufs();

        // [LRAG-WARMUP-OVERWRITE-FIX]
        // After a stop/start, integrated is not valid immediately (warm-up).
        // If we're overwriting an already-existing timeline frame during warm-up,
        // do NOT overwrite its previously stored gate curve value.
        const float I = runningStats.getIntegratedLufs();
        const bool integratedValid = (I > -199.0f);

        if (! integratedValid && historyModel.frameExists (frameIndex))
        {
            float existingGate = -200.0f;
            if (historyModel.getLraGateLufsAtFrameIndex (frameIndex, existingGate))
                gateLufs = existingGate;
        }

        historyModel.pushEnergyFrame (frameIndex,
                                     energyMS,
                                     momentaryFrames,
                                     shortTermFrames,
                                     currentBlockIsPlaying,
                                     gateLufs);
    }
}

void LevelScopeAudioProcessor::processBlock (juce::AudioBuffer<float>& buffer,
                                            juce::MidiBuffer& midiMessages)
{
    juce::ScopedNoDenormals noDenormals;
    juce::ignoreUnused (midiMessages);

    const int numChannels = getTotalNumInputChannels();
    const int numSamples  = buffer.getNumSamples();

    processBlockSequence.fetch_add (1u, std::memory_order_relaxed);
    transportStateFresh.store (true, std::memory_order_relaxed);

    // [BEGIN LS-IO-METERING-IN-COMPUTE]
    // Compute INPUT metering before ProcessorCore modifies the buffer.
    const auto inMet = computeBlockPeakRmsMaxAcrossChannels (buffer, numChannels, numSamples);

    ioInPeakDbCurrent.store (inMet.peakDb, std::memory_order_relaxed);
    ioInRmsDbCurrent.store  (inMet.rmsDb,  std::memory_order_relaxed);

    {
        const float hold = ioInPeakDbHold.load (std::memory_order_relaxed);
        const float dt = (currentSampleRate > 1.0e-6 ? (float) numSamples / (float) currentSampleRate : 0.0f);
        const float dec = ioPeakHoldDecayDbPerSec * dt;
        ioInPeakDbHold.store (decayHoldDb (hold, inMet.peakDb, dec), std::memory_order_relaxed);
    }
    // [END LS-IO-METERING-IN-COMPUTE]

    for (int ch = numChannels; ch < getTotalNumOutputChannels(); ++ch)
        buffer.clear (ch, 0, numSamples);

    if (numChannels <= 0 || numSamples <= 0)
        return;

    // [TIMEBASE-PLAYHEAD]
    juce::int64 blockStartProjectSample = 0;
    int blockIsPlaying = 0;
    bool haveTimeInSamples = false;
    bool haveTimeInSeconds = false;

    juce::int64 hostSamplesValue = 0;
    double hostSecondsValue = 0.0;

    if (auto* ph = getPlayHead())
    {
        if (auto pos = ph->getPosition())
        {
            const auto tSamp = pos->getTimeInSamples();
            if (optHasValue (tSamp))
            {
                hostSamplesValue = (juce::int64) optValue (tSamp);
                haveTimeInSamples = true;
                blockStartProjectSample = hostSamplesValue;
            }

            const auto tSec = pos->getTimeInSeconds();
            if (optHasValue (tSec))
            {
                hostSecondsValue = (double) optValue (tSec);
                haveTimeInSeconds = true;
            }

            const auto pl = pos->getIsPlaying();
            if (optHasValue (pl))
                blockIsPlaying = (optValue (pl) ? 1 : 0);

            haveHostTimeSamples.store (haveTimeInSamples ? 1 : 0, std::memory_order_relaxed);
            haveHostTimeSeconds.store (haveTimeInSeconds ? 1 : 0, std::memory_order_relaxed);
            if (haveTimeInSamples)
                lastHostTimeSamples.store (hostSamplesValue, std::memory_order_relaxed);
            if (haveTimeInSeconds)
                lastHostTimeSeconds.store (hostSecondsValue, std::memory_order_relaxed);

            if (haveTimeInSamples && haveTimeInSeconds && currentSampleRate > 1.0e-12)
            {
                const double sampleSeconds = (double) hostSamplesValue / currentSampleRate;
                timecodeOffsetSeconds.store (hostSecondsValue - sampleSeconds, std::memory_order_relaxed);
            }
        }
    }
    transportIsPlaying.store (blockIsPlaying == 1, std::memory_order_relaxed);

    const bool shouldAnalyse = (blockIsPlaying == 1 && haveTimeInSamples);

    // [FIX-RESTART-PARTIAL-FRAME]
    const bool discontinuity =
        (! haveLastBlockEnd) ||
        (lastBlockIsPlaying == 0) ||
        (blockStartProjectSample != lastBlockEndProjectSample);

    if (shouldAnalyse && discontinuity && frameSamples > 0)
    {
        frameEnergyAccum = 0.0;
        kWeight.reset();
        runningStats.reset();

        juce::int64 mod = blockStartProjectSample % (juce::int64) frameSamples;
        if (mod < 0) mod += (juce::int64) frameSamples;

        samplesUntilNextFrame = frameSamples - (int) mod;
        if (samplesUntilNextFrame <= 0)
            samplesUntilNextFrame += frameSamples;

        skipNextPartialFrameWrite = (mod != 0);
    }

    const bool transportStart = (blockIsPlaying == 1 && lastBlockIsPlaying == 0);
    if (transportStart)
        transportStartOverwriteGuardFrames = 6;

    haveLastBlockEnd = true;
    lastBlockEndProjectSample = blockStartProjectSample + (juce::int64) numSamples;
    lastBlockIsPlaying = blockIsPlaying;

    // [BEGIN LS-PROCESSORCORE-SKIP-ANALYSIS-STILL-PROCESS]
        if (! shouldAnalyse)
        {
            // Stage A: still run ProcessorCore (empty graph => no-op, buffer unchanged).
            levelscope::ProcessContext ctx { buffer,
                                            &midiMessages,
                                            currentSampleRate,
                                            numSamples,
                                            getBusesLayout().getMainInputChannelSet() };

            ctx.isRealtime        = ! isNonRealtime();
            ctx.isPlaying         = (blockIsPlaying == 1);
            ctx.isDiscontinuity   = discontinuity;
            ctx.freezeAnalysis    = true; // matches existing "stop-time silence freeze" behavior intent
            ctx.absoluteSampleIndex = (int64_t) blockStartProjectSample;

            if (frameSamples > 0 && haveTimeInSamples)
            {
                ctx.hasFrameIndex60Hz = true;
                ctx.absoluteFrameIndex60Hz = (int64_t) floorDivInt64 (blockStartProjectSample,
                                                                     (juce::int64) frameSamples);
            }

            processorCore.process (ctx);

            // [BEGIN LS-IO-METERING-OUT-COMPUTE-STOPPATH]
            // Compute OUTPUT metering after ProcessorCore (even if analysis is frozen).
            const auto outMet = computeBlockPeakRmsMaxAcrossChannels (buffer, numChannels, numSamples);

            ioOutPeakDbCurrent.store (outMet.peakDb, std::memory_order_relaxed);
            ioOutRmsDbCurrent.store  (outMet.rmsDb,  std::memory_order_relaxed);

            {
                const float hold = ioOutPeakDbHold.load (std::memory_order_relaxed);
                const float dt = (currentSampleRate > 1.0e-6 ? (float) numSamples / (float) currentSampleRate : 0.0f);
                const float dec = ioPeakHoldDecayDbPerSec * dt;
                ioOutPeakDbHold.store (decayHoldDb (hold, outMet.peakDb, dec), std::memory_order_relaxed);
            }
            // [END LS-IO-METERING-OUT-COMPUTE-STOPPATH]

            return;
        }
    // [END LS-PROCESSORCORE-SKIP-ANALYSIS-STILL-PROCESS]

    // [BEGIN LS-PROCESSORCORE-PROCESSBLOCK]
    currentBlockStartProjectSample = blockStartProjectSample;
        currentBlockIsPlaying = blockIsPlaying;

        // RT: avoid per-block heap allocation (was juce::HeapBlock)
        const float* const* channelData = buffer.getArrayOfReadPointers();

        for (int i = 0; i < numSamples; ++i)
        {
            currentBlockSampleIndex = i;
            processSampleForLoudness (channelData, numChannels, i);
        }

        // Stage A: run ProcessorCore with an empty/default graph (no-op => no audible change).
        {
            levelscope::ProcessContext ctx { buffer,
                                            &midiMessages,
                                            currentSampleRate,
                                            numSamples,
                                            getBusesLayout().getMainInputChannelSet() };

            ctx.isRealtime         = ! isNonRealtime();
            ctx.isPlaying          = (blockIsPlaying == 1);
            ctx.isDiscontinuity    = discontinuity;
            ctx.freezeAnalysis     = false;
            ctx.absoluteSampleIndex = (int64_t) blockStartProjectSample;

            if (frameSamples > 0 && haveTimeInSamples)
            {
                ctx.hasFrameIndex60Hz = true;
                ctx.absoluteFrameIndex60Hz = (int64_t) floorDivInt64 (blockStartProjectSample,
                                                                     (juce::int64) frameSamples);
            }

            processorCore.process (ctx);
        }
        // [BEGIN LS-IO-METERING-OUT-COMPUTE]
        // Compute OUTPUT metering after ProcessorCore.
        const auto outMet = computeBlockPeakRmsMaxAcrossChannels (buffer, numChannels, numSamples);

        ioOutPeakDbCurrent.store (outMet.peakDb, std::memory_order_relaxed);
        ioOutRmsDbCurrent.store  (outMet.rmsDb,  std::memory_order_relaxed);

        {
            const float hold = ioOutPeakDbHold.load (std::memory_order_relaxed);
            const float dt = (currentSampleRate > 1.0e-6 ? (float) numSamples / (float) currentSampleRate : 0.0f);
            const float dec = ioPeakHoldDecayDbPerSec * dt;
            ioOutPeakDbHold.store (decayHoldDb (hold, outMet.peakDb, dec), std::memory_order_relaxed);
        }
        // [END LS-IO-METERING-OUT-COMPUTE]

        // Audio passthrough unchanged (empty graph => buffer is unchanged)
    // [END LS-PROCESSORCORE-PROCESSBLOCK]
}

//==============================================================================
// [CORE->UI] forwarding
//==============================================================================

int LevelScopeAudioProcessor::readLoudnessFromFifoEx (float* momentaryDest,
                                                     float* shortTermDest,
                                                     float* lraGateDest,
                                                     juce::int64* frameIndexDest,
                                                     int* isPlayingDest,
                                                     int maxNumToRead) noexcept
{
    return historyModel.readLoudnessFromFifo (momentaryDest,
                                             shortTermDest,
                                             lraGateDest,
                                             frameIndexDest,
                                             isPlayingDest,
                                             maxNumToRead);
}

int LevelScopeAudioProcessor::readLoudnessFromFifo (float* momentaryDest,
                                                   float* shortTermDest,
                                                   juce::int64* frameIndexDest,
                                                   int* isPlayingDest,
                                                   int maxNumToRead) noexcept
{
    // Backward-compatible wrapper (no gate)
    return historyModel.readLoudnessFromFifo (momentaryDest,
                                             shortTermDest,
                                             frameIndexDest,
                                             isPlayingDest,
                                             maxNumToRead);
}

//==============================================================================
// [STATE-PERSIST] forward to core (format unchanged)
//==============================================================================

// [BEGIN LS-C2-STATE-GETSET]
// [BEGIN LS-SETTINGS-PRESET-HELPERS-IMPL]
void LevelScopeAudioProcessor::buildApvsStateChunk (juce::MemoryBlock& destData)
{
    destData.setSize (0);

    juce::MemoryOutputStream os (destData, true);
    os.writeInt (1); // APVS schema version

    auto vt = apvts.copyState();
    vt.writeToStream (os);
}

// [BEGIN LVLR-MODG-BUILD-STATE-CHUNK]
void LevelScopeAudioProcessor::buildModgStateChunk (juce::MemoryBlock& destData)
{
    destData.setSize (0);

    juce::MemoryOutputStream os (destData, true);
    os.writeInt (1); // MODG schema version

    // Persist default graph order + bypass byte.
    const juce::String levelerId = juce::String (levelscope::LevelerModule().getModuleID());
    const juce::String mtdmId    = juce::String (levelscope::MultiThresholdDynamicsModule().getModuleID());

    os.writeInt (2); // num modules: Leveler -> MTDM

    os.writeString (levelerId);
    {
        const auto* enabled01 = apvts.getRawParameterValue (levelscope::lvlr::ParamIDs::enabled);
        const bool bypassed = (enabled01 != nullptr ? (enabled01->load() < 0.5f) : true);
        os.writeByte ((char) (bypassed ? 1 : 0));
    }

    os.writeString (mtdmId);
    {
        const auto* enabled01 = apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::enabled);
        const bool bypassed = (enabled01 != nullptr ? (enabled01->load() < 0.5f) : true);
        os.writeByte ((char) (bypassed ? 1 : 0));
    }
}
// [END LVLR-MODG-BUILD-STATE-CHUNK]

bool LevelScopeAudioProcessor::applyApvsStateChunk (const juce::MemoryBlock& apvsChunk)
{
    if (apvsChunk.getSize() == 0)
        return false;

    juce::MemoryInputStream in (apvsChunk.getData(), apvsChunk.getSize(), false);
    const int schema = in.readInt();

    if (schema != 1)
        return false;

    auto vt = juce::ValueTree::readFromStream (in);
    if (! vt.isValid())
        return false;

    apvts.replaceState (vt);
    // [BEGIN LVLR-HOSTGAIN-CAPTURE-DISARM-ON-LOAD]
    // Safety: capture arm must never come back armed after state/preset load.
    forceDisarmLevelerHostGainCapture_NonRT();
    // [END LVLR-HOSTGAIN-CAPTURE-DISARM-ON-LOAD]
    return true;
}

bool LevelScopeAudioProcessor::parseSettingsPresetBlob (const void* data, int sizeInBytes,
                                                        juce::MemoryBlock& apvsChunkOut,
                                                        juce::MemoryBlock& modgChunkOut) const
{
    apvsChunkOut.setSize (0);
    modgChunkOut.setSize (0);

    if (data == nullptr || sizeInBytes <= 8)
        return false;

    juce::MemoryInputStream in (data, (size_t) sizeInBytes, false);

    char magic[4] = {};
    if (in.read (magic, 4) != 4)
        return false;

    if (std::memcmp (magic, "LSSP", 4) != 0)
        return false;

    const int version = in.readInt();
    if (version != 1)
        return false;

    while (in.getNumBytesRemaining() >= 8)
    {
        char chunkId[4] = {};
        if (in.read (chunkId, 4) != 4)
            return false;

        const int chunkSize = in.readInt();
        if (chunkSize < 0 || in.getNumBytesRemaining() < (juce::int64) chunkSize)
            return false;

        juce::MemoryBlock chunk;
        chunk.setSize ((size_t) chunkSize);

        if (chunkSize > 0)
        {
            if (in.read (chunk.getData(), chunkSize) != chunkSize)
                return false;
        }

        if (std::memcmp (chunkId, "APVS", 4) == 0)
            apvsChunkOut = chunk;
        else if (std::memcmp (chunkId, "MODG", 4) == 0)
            modgChunkOut = chunk;
        // Unknown chunks ignored for forward compatibility
    }

    return (apvsChunkOut.getSize() > 0 || modgChunkOut.getSize() > 0);
}
// [END LS-SETTINGS-PRESET-HELPERS-IMPL]

// [BEGIN LS-SETTINGS-PRESET-GETSET]
void LevelScopeAudioProcessor::getSettingsPresetInformation (juce::MemoryBlock& destData)
{
    juce::MemoryBlock apvsChunk;
    juce::MemoryBlock modgChunk;

    buildApvsStateChunk (apvsChunk);
    buildModgStateChunk (modgChunk);

    destData.setSize (0);
    juce::MemoryOutputStream os (destData, true);

    // Simple versioned settings-preset container:
    // magic "LSSP" + version + chunk stream
    os.write ("LSSP", 4);
    os.writeInt (1);

    auto writeChunk = [&os] (const char* chunkId4, const juce::MemoryBlock& chunk)
    {
        os.write (chunkId4, 4);
        os.writeInt ((int) chunk.getSize());

        if (chunk.getSize() > 0)
            os.write (chunk.getData(), chunk.getSize());
    };

    writeChunk ("APVS", apvsChunk);
    writeChunk ("MODG", modgChunk);
}

bool LevelScopeAudioProcessor::setSettingsPresetInformation (const void* data, int sizeInBytes)
{
    juce::MemoryBlock apvsChunk;
    juce::MemoryBlock modgChunk;

    if (! parseSettingsPresetBlob (data, sizeInBytes, apvsChunk, modgChunk))
        return false;

    if (apvsChunk.getSize() > 0)
        (void) applyApvsStateChunk (apvsChunk);

    rebuildModuleGraphFromState (modgChunk.getSize() > 0 ? &modgChunk : nullptr);

    updateLatencyFromAPVTS_NonRT();
    return true;
}
// [END LS-SETTINGS-PRESET-GETSET]

void LevelScopeAudioProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    // Full snapshot: baseline history/timeline blob + APVS + MODG + UIST.
    juce::MemoryBlock apvsChunk;
    juce::MemoryBlock modgChunk;
    juce::MemoryBlock uistChunk;

    buildApvsStateChunk (apvsChunk);
    buildModgStateChunk (modgChunk);
    buildUistStateChunk (uistChunk);

    LevelScopeHistoryModel::ExtraStateChunksIn extra;
    extra.apvs = &apvsChunk;
    extra.modg = &modgChunk;
    extra.uist = &uistChunk;

    historyModel.saveState (destData, extra);
}

void LevelScopeAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    juce::MemoryBlock apvsChunk;
    juce::MemoryBlock modgChunk;
    juce::MemoryBlock uistChunk;

    LevelScopeHistoryModel::ExtraStateChunksOut extra;
    extra.apvs = &apvsChunk;
    extra.modg = &modgChunk;
    extra.uist = &uistChunk;

    historyModel.loadState (data, sizeInBytes, extra);

    if (apvsChunk.getSize() > 0)
        (void) applyApvsStateChunk (apvsChunk);

    rebuildModuleGraphFromState (modgChunk.getSize() > 0 ? &modgChunk : nullptr);

    if (uistChunk.getSize() > 0)
        (void) applyUistStateChunk (uistChunk);

    // [BEGIN LS-LATENCY-STATE-UPDATE]
    updateLatencyFromAPVTS_NonRT();
    // [END LS-LATENCY-STATE-UPDATE]
}
// [END LS-C2-STATE-GETSET]

//==============================================================================

juce::AudioProcessorEditor* LevelScopeAudioProcessor::createEditor()
{
    return new LevelScopeAudioProcessorEditor (*this);
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new LevelScopeAudioProcessor();
}