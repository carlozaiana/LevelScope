#include "PluginProcessor.h"
#include "PluginEditor.h"

// [BEGIN MTDM-PLUGIN-INCLUDE]
#include "Core/Processing/Modules/MultiThresholdDynamicsModule.h"
// [END MTDM-PLUGIN-INCLUDE]

// [BEGIN MTDM-PARAMIDS-INCLUDE]
#include "Core/Processing/Modules/MultiThresholdDynamicsParamIDs.h"
// [END MTDM-PARAMIDS-INCLUDE]

#include <cmath>
#include <algorithm>
#include <type_traits>

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
}

// [BEGIN MTDM-APVTS-PARAM-LAYOUT]
juce::AudioProcessorValueTreeState::ParameterLayout LevelScopeAudioProcessor::createParameterLayout()
{
    juce::AudioProcessorValueTreeState::ParameterLayout layout;

    using namespace levelscope::mtdm;

    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::enabled, 1 },
        "MTDM Enabled",
        (Defaults::enabled01 >= 0.5f)));

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
        Defaults::sucFftSizeChoice));

    layout.add (std::make_unique<juce::AudioParameterChoice> (
        juce::ParameterID { ParamIDs::sucBandsPerOctChoice, 1 },
        "SUC Bands/Oct",
        juce::StringArray { "1", "2", "3", "6" },
        Defaults::sucBandsPerOctChoice));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::sucMinFreqHz, 1 },
        "SUC Min Freq (Hz)",
        juce::NormalisableRange<float> (Ranges::sucMinFreqMinHz, Ranges::sucMinFreqMaxHz, 1.0f),
        Defaults::sucMinFreqHz));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::sucMaxFreqHz, 1 },
        "SUC Max Freq (Hz)",
        juce::NormalisableRange<float> (Ranges::sucMaxFreqMinHz, Ranges::sucMaxFreqMaxHz, 1.0f),
        Defaults::sucMaxFreqHz));
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
        Defaults::upwardModeChoice));
    // [END MTDM-APVTS-PARAM-LAYOUT-UPWARD-MODE]

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

    // [BEGIN MTDM-APVTS-PARAM-LAYOUT-LIMITER]
    layout.add (std::make_unique<juce::AudioParameterBool> (
        juce::ParameterID { ParamIDs::limEnabled01, 1 },
        "Limiter Enabled",
        (Defaults::limEnabled01 >= 0.5f)));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::limCeilingDb, 1 },
        "Limiter Ceiling (dBFS)",
        juce::NormalisableRange<float> (Ranges::limCeilingMinDb, Ranges::limCeilingMaxDb, 0.1f),
        Defaults::limCeilingDb));

    layout.add (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ParamIDs::limLookaheadMs, 1 },
        "Limiter Lookahead (ms)",
        juce::NormalisableRange<float> (Ranges::limLookaheadMinMs, Ranges::limLookaheadMaxMs, 0.1f),
        Defaults::limLookaheadMs));

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
        Defaults::limOversamplingChoice));
    // [END MTDM-APVTS-PARAM-LAYOUT-LIMITER-TP-ATTACK-DRIVE]
    // [END MTDM-APVTS-PARAM-LAYOUT-LIMITER]

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
    const juce::String mtdmId = juce::String (levelscope::MultiThresholdDynamicsModule().getModuleID());

    std::vector<juce::String> orderedModuleIds;
    orderedModuleIds.reserve (4);

    // Default graph if no MODG chunk (or invalid)
    orderedModuleIds.push_back (mtdmId);

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
                    orderedModuleIds.push_back (mtdmId);
            }
        }
    }

    auto graph = std::make_shared<levelscope::ModuleGraph>();
    graph->revision = 3;

    for (const auto& id : orderedModuleIds)
    {
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
                                  apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::limOversamplingChoice));

            // [BEGIN LS-MTDM-UI-HANDLE-STORE]
            std::atomic_store_explicit (&mtdmForUI, mtdm, std::memory_order_release);
            // [END LS-MTDM-UI-HANDLE-STORE]                  

            graph->modules.push_back (mtdm);
        }
        // [END MTDM-BINDPARAMS-CALL-FULL-D1C]
        // Unknown modules are ignored (forward-compat)
    }

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

// [BEGIN LS-LATENCY-HELPER-IMPL]
void LevelScopeAudioProcessor::updateLatencyFromAPVTS_NonRT()
{
    int totalLatency = 0;

    const auto* enabled01 = apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::enabled);
    const bool mtdmEnabled = (enabled01 != nullptr && enabled01->load (std::memory_order_relaxed) >= 0.5f);

    if (! mtdmEnabled)
    {
        setLatencySamples (0);
        return;
    }

    // --- Upward latency
    int upwardLatency = 0;

    const auto* upwardModeP = apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::upwardModeChoice);
    const int upwardMode = (int) std::lround (upwardModeP != nullptr
                                              ? upwardModeP->load (std::memory_order_relaxed)
                                              : (float) levelscope::mtdm::Defaults::upwardModeChoice);

    if (upwardMode == 0) // Spectral
    {
        const auto* choiceP = apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::sucFftSizeChoice);
        const int choice = (int) std::lround (choiceP != nullptr ? choiceP->load (std::memory_order_relaxed)
                                                                 : (float) levelscope::mtdm::Defaults::sucFftSizeChoice);

        const int fftSizes[] = { 1024, 2048, 4096, 8192 };
        upwardLatency = fftSizes[juce::jlimit (0, 3, choice)];
    }

    // [BEGIN LS-LATENCY-LIMITER-INCLUDES-FIR]
    int limiterLatency = 0;

    const auto* limEnabledP = apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::limEnabled01);
    const bool limEnabled = (limEnabledP != nullptr && limEnabledP->load (std::memory_order_relaxed) >= 0.5f);

    if (limEnabled)
    {
        const auto* lookMsP = apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::limLookaheadMs);
        const float lookMs = (lookMsP != nullptr ? lookMsP->load (std::memory_order_relaxed)
                                                 : levelscope::mtdm::Defaults::limLookaheadMs);

        const double sr = (currentSampleRate > 0.0 ? currentSampleRate : 48000.0);
        const double ms = juce::jlimit (0.0, 50.0, (double) lookMs);

        const int lookSamples = (int) std::lround (sr * ms * 0.001);

        // FIR oversampling detector delay (base-rate samples)
        int detDelay = 0;
        const auto* osP = apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::limOversamplingChoice);
        const int osChoice = (int) std::lround (osP != nullptr ? osP->load (std::memory_order_relaxed)
                                                               : (float) levelscope::mtdm::Defaults::limOversamplingChoice);

        const int clamped = juce::jlimit (0, 2, osChoice);
        if (clamped > 0)
        {
            const int stages = (clamped == 1 ? 1 : 2); // 2x=1 stage, 4x=2 stages
            juce::dsp::Oversampling<float> os (1, stages,
                juce::dsp::Oversampling<float>::filterHalfBandFIREquiripple, true);
            os.initProcessing (1);
            detDelay = os.getLatencyInSamples();
        }

        limiterLatency = lookSamples + detDelay;
    }
    // [END LS-LATENCY-LIMITER-INCLUDES-FIR]
}
// [END LS-LATENCY-HELPER-IMPL]

// [BEGIN LS-LATENCY-LISTENER-IMPL]

LevelScopeAudioProcessor::~LevelScopeAudioProcessor()
{
    stopTimer();
    unregisterLatencyParamListeners();
}

void LevelScopeAudioProcessor::registerLatencyParamListeners()
{
    using namespace levelscope::mtdm::ParamIDs;

    // Only params that can change total latency.
    apvts.addParameterListener (enabled, this);
    apvts.addParameterListener (upwardModeChoice, this);
    apvts.addParameterListener (sucFftSizeChoice, this);

    apvts.addParameterListener (limEnabled01, this);
    apvts.addParameterListener (limLookaheadMs, this);
    apvts.addParameterListener (limOversamplingChoice, this);

    // 10 Hz polling on message thread to apply latency updates non-RT.
    startTimerHz (10);
}

void LevelScopeAudioProcessor::unregisterLatencyParamListeners()
{
    using namespace levelscope::mtdm::ParamIDs;

    apvts.removeParameterListener (enabled, this);
    apvts.removeParameterListener (upwardModeChoice, this);
    apvts.removeParameterListener (sucFftSizeChoice, this);

    apvts.removeParameterListener (limEnabled01, this);
    apvts.removeParameterListener (limLookaheadMs, this);
    apvts.removeParameterListener (limOversamplingChoice, this);
}

void LevelScopeAudioProcessor::parameterChanged (const juce::String& parameterID, float newValue)
{
    juce::ignoreUnused (parameterID, newValue);

    // RT-safe: just mark dirty. Do NOT call setLatencySamples() here.
    latencyDirty.store (true, std::memory_order_release);
}

void LevelScopeAudioProcessor::timerCallback()
{
    if (! latencyDirty.exchange (false, std::memory_order_acq_rel))
        return;

    updateLatencyFromAPVTS_NonRT();

    // Helps some hosts refresh; safe on message thread.
    updateHostDisplay();
}
// [END LS-LATENCY-LISTENER-IMPL]

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

//==============================================================================

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

    // RT-safe reset; Stage A graph is empty => no audible change
    processorCore.reset();
}
// [END LS-PROCESSORCORE-RESETLOUDNESSSTATE]

//==============================================================================

bool LevelScopeAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
    const auto& mainIn  = layouts.getMainInputChannelSet();
    const auto& mainOut = layouts.getMainOutputChannelSet();

    if (mainIn.isDisabled() || mainOut.isDisabled())
        return false;

    // Prototype: mono or stereo only
    if (mainIn.size() > 2 || mainOut.size() > 2)
        return false;

    if (mainIn != mainOut)
        return false;

    return true;
}

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
void LevelScopeAudioProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    // Baseline chunks are written exactly as before, then we append APVS + MODG.
    juce::MemoryBlock apvsChunk;
    {
        juce::MemoryOutputStream os (apvsChunk, true);
        os.writeInt (1); // APVS schema version

        // Store the APVTS ValueTree in binary form (non-audio-thread only).
        auto vt = apvts.copyState();
        vt.writeToStream (os);
    }

    juce::MemoryBlock modgChunk;
    {
        juce::MemoryOutputStream os (modgChunk, true);
        os.writeInt (1); // MODG schema version

        // For Stage C2, we persist: module ID list (order) + bypass flag.
        // Current graph is effectively fixed (MTDM only), but this chunk establishes the additive format.
        const juce::String mtdmId = juce::String (levelscope::MultiThresholdDynamicsModule().getModuleID());

        os.writeInt (1); // num modules
        os.writeString (mtdmId);

        // "bypass" persisted here as a graph-level concept.
        // For now (Stage C2), we mirror the MTDM Enabled param for persistence.
        const auto* enabled01 = apvts.getRawParameterValue (levelscope::mtdm::ParamIDs::enabled);
        const bool bypassed = (enabled01 != nullptr ? (enabled01->load() < 0.5f) : true);
        os.writeByte ((char) (bypassed ? 1 : 0));
    }

    historyModel.saveState (destData, &apvsChunk, &modgChunk);
}

void LevelScopeAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    juce::MemoryBlock apvsChunk;
    juce::MemoryBlock modgChunk;

    historyModel.loadState (data, sizeInBytes, &apvsChunk, &modgChunk);

    // Restore APVTS if present (old sessions won't have APVS).
    if (apvsChunk.getSize() > 0)
    {
        juce::MemoryInputStream in (apvsChunk.getData(), apvsChunk.getSize(), false);
        const int schema = in.readInt();

        if (schema == 1)
        {
            auto vt = juce::ValueTree::readFromStream (in);
            if (vt.isValid())
                apvts.replaceState (vt);
        }
    }

    // Rebuild module graph safely (old sessions won't have MODG).
    rebuildModuleGraphFromState (modgChunk.getSize() > 0 ? &modgChunk : nullptr);

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