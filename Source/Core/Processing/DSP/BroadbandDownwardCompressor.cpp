#include "BroadbandDownwardCompressor.h"

namespace levelscope::dsp
{
void BroadbandDownwardCompressor::prepare (double sampleRate,
                                          int numChannels,
                                          const juce::AudioChannelSet& channelSet,
                                          int maxBlockSize)
{
    fs = (sampleRate > 0.0 ? sampleRate : 48000.0);
    preparedNumChannels = std::max (1, numChannels);
    preparedChannelSet = channelSet;
    preparedMaxBlockSize = std::max (0, maxBlockSize);

    // Prebuild default channel lists (exclude LFE)
    nonLfeChannels.clear();
    allChannels.clear();

    nonLfeChannels.reserve ((size_t) preparedNumChannels);
    allChannels.reserve ((size_t) preparedNumChannels);

    for (int ch = 0; ch < preparedNumChannels; ++ch)
    {
        allChannels.push_back (ch);

        const bool isLfe = (preparedChannelSet.getTypeOfChannel (ch) == juce::AudioChannelSet::LFE);
        if (! isLfe)
            nonLfeChannels.push_back (ch);
    }

    // fallback if layout is odd
    if (nonLfeChannels.empty())
        nonLfeChannels = allChannels;

    // default pointers
    detectPtr = &nonLfeChannels;
    applyPtr  = &nonLfeChannels;

    updateCoefficientsIfNeeded();
    reset();
}

void BroadbandDownwardCompressor::reset() noexcept
{
    envMS = 0.0f;
    gainZ = 1.0f;
}

void BroadbandDownwardCompressor::setParametersAudioThread (const Parameters& p) noexcept
{
    params = p;
    updateCoefficientsIfNeeded();
}

void BroadbandDownwardCompressor::setChannelListsAudioThread (const std::vector<int>* detectList,
                                                             const std::vector<int>* applyList) noexcept
{
    // Pointers must refer to prebuilt vectors owned by this object (nonLfeChannels/allChannels).
    // If null, fall back to default.
    detectPtr = (detectList != nullptr ? detectList : &nonLfeChannels);
    applyPtr  = (applyList  != nullptr ? applyList  : &nonLfeChannels);
}

// [BEGIN LS-BDC-LFE-POLICY-IMPL]
void BroadbandDownwardCompressor::setLfePolicyAudioThread (bool lfeInDetector, bool lfeInApply) noexcept
{
    detectPtr = (lfeInDetector ? &allChannels : &nonLfeChannels);
    applyPtr  = (lfeInApply    ? &allChannels : &nonLfeChannels);
}
// [END LS-BDC-LFE-POLICY-IMPL]

void BroadbandDownwardCompressor::updateCoefficientsIfNeeded() noexcept
{
    if (params.attackMs == lastAttackMs && params.releaseMs == lastReleaseMs)
        return;

    lastAttackMs  = params.attackMs;
    lastReleaseMs = params.releaseMs;

    const float attackS  = std::max (1.0e-4f, params.attackMs  * 0.001f);
    const float releaseS = std::max (1.0e-4f, params.releaseMs * 0.001f);

    const float sr = (float) std::max (1.0, fs);

    aA = std::exp (-1.0f / (attackS  * sr));
    aR = std::exp (-1.0f / (releaseS * sr));
}

float BroadbandDownwardCompressor::computeCompressionGain (float levelLufs) const noexcept
{
    const float ratio = std::max (1.0f, params.ratio);
    const float thr   = params.t2Lufs;
    const float knee  = std::max (0.0f, params.kneeDb);

    if (ratio <= 1.0001f)
        return 1.0f;

    const float x = levelLufs - thr;

    // Soft knee around threshold (width = kneeDb, symmetric)
    float overDb = 0.0f;

    if (knee <= 1.0e-3f)
    {
        overDb = std::max (0.0f, x);
    }
    else
    {
        const float half = 0.5f * knee;

        if (x <= -half)       overDb = 0.0f;
        else if (x >= half)   overDb = x;
        else
        {
            // quadratic ease-in within knee region
            const float t = (x + half) / knee; // 0..1
            overDb = t * t * half;             // 0..half (smooth)
        }
    }

    if (overDb <= 0.0f)
        return 1.0f;

    // Gain reduction dB = overDb * (1 - 1/ratio)
    const float reductionDb = overDb * (1.0f - 1.0f / ratio);
    const float g = dbToLin (-reductionDb);
    return juce::jlimit (0.0f, 1.0f, g);
}

void BroadbandDownwardCompressor::process (juce::AudioBuffer<float>& buffer) noexcept
{
    juce::ScopedNoDenormals noDenormals;

    if (! params.enabled)
        return;

    const int numSamples = buffer.getNumSamples();
    const int numChInBuf = buffer.getNumChannels();
    if (numSamples <= 0 || numChInBuf <= 0)
        return;

    float* const* chans = buffer.getArrayOfWritePointers();

    const auto& detectList = (detectPtr != nullptr ? *detectPtr : nonLfeChannels);
    const auto& applyList  = (applyPtr  != nullptr ? *applyPtr  : nonLfeChannels);

    const int chToProcess = std::min (preparedNumChannels > 0 ? preparedNumChannels : numChInBuf, numChInBuf);

    const int numDetect = (int) detectList.size();
    const int numApply  = (int) applyList.size();

    const float t2 = std::min (params.t2Lufs, params.t3Lufs);
    const float t3 = std::max (params.t2Lufs, params.t3Lufs);
    const float zoneWidth = std::max (1.0f, t3 - t2);

    const float makeupLin = dbToLin (params.makeupDb);

    // [BEGIN LS-BDC-METERING-INIT]
    float blockMinG  = 1.0f;
    float blockLastG = 1.0f;
    // [END LS-BDC-METERING-INIT]

    for (int i = 0; i < numSamples; ++i)
    {
        // Linked detector energy (default excludes LFE)
        double sumSq = 0.0;
        for (int di = 0; di < numDetect; ++di)
        {
            const int ch = detectList[(size_t) di];
            if (ch >= 0 && ch < chToProcess)
            {
                const float x = chans[ch][i];
                sumSq += (double) x * (double) x;
            }
        }

        const float e = (float) (sumSq / (double) std::max (1, numDetect));

        // detector smoothing (mean-square)
        const float a = (e > envMS ? aA : aR);
        envMS = a * envMS + (1.0f - a) * e;

        // LUFS-ish proxy
        const float L = (float) (-0.691 + 10.0 * std::log10 ((double) envMS + 1.0e-12));

        // T2–T3 engagement ramp: 0 below T2, 1 at/above T3
        float pos = (L - t2) / zoneWidth;
        const float zone01 = smoothstep01 (pos);

        // Full compression gain (<=1) based on threshold at T2
        const float gComp = computeCompressionGain (L);

        // Crossfade from unity to compression by zone amount
        const float gTarget = 1.0f + (gComp - 1.0f) * zone01;

        // Gain smoothing (use same A/R as detector for minimal control set)
        const float aG = (gTarget < gainZ ? aA : aR); // reducing gain uses "attack"
        gainZ = aG * gainZ + (1.0f - aG) * gTarget;

        const float gOut = gainZ * makeupLin;

        // [BEGIN LS-BDC-METERING-UPDATE]
            // Meter compressor gain only (exclude makeup)
            blockMinG  = std::min (blockMinG, gainZ);
            blockLastG = gainZ;
        // [END LS-BDC-METERING-UPDATE]

        // Apply to apply channels (default excludes LFE)
        for (int ai = 0; ai < numApply; ++ai)
        {
            const int ch = applyList[(size_t) ai];
            if (ch >= 0 && ch < chToProcess)
                chans[ch][i] *= gOut;
        }
    }
    // [BEGIN LS-BDC-METERING-STORE]
    lastBlockMinCompGain  = blockMinG;
    lastBlockLastCompGain = blockLastG;
    // [END LS-BDC-METERING-STORE]
}
} // namespace levelscope::dsp