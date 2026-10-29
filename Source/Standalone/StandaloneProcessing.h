#pragma once
#include "StandaloneServices.h"
#include <cmath>

namespace levelscope::standalone
{
// The qualified static-gain-v1 operation. No channel-role weighting, clipping,
// monitor controls or plugin processing belongs here.
class StaticGainProcessor
{
public:
    static bool supports (const ProcessingPlan& p) noexcept
    {
        return p.schemaVersion == 1 && p.engineVersion == "static-gain-v1"
            && p.graphJson.isEmpty() && !p.sectionGains && std::isfinite (p.staticGainDb)
            && p.staticGainDb >= -60 && p.staticGainDb <= 24;
    }
    explicit StaticGainProcessor (double db) : coefficient (std::pow (10.0, db / 20)) {}
    bool process (float* const* audio, int channels, int frames) const noexcept
    {
        for (int channel = 0; channel < channels; ++channel)
            for (int i = 0; i < frames; ++i)
            {
                const float x = static_cast<float> (audio[channel][i] * coefficient);
                if (! std::isfinite (x)) return false;
                audio[channel][i] = x;
            }
        return true;
    }
private:
    double coefficient;
};

juce::String sectionInputIdentity(const ServiceContext&);

// Shared worker-side processing. The callback, routing and monitor stage remain
// unchanged. Frame coordinates come from the reader, including after a seek.
class ProgrammeGainProcessor
{
public:
    static juce::String validationError(const ServiceContext& c)
    {
        const auto& p=c.processing;
        if(StaticGainProcessor::supports(p)) return {};
        if(p.schemaVersion!=2 || p.engineVersion!=processing::SectionGainPlan::version || !p.sectionGains
            || p.graphJson.isNotEmpty() || !std::isfinite(p.staticGainDb) || p.staticGainDb < -60 || p.staticGainDb > 24)
            return "Unsupported Current processing plan";
        const auto& plan=*p.sectionGains;
        const auto error=plan.validationError(); if(error.isNotEmpty()) return error;
        // Product boundary, not a six-channel assumption in the sample evaluator.
        juce::AudioChannelSet roles;
        for(int n=0;n<c.fileLayout.count;++n)
        {
            auto role=c.fileLayout.roles[n];
            if(role==juce::AudioChannelSet::leftSurroundSide) role=juce::AudioChannelSet::leftSurround;
            if(role==juce::AudioChannelSet::rightSurroundSide) role=juce::AudioChannelSet::rightSurround;
            roles.addChannel(role);
        }
        const bool mono=c.fileLayout.count==1 && (c.fileLayout.roles[0]==juce::AudioChannelSet::left
            || c.fileLayout.roles[0]==juce::AudioChannelSet::right || c.fileLayout.roles[0]==juce::AudioChannelSet::centre);
        if(!c.fileLayout.confirmed || !c.fileLayout.isValid() || (!mono
            && ((c.fileLayout.count!=2 && c.fileLayout.count!=6)
                || (roles!=juce::AudioChannelSet::stereo() && roles!=juce::AudioChannelSet::create5point1()))))
            return "Section processing supports confirmed mono (Left/Right/Centre), stereo and conventional 5.1 only";
        if(plan.inputIdentity!=sectionInputIdentity(c) || plan.frames!=c.media.frames || plan.sampleRate!=c.media.sampleRate)
            return "Section plan does not describe this exact programme input";
        if(plan.isUnity()) return "Unity section plans must use the legacy static-gain identity";
        return {};
    }
    explicit ProgrammeGainProcessor(const ProcessingPlan& p) : global(p.staticGainDb)
    { if(p.sectionGains) sections=p.sectionGains->canonical(); }
    bool process(float* const* audio,int channels,int frames,int64_t firstFrame) const noexcept
    {
        if(channels<1 || channels>maxChannels || frames<0 || firstFrame<0) return false;
        if(sections)
        {
            if(firstFrame>sections->frames || frames>sections->frames-firstFrame) return false;
            for(int n=0;n<frames;++n)
            {
                const auto db=sections->gainAt(firstFrame+n);
                if(db==0) continue;
                const double coefficient=std::pow(10.0,db/20.0);
                for(int ch=0;ch<channels;++ch)
                {
                    const float x=static_cast<float>(audio[ch][n]*coefficient);
                    if(!std::isfinite(x)) return false;
                    audio[ch][n]=x;
                }
            }
        }
        return global.process(audio,channels,frames);
    }
private:
    StaticGainProcessor global;
    std::optional<processing::SectionGainPlan> sections;
};

// Injectable test sink; caller owns storage/lifetime. RT implementations must be
// bounded, allocation-free and noexcept. No observer is installed in the product.
struct ProcessedAudioObserver
{
    virtual ~ProcessedAudioObserver() = default;
    virtual void observe (const float* const*, int channels, int frames, int64_t sourceFrame) noexcept = 0;
};
}
