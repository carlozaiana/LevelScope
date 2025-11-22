#pragma once

#include <JuceHeader.h>
#include <vector>

class LevelScopeAudioProcessor;

//==============================================================================
// Displays incoming audio volume as a scrolling, zoomable history curve.
// Using sample-based spacing on the X axis.
//==============================================================================

class VolumeHistoryComponent : public juce::Component,
                               private juce::Timer
{
public:
    explicit VolumeHistoryComponent (LevelScopeAudioProcessor& processor);
    ~VolumeHistoryComponent() override;

    void paint (juce::Graphics& g) override;
    void resized() override;

    void mouseWheelMove (const juce::MouseEvent& event,
                         const juce::MouseWheelDetails& wheel) override;

private:
    // juce::Timer
    void timerCallback() override;

    // Internal helpers
    void drainProcessorFifo();
    void pushRmsBatchToHistory (const float* values, int numValues);

    float sampleIndexToDbNearest (juce::int64 sampleIndex) const noexcept;
    float dbToY (float db, float height) const noexcept;

    void applyHorizontalZoom (float wheelDelta, float mouseX);
    void applyVerticalZoom (float wheelDelta);

    LevelScopeAudioProcessor& processor;

    const double visualSampleRate;       // RMS frames per second
    const double historyLengthSeconds;   // total history buffer length in seconds

    const float minDb;                   // bottom of world dB range
    const float maxDb;                   // top of world dB range (0 dB)
    const float baseDbRange;             // |minDb| (absolute full dB span)

    const int historyCapacitySamples;    // number of RMS samples we store

    std::vector<float> historyDb;        // ring buffer of dB values

    // Total number of RMS samples written so far (monotonic)
    juce::int64 historySampleCount = 0;

    // X-axis offset in samples behind "now" (newest sample)
    // 0 = right edge shows the newest sample,
    // >0 = right edge shows an older sample.
    double viewOffsetSamples  = 0.0;

    // X-axis zoom (sample spacing)
    // zoomX = pixels per RMS sample.
    double zoomX      = 1.0;
    double minZoomX   = 0.01;
    double maxZoomX   = 50.0;
    bool   hasCustomZoomX = false;   // to keep "initial ~10s" default until user changes zoom

    // Y-axis zoom (amplitude)
    double yZoom      = 1.0;   // >1 = zoom in (smaller dB span)
    double minYZoom   = 0.25;
    double maxYZoom   = 4.0;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VolumeHistoryComponent)
};