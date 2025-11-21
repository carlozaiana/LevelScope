#pragma once

#include <JuceHeader.h>
#include <vector>

class LevelScopeAudioProcessor;

//==============================================================================
// Displays incoming audio volume as a scrolling, zoomable history curve.
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

    float sampleIndexToDb (double sampleIndex) const noexcept;
    float dbToY (float db, float height) const noexcept;

    void applyHorizontalZoom (float wheelDelta, float mouseX);
    void applyVerticalZoom (float wheelDelta);

    LevelScopeAudioProcessor& processor;

    const double visualSampleRate;       // RMS frames per second
    const double historyLengthSeconds;   // total history buffer length in seconds

    const float minDb;                   // bottom of world dB range
    const float maxDb;                   // top of world dB range (0 dB)
    const float baseDbRange;             // maxDb - minDb (absolute dB span)

    const int historyCapacitySamples;    // number of RMS samples we store

    std::vector<float> historyDb;        // ring buffer of dB values

    juce::int64 historySampleCount = 0;  // total RMS samples written so far
    double      viewOffsetSamples  = 0;  // how many RMS samples behind "now" the right edge is

    // X-axis zoom (time)
    double visibleDurationSeconds;       // currently visible time window (seconds)
    double minVisibleSeconds;
    double maxVisibleSeconds;

    // Y-axis zoom (amplitude)
    double yZoom      = 1.0;             // >1 = zoom in (smaller dB span)
    double minYZoom   = 0.25;
    double maxYZoom   = 4.0;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VolumeHistoryComponent)
};