#pragma once

#include <JuceHeader.h>
#include <vector>

class LevelScopeAudioProcessor;

//==============================================================================
// Displays momentary & short-term "loudness" as a scrolling, zoomable history.
// Using sample-based spacing on the X axis (frames from newest).
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

    void mouseDown (const juce::MouseEvent& event) override;

private:
    // juce::Timer
    void timerCallback() override;

    // Internal helpers
    void drainProcessorFifo();
    void pushLoudnessBatchToHistory (const float* momentaryValues,
                                     const float* shortTermValues,
                                     int numValues);

    enum class CurveKind
    {
        Momentary,
        ShortTerm
    };

    float getDbAtFrame (CurveKind kind,
                        juce::int64 frameIndex) const noexcept;

    float dbToY (float db, float height) const noexcept;

    void applyHorizontalZoom (float wheelDelta, float mouseX);
    void applyVerticalZoom (float wheelDelta);

    LevelScopeAudioProcessor& processor;

    const double visualFrameRate;        // loudness frames per second
    const double historyLengthSeconds;   // total history buffer length (seconds)

    const float minDb;                   // bottom of world dB range
    const float maxDb;                   // top of world dB range (0 dB)
    const float baseDbRange;             // |minDb| (absolute full dB span)

    const int historyCapacityFrames;     // number of loudness frames we store

    // Histories for momentary & short-term, in dB
    std::vector<float> historyDbMomentary;
    std::vector<float> historyDbShortTerm;

    // Total number of frames written so far (monotonic)
    juce::int64 historyFrameCount = 0;

    // X-axis offset in frames behind "now" (newest frame)
    double viewOffsetFrames  = 0.0;

    // X-axis zoom (frame spacing)
    // zoomX = pixels per loudness frame.
    double zoomX      = 1.0;
    double minZoomX   = 0.05;
    double maxZoomX   = 50.0;
    bool   hasCustomZoomX = false;   // to keep "initial ~10s" default until user changes zoom

    // Y-axis zoom (amplitude)
    double yZoom      = 1.0;   // >1 = zoom in (smaller dB span)
    double minYZoom   = 0.25;
    double maxYZoom   = 4.0;

    // Auto-fit mode: keep entire history visible horizontally
    bool  autoFitEnabled   = false;
    int   autoFitBarWidth  = 5;   // pixels on the left side

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VolumeHistoryComponent)
};