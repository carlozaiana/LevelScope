#pragma once

#include <JuceHeader.h>
#include <vector>

class LevelScopeAudioProcessor;

//==============================================================================
// Displays momentary & short-term loudness as a scrolling, zoomable history.
// Uses a RAW history and an OVERVIEW history (decimated) with fixed grouping.
// OVERVIEW stores per-group min+max to preserve both peaks and valleys.
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
    void pushLoudnessBatchToHistory (const float* momentaryValues,
                                     const float* shortTermValues,
                                     int numValues);

    float dbToY (float db, float height) const noexcept;

    void applyHorizontalZoom (float wheelDelta);
    void applyVerticalZoom (float wheelDelta);

    // Data access helpers
    struct Frame
    {
        float momentaryMinDb = -90.0f;
        float momentaryMaxDb = -90.0f;
        float shortTermMinDb = -90.0f;
        float shortTermMaxDb = -90.0f;
    };

    Frame getRawFrameAgo (int framesAgo) const noexcept;
    Frame getOverviewFrameAgo (int groupsAgo) const noexcept;

    LevelScopeAudioProcessor& processor;

    // Loudness frame rate (frames per second)
    const double visualFrameRate;

    // Total history length in seconds (RAW buffer capacity)
    const double historyLengthSeconds;

    // dB range
    const float minDb;
    const float maxDb;
    const float baseDbRange;

    // RAW history
    int                rawCapacityFrames = 0;   // number of RAW frames stored
    std::vector<Frame> rawHistory;
    int                rawWriteIndex   = 0;
    juce::int64        totalRawFrames  = 0;     // total RAW frames written since start

    // OVERVIEW history (decimated)
    static constexpr int decimationFactor = 64;   // RAW->OVERVIEW grouping size

    int                overviewCapacityFrames = 0; // rawCapacityFrames / decimationFactor
    std::vector<Frame> overviewHistory;
    int                overviewWriteIndex  = 0;
    juce::int64        totalOverviewFrames = 0;   // total OVERVIEW frames written

    // Accumulator for current overview group
    Frame currentOverview;
    int   currentOverviewCount = 0;

    // Zoom parameters
    double zoomX      = 5.0;   // pixels per RAW frame
    double minZoomX   = 0.0005;
    double maxZoomX   = 50.0;
    double zoomY      = 1.0;   // vertical zoom in dB
    double minZoomY   = 0.25;
    double maxZoomY   = 4.0;

    bool   hasCustomZoomX = false;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VolumeHistoryComponent)
};