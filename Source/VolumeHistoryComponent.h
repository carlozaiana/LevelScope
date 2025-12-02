#pragma once

#include <JuceHeader.h>
#include <vector>

class LevelScopeAudioProcessor;

//==============================================================================
// Displays momentary & short-term loudness as a scrolling, zoomable history.
//
// Data model:
//   - RAW history at full frame rate.
//   - MID history at decimationFactorMid (min+max per group).
//   - OVERVIEW history at decimationFactorHigh (min+max per group).
//
// Drawing model:
//   - Always follows "now" (right edge = newest frame).
//   - Choose RAW/MID/OVERVIEW layer based on zoomX.
//   - Use fractional frame offsets so scroll is smooth even in MID/OVERVIEW.
//   - 60 Hz repaint in RAW, ~30 Hz in MID/OVERVIEW for performance.
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
    Frame getMidFrameAgo (int groupsAgo) const noexcept;
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

    // RAW history (full resolution)
    int                rawCapacityFrames = 0;   // number of RAW frames stored
    std::vector<Frame> rawHistory;
    int                rawWriteIndex   = 0;
    juce::int64        totalRawFrames  = 0;     // total RAW frames written since start

    // MID history (decimated, factor 16)
    static constexpr int decimationFactorMid  = 16;
    int                midCapacityFrames      = 0;
    std::vector<Frame> midHistory;
    int                midWriteIndex   = 0;
    juce::int64        totalMidFrames  = 0;
    Frame              currentMid;
    int                currentMidCount = 0;     // how many RAW frames in current MID group

    // OVERVIEW history (decimated, factor 64)
    static constexpr int decimationFactorHigh = 64;
    int                overviewCapacityFrames = 0;
    std::vector<Frame> overviewHistory;
    int                overviewWriteIndex  = 0;
    juce::int64        totalOverviewFrames = 0;
    Frame              currentOverview;
    int                currentOverviewCount = 0; // how many RAW frames in current OVERVIEW group

    // Zoom parameters
    double zoomX      = 5.0;   // pixels per RAW frame
    double minZoomX   = 0.0005;
    double maxZoomX   = 50.0;
    double zoomY      = 1.0;   // vertical zoom in dB
    double minZoomY   = 0.25;
    double maxZoomY   = 4.0;

    bool   hasCustomZoomX = false;

    // Repaint decimation for MID/OVERVIEW: keep track of timer ticks
    int repaintDecimator = 0;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VolumeHistoryComponent)
};