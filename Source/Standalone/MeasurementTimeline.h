#pragma once
#include "Core/Proposal/AnalysisEngine.h"
#include <juce_gui_basics/juce_gui_basics.h>

namespace levelscope::standalone
{
// Reusable measured-data inspector. Deliberately has no session, player, device,
// processing or profile reference. All navigation is in audio-source frames.
class MeasurementTimeline final : public juce::Component, private juce::ScrollBar::Listener
{
public:
    MeasurementTimeline();
    void setAnalysis (std::shared_ptr<const proposal::ProgramAnalysis>);
    void setViewRange (double beginFrame, double endFrame);
    void highlight (double beginFrame, double endFrame, bool reveal = true);
    void clearHighlight() { selection.reset(); repaint(); }
    bool inspectSeconds (double seconds);
    juce::String inspectionText() const { return inspection; }
    const proposal::TimelinePlot& displayedPlot() const { return plot; }
    void paint (juce::Graphics&) override;
    void resized() override;
    void mouseDown (const juce::MouseEvent&) override;
    void mouseWheelMove (const juce::MouseEvent&, const juce::MouseWheelDetails&) override;
    std::function<void(const juce::String&)> onInspect;
private:
    std::shared_ptr<const proposal::ProgramAnalysis> analysis;
    proposal::TimelinePlot plot;
    double viewBegin = 0, viewEnd = 1;
    std::optional<std::pair<double, double>> selection;
    juce::Rectangle<int> graph;
    juce::TextButton zoomIn { "Zoom +" }, zoomOut { "Zoom -" }, fit { "Whole programme" }, inspect { "Inspect seconds" };
    juce::TextEditor seconds;
    juce::ScrollBar scroll { false };
    juce::String inspection;
    bool updatingScroll = false;
    void rebuild();
    void zoom (double factor);
    void scrollBarMoved (juce::ScrollBar*, double) override;
};
juce::String diagnosticWindowText (const proposal::DiagnosticWindow&, int sampleRate, int64_t programmeFrames, bool ungated = false);
juce::String diagnosticMetricText (const measurement::Metric&, const char* unit);
juce::String diagnosticLraInterpretation (const proposal::ProgramAnalysis&);
}
