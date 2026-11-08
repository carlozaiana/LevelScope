#pragma once
#include "Core/Proposal/AnalysisEngine.h"
#include "Core/Waveform/WaveformOverview.h"
#include "StandaloneLoudnessPreview.h"
#include "TimelineAppearance.h"
#include "SectionEditPrompt.h"
#include <juce_gui_basics/juce_gui_basics.h>

namespace levelscope::standalone
{
// Reusable measured-data inspector. Deliberately has no session, player, device,
// processing or profile reference. All navigation is in audio-source frames.
class TimelineGainLane;
class StandaloneSessionModel;
class MeasurementTimeline final : public juce::Component, public juce::TooltipClient, private juce::ScrollBar::Listener
{
public:
    MeasurementTimeline();
    void attachGainLane(TimelineGainLane&);
    void bindSectionEditing(StandaloneSessionModel&,std::function<void(uint64_t)>);
    void refreshSectionEditing();
    void changeSectionEnd(uint64_t,int64_t);
    SectionEditPrompt* editPrompt(){return sectionEdit.get();}
    juce::String peakHoverText(juce::Point<float>) const;
    juce::String getTooltip() override {return peakHoverText(getMouseXYRelative().toFloat());}
    void mouseMove(const juce::MouseEvent&) override;
    void focusForTransport();
    bool keyPressed(const juce::KeyPress&) override;
    std::function<void()> onTogglePlayback;
    std::function<void()> onToggleInspectionAudition;
    void setProgramme(int64_t frames,int sampleRate,const juce::String& identity);
    void setAnalysis (std::shared_ptr<const proposal::ProgramAnalysis>);
    void setWaveform(std::shared_ptr<const waveform::Overview>,bool accepted,const juce::StringArray& channels);
    void setPreview(std::shared_ptr<const LoudnessPreview>);
    void setReferenceMeasurement(std::shared_ptr<const measurement::Result>,bool inspectingCurrent);
    juce::String rangePeakText(int64_t begin,int64_t end) const;
    TimelineBlend& visibilityBlend(){return blend;}
    void setPeakTarget(std::optional<double> dbtp,bool measured);
    juce::ToggleButton& peakVisibility(){return overTargetPeaks;}
    bool sampleExceedsTarget(double value,int64_t frame) const;
    const std::vector<proposal::PlotBucket>& shortTermGatePlot() const {return gatedShortTerm;}
    const std::vector<bool>& changedMomentaryColumns() const {return changedMomentary;}
    void setTransport(int64_t frame,bool seekAllowed,const juce::String& reason,const juce::String& identity);
    std::function<juce::String(int64_t,const juce::String&)> onRequestSeek;
    std::pair<double,double> viewRange() const {return {viewBegin,viewEnd};}
    juce::Rectangle<int> rulerBounds() const {return ruler;}
    juce::Rectangle<int> graphBounds() const {return graph;}
    static juce::String timecode(int64_t frame,int sampleRate);
    const std::vector<waveform::Column>& displayedWaveform() const {return waveformColumns;}
    void setViewRange (double beginFrame, double endFrame);
    void highlight (double beginFrame, double endFrame, bool reveal = true);
    void clearHighlight() { selection.reset();inspectedRange.reset();fitSelection.setEnabled(false);repaint(); }
    void setSectionBoundaries (std::vector<int64_t> frames) { sectionBoundaries = std::move (frames); repaint(); }
    bool inspectSeconds (double seconds);
    juce::String inspectionText() const { return inspection; }
    juce::String navigationStatus() const {return navigationFeedback.getText();}
    const proposal::TimelinePlot& displayedPlot() const { return plot; }
    void paint (juce::Graphics&) override;
    void resized() override;
    void mouseDown (const juce::MouseEvent&) override;
    void mouseDrag (const juce::MouseEvent&) override;
    void mouseUp (const juce::MouseEvent&) override;
    void mouseWheelMove (const juce::MouseEvent&, const juce::MouseWheelDetails&) override;
    std::function<void(const juce::String&)> onInspect;
private:
    StandaloneSessionModel* sectionSession=nullptr;
    std::unique_ptr<SectionEditPrompt> sectionEdit;
    juce::TooltipWindow hoverTips{this,350};
    std::function<void(uint64_t)> sectionEdited;
    juce::String sectionInput,sectionRuntime;
    uint64_t sectionRevision=0, audioRevision=0;
    bool sectionEditable=false, sectionGesture=false;
    std::optional<size_t> boundaryHandle(float x) const;
    std::optional<size_t> movingBoundary;
    int64_t dragFrame=0;
    void finishSectionGesture(int64_t end);
    TimelineGainLane* gainLane=nullptr; // borrowed from the main shell
    std::shared_ptr<const proposal::ProgramAnalysis> analysis;
    std::shared_ptr<const LoudnessPreview> preview;
    int64_t previewCoverage=0;
    int64_t knownFrames=0;
    int knownRate=1;
    juce::String programmeIdentity;
    proposal::TimelinePlot plot;
    TimelineBlend blend;
    std::shared_ptr<const measurement::Result> originalMeasurement;
    std::optional<double> originalIntegratedGateEnergy;
    bool currentInspection=false;
    std::vector<proposal::PlotBucket> gatedShortTerm;
    std::vector<bool> changedMomentary;
    std::shared_ptr<const waveform::Overview> waveformData;
    std::vector<waveform::Column> waveformColumns;
    juce::StringArray channelNames;
    uint64_t waveformPublication=0;
    bool waveformAccepted=false,waveformEnabled=false;
    int inspectedChannel=0;
    int64_t programmeFrames() const {return waveformData?waveformData->spec.frames:preview?preview->frames:analysis?analysis->frames:knownFrames;}
    int programmeRate() const {return waveformData?waveformData->spec.sampleRate:preview?preview->sampleRate:analysis?analysis->sampleRate:knownRate;}
    void paintWaveform(juce::Graphics&);
    double viewBegin = 0, viewEnd = 1;
    std::optional<std::pair<double, double>> selection;
    std::optional<std::pair<double,double>> inspectedRange;
    std::optional<int64_t> inspectedFrame;
    juce::String inspectedEntry;
    std::vector<int64_t> sectionBoundaries; // read-only annotations, never seeks or edits
    juce::Rectangle<int> graph,ruler;
    juce::TextButton zoomIn { "Zoom +" }, zoomOut { "Zoom -" }, fit { "Whole programme" }, fitSelection {"Fit selection"}, inspect { "Inspect seconds" };
    juce::TextButton seek {"Seek"};
    juce::ToggleButton follow {"Follow"};
    juce::ToggleButton overTargetPeaks {"Peaks over target"};
    std::optional<double> peakTarget;
    double peakCeiling=1;
    bool peakMeasurementCurrent=false;
    bool peakOverlayActive() const {return overTargetPeaks.getToggleState()&&peakTarget&&peakMeasurementCurrent&&!preview;}
    juce::Label navigationFeedback;
    int64_t playhead=0;
    bool canSeek=false,dragging=false,gesture=false,gestureRuler=false,gestureGain=false;
    juce::String seekReason,transportIdentity,gestureIdentity;
    uint64_t navigationEpoch=0,gestureEpoch=0;
    juce::Point<float> downPoint;
    double dragBegin=0,dragEnd=1;
    juce::TextEditor seconds;
    juce::ScrollBar scroll { false };
    juce::String inspection;
    bool updatingScroll = false;
    void rebuild();
    bool inspectFrame(double);
    int64_t frameAt(float x) const;
    int64_t enteredFrame() const;
    void requestSeek(int64_t,const juce::String&);
    void navigationMessage(const juce::String&);
    void zoom (double factor);
    void zoomAt(double factor,float x);
    void scrollBarMoved (juce::ScrollBar*, double) override;
};
juce::String diagnosticWindowText (const proposal::DiagnosticWindow&, int sampleRate, int64_t programmeFrames, bool ungated = false);
juce::String diagnosticMetricText (const measurement::Metric&, const char* unit);
juce::String diagnosticLraInterpretation (const proposal::ProgramAnalysis&);
}
