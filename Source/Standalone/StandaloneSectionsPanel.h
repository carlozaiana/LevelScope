#pragma once
#include "StandaloneSections.h"
#include "StandaloneSessionModel.h"
#include "MeasurementTimeline.h"
#include "StandaloneMarkerImportPanel.h"

namespace levelscope::standalone
{
class StandaloneSectionsPanel final : public juce::Component, private juce::Timer, private juce::ListBoxModel
{
public:
    explicit StandaloneSectionsPanel (StandaloneSessionModel&);
    ~StandaloneSectionsPanel() override;
    void synchronize() { timerCallback(); }
    void suspend();
    void resume();
    MeasurementTimeline& timelineView() {return timeline;}
    void refreshMeasuredView() {timeline.setAnalysis(measured?measured->analysis:nullptr);}
    MeasurementDomain inspectionDomain() const {return domainBox.getSelectedId()==2?MeasurementDomain::currentState:MeasurementDomain::source;}
    void setPlaybackState(bool canPlay,bool canStop,const juce::String& domain);
    bool keyPressed(const juce::KeyPress& key) override {return timeline.keyPressed(key);}
    std::function<void()> onPlay,onStop;
    void paint (juce::Graphics&) override;
    void resized() override;
    std::function<void()> onChanged;
    std::function<bool()> canImportMarkers;
private:
    StandaloneSessionModel& session;
    proposal::SectionDocument& document;
    DiagnosticsController diagnostics;
    SectionsController sections;
    std::shared_ptr<const MeasuredDiagnostics> measured;
    std::shared_ptr<const proposal::SectionAnalysis> displayed;
    juce::String mapKey, annotationKey, editNotice, appliedGainKey, intentViewKey;
    bool suspended = false, refreshing = false;
    uint64_t selectedId = 0;
    juce::TextButton play{"Play"},stop{"Stop"};
    juce::Label auditionDomain;
    juce::ComboBox domainBox, windowKind, windowScope;
    juce::TextButton splitButton { "Split at frame" }, moveButton { "Move end" }, mergeButton { "Merge next" }, renameButton { "Rename" };
    juce::TextButton rebuildButton { "Build sections" }, cancelButton { "Cancel sections" };
    juce::TextButton draftGainButton { "Draft section gains" }, setGainButton { "Set draft gain" },
        applyGainButton { "Apply section trims" }, discardGainButton { "Discard draft" };
    juce::TextButton importMarkers { "Import markers..." }, pasteMarkers { "Paste positions..." };
    std::unique_ptr<MarkerImportWindow> markerWindow;
    juce::TextEditor gainEditor;
    juce::Label gainStatus;
    struct GainView final : juce::Component
    {
        StandaloneSessionModel& session;
        explicit GainView(StandaloneSessionModel& s):session(s){}
        void paint(juce::Graphics&) override;
    } gainView { session };
    juce::TextEditor boundary, name, summary, inspection;
    juce::Label boundaryLabel, nameLabel, status;
    MeasurementTimeline timeline;
    juce::ListBox sectionList { "Programme sections", this };
    struct Windows final : juce::ListBoxModel
    {
        StandaloneSectionsPanel& panel;
        explicit Windows (StandaloneSectionsPanel& p) : panel (p) {}
        int getNumRows() override;
        void paintListBoxItem (int, juce::Graphics&, int, int, bool) override;
        void selectedRowsChanged (int) override;
    } windowsModel { *this };
    juce::ListBox windowList { "Section full-resolution windows", &windowsModel };
    void timerCallback() override;
    void refreshDetails();
    void performEdit (int operation);
    const proposal::ProgrammeSection* selectedSection() const;
    const proposal::SectionSummary* selectedSummary() const;
    const proposal::SectionWindowSummary* windowSummary() const;
    const proposal::WindowSubset* subset() const;
    const proposal::DiagnosticWindow* windowAt (int row) const;
    int getNumRows() override;
    void paintListBoxItem (int, juce::Graphics&, int, int, bool) override;
    void selectedRowsChanged (int) override;
};
}
