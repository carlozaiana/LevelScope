#pragma once
#include "StandaloneSections.h"
#include "StandaloneSessionModel.h"
#include "MeasurementTimeline.h"

namespace levelscope::standalone
{
class StandaloneSectionsPanel final : public juce::Component, private juce::Timer, private juce::ListBoxModel
{
public:
    explicit StandaloneSectionsPanel (const StandaloneSessionModel&);
    ~StandaloneSectionsPanel() override;
    void synchronize() { timerCallback(); }
    void suspend();
    void resume();
    void paint (juce::Graphics&) override;
    void resized() override;
private:
    const StandaloneSessionModel& session;
    proposal::SectionDocument document;
    DiagnosticsController diagnostics;
    SectionsController sections;
    std::shared_ptr<const MeasuredDiagnostics> measured;
    std::shared_ptr<const proposal::SectionAnalysis> displayed;
    juce::String mapKey, annotationKey, editNotice;
    bool suspended = false, refreshing = false;
    uint64_t selectedId = 0;
    juce::ComboBox domainBox, windowKind, windowScope;
    juce::TextButton splitButton { "Split at frame" }, moveButton { "Move end" }, mergeButton { "Merge next" }, renameButton { "Rename" };
    juce::TextButton rebuildButton { "Build sections" }, cancelButton { "Cancel sections" };
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
