#pragma once
#include "StandaloneDiagnostics.h"
#include "StandaloneSessionModel.h"
#include "MeasurementTimeline.h"

namespace levelscope::standalone
{
class StandaloneDiagnosticsPanel final : public juce::Component, private juce::Timer, private juce::ListBoxModel
{
public:
    explicit StandaloneDiagnosticsPanel (const StandaloneSessionModel&);
    ~StandaloneDiagnosticsPanel() override;
    void resume();
    void suspend();
    void synchronize() { timerCallback(); }
    void paint (juce::Graphics&) override;
    void resized() override;
    std::function<void()> onOpenSections;
private:
    const StandaloneSessionModel& session; // observation only
    DiagnosticsController controller;
    std::shared_ptr<const MeasuredDiagnostics> displayed;
    DiagnosticAnnotations annotations;
    juce::String annotationKey;
    juce::ComboBox domainBox, rankingBox;
    juce::TextButton rebuildButton { "Build diagnostics" }, cancelButton { "Cancel diagnostics" };
    juce::TextButton sectionsButton { "Programme sections..." };
    juce::Label status;
    juce::TextEditor summary, details;
    MeasurementTimeline timeline;
    juce::ListBox ranking { "Full-resolution ranked diagnostic windows", this };
    juce::Rectangle<int> distributionBounds;
    bool suspended = false;
    void timerCallback() override;
    void refreshSummary();
    int getNumRows() override;
    void paintListBoxItem (int, juce::Graphics&, int, int, bool) override;
    void selectedRowsChanged (int) override;
    const proposal::DiagnosticWindow* windowForRow (int) const;
    juce::String rowText (int) const;
    void paintDistribution (juce::Graphics&);
};
}
