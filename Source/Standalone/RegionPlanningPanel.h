#pragma once
#include "TimelineAppearance.h"
#include "StandaloneProposalPreferences.h"

namespace levelscope::standalone
{
// Private exact-frame row editor. A pin is a closed range with equal endpoints;
// protected material instead uses half-open, nonempty intervals.
class RegionFrameRows final : public juce::Component, private juce::ListBoxModel
{
public:
    RegionFrameRows(juce::String title, bool pins, size_t maximum);
    std::vector<proposal::RegionEdgeRange> values;
    void load(std::vector<proposal::RegionEdgeRange>);
    juce::String pendingError() const;
    std::function<void()> onChanged;
    void resized() override;
private:
    bool allowPins,cleanPin=false;size_t limit;
    juce::String cleanFirst,cleanLast;
    void markClean();
    juce::Label title,status;
    juce::ListBox rows{"Exact source-frame rows",this};
    juce::TextEditor first,last;
    juce::ToggleButton pin{"Exact pin"};
    juce::TextButton add{"Add row"},update{"Update row"},remove{"Remove row"};
    int getNumRows() override {return int(values.size());}
    void paintListBoxItem(int,juce::Graphics&,int,int,bool) override;
    void selectedRowsChanged(int) override;
    void store(bool append);
};

class RegionPlanningPanel final : public juce::Component, private juce::Timer
{
public:
    explicit RegionPlanningPanel(StandaloneSessionModel&);
    static void showWindow(std::unique_ptr<ProposalWindow>&,StandaloneSessionModel&,juce::Component*,
                           std::function<void()>,std::function<void()>);
    std::function<void()> onChanged,onAssess;
    void resized() override;
    static std::optional<int64_t> exactFrame(juce::String);
    void paint(juce::Graphics& g) override {g.fillAll(palette::panel());}
    // Parsing/validation is shared by explicit Confirm and Confirm and assess.
    juce::String readDraft(std::optional<proposal::RegionPlanningSettings>&) const;
private:
    StandaloneSessionModel& session;
    juce::String runtime,input,baseDigest,feedback;
    proposal::SectionMap annotations;
    proposal::RegionPlanningSettings loadedSettings;
    juce::Viewport viewport;juce::Component body;
    juce::Label heading,methodLabel,modeLabel,sectionLabel,selectionLabel,edgeLabel,footnote;
    juce::ComboBox method,mode,section;
    juce::TextEditor selectionBegin,selectionEnd,startMin,startMax,endMin,endMax;
    juce::ToggleButton startRange{"Restrict processing start"},endRange{"Restrict processing end"},
        startPin{"Pin start"},endPin{"Pin end"},splitting{"Allow splitting"},refinement{"Allow boundary refinement"},
        fillGap{"Consider filling a short gap"},global{"Allow global trim of the whole programme"},compare{"Assess alternatives after a passing candidate"};
    juce::ToggleButton automaticMargins{"Allow automatic outward margins"};
    juce::Label eventNote;
    std::array<juce::Label,7> labels;std::array<juce::TextEditor,7> fields;
    RegionFrameRows cuts{"Internal nominal cuts - source frames; ranges or exact pins",true,4};
    RegionFrameRows protection{"No additional local-envelope change - [begin, end) source frames",false,128};
    juce::TextEditor status;
    juce::TextButton confirm{"Confirm settings"},assess{"Confirm and assess"},discard{"Discard edits / reload"},cancel{"Cancel"};
    void load();
    void commit(bool assessAfter);
    void timerCallback() override;
    void updateEnabled();
    juce::String geometryText(const proposal::RegionPlanningSettings&) const;
};
}
