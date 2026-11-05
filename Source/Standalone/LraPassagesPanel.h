#pragma once
#include "StandaloneSessionModel.h"
#include "TimelineGainOutline.h"
#include <juce_gui_basics/juce_gui_basics.h>
#include <future>

namespace levelscope::standalone
{
class LraPassagesPanel final : public juce::Component,private juce::Timer,private juce::ListBoxModel
{
public:
    LraPassagesPanel(StandaloneSessionModel&,MeasurementDomain);
    ~LraPassagesPanel() override;
    std::function<void()> onAssess;
    std::function<void(MeasurementDomain,int64_t,int64_t)> onInspect;
    void resized() override;
    void paint(juce::Graphics&) override;
private:
    struct Built {std::shared_ptr<const proposal::LraPassages> passages;std::optional<proposal::LraTrimSelection> selection;juce::String error;};
    StandaloneSessionModel& session;MeasurementDomain domain;
    juce::String key,runningKey,summary;
    std::shared_ptr<JobControl> control;
    std::future<Built> future;
    Built result;
    int page=0;
    juce::TextButton source{"Inspect Source"},current{"Inspect Current"},refresh{"Refresh passages"},cancel{"Cancel view job"},assess{"Assess candidate"},previous{"Previous 200"},next{"Next 200"};
    juce::ToggleButton quiet{"Quiet population"};
    juce::Label status;
    juce::TextEditor details;
    juce::ListBox list{"Measured LRA passages",this};
    juce::Rectangle<int> preview;
    void timerCallback() override;
    juce::String selectedKey() const;
    const std::vector<proposal::LraPassage>& rows() const;
    int getNumRows() override;
    void paintListBoxItem(int,juce::Graphics&,int,int,bool) override;
    void selectedRowsChanged(int) override;
    void updateList();
};
}
