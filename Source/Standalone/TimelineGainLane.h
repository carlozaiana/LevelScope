#pragma once
#include "StandaloneSessionModel.h"
#include <juce_gui_basics/juce_gui_basics.h>

namespace levelscope::standalone
{
// Message-thread view of the existing M12 section draft. Owns no audio or transport.
class TimelineGainLane final : public juce::Component
{
public:
    static constexpr int height=128;
    explicit TimelineGainLane(StandaloneSessionModel&);
    void synchronize(const juce::String& confirmedProgramme);
    void setViewport(double begin,double end,int left,int width);
    void selectFrame(int64_t frame,bool notify);
    juce::Rectangle<int> plotBounds() const {return {plotLeft,47,plotWidth,48};}
    juce::String inspectionText() {return detail.getTooltip();}
    std::function<void(int64_t,int64_t,int64_t)> onSelection;
    std::function<void()> onChanged;
    void resized() override;
    void paint(juce::Graphics&) override;
private:
    StandaloneSessionModel& session;
    juce::Label legend,notice,detail,dbLabel;
    juce::TextButton previous{"<"},next{">"},stage{"Set draft"},apply{"Apply gains"},discard{"Discard gains"};
    juce::TextEditor gain;
    juce::String programme,runtime,input,supportError,draftError,feedback;
    uint64_t audioRevision=0,mapRevision=0,editSerial=0;
    bool initialized=false,hasDraft=false,fresh=false,geometryMatches=true,dirty=false,usableMap=false,measured=false;
    size_t row=0;
    int64_t cursor=0;
    double viewBegin=0,viewEnd=1;
    int plotLeft=0,plotWidth=0;
    std::optional<processing::SectionGainPlan> applied,draft;
    std::shared_ptr<const SectionProposalResult> automatic;
    std::optional<ProcessingPlan> proposed;
    std::vector<proposal::ProgrammeSection> lraRows;
    bool automaticStale=false;
    const std::vector<proposal::ProgrammeSection>& rows() const;
    void adjacent(int direction);
    const proposal::ProgrammeSection* selected() const;
    void showSelection();
    void updateControls();
    void stageValue();
    bool bindingMatches() const;
    void afterAction();
};
}
