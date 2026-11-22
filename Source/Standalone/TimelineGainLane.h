#pragma once
#include "StandaloneSessionModel.h"
#include "StandaloneRegionProposal.h"
#include "SectionEditPrompt.h"
#include <limits>
#include <juce_gui_basics/juce_gui_basics.h>

namespace levelscope::standalone
{
// Message-thread view of the existing M12 section draft. Owns no audio or transport.
class TimelineGainLane final : public juce::Component
{
public:
    static constexpr int height=170;
    explicit TimelineGainLane(StandaloneSessionModel&);
    void synchronize(const juce::String& confirmedProgramme);
    void refreshProcessing(){const auto identity=programme;synchronize(identity);}
    void setViewport(double begin,double end,int left,int width);
    void selectFrame(int64_t frame,bool notify);
    void selectAt(int64_t frame,int y);
    void toggleSection(int64_t frame);
    std::vector<std::pair<int64_t,int64_t>> selectedSpans() const;
    bool deleteSelected();
    bool removeSelectedBoundaries();
    void selectRange(int64_t begin,int64_t end);
    void editReadout(int64_t frame);
    size_t selectedCount() const {return selectedIds.size();}
    std::optional<size_t> appliedBoundaryHandle(juce::Point<int>) const;
    juce::String appliedBoundaryHelp() const;
    bool selectsSections() const {return selectionKind.getSelectedId()==1;}
    bool commitPending(){return stageValue();}
    SectionEditPrompt& editPrompt(){return sectionEdit;}
    void clearDeleteSelection(){deleteArmed=false;}
    void inspectRegionAttempt(size_t);
    void inspectApplied(int64_t frame){selectionKind.setSelectedId(2,juce::sendNotificationSync);selectFrame(frame,true);}
    bool diagnosticRegionInspection() const;
    int preferredHeight() const {return height+(hasRegionLayers()?76:0);}
    juce::Rectangle<int> plotBounds() const {return {plotLeft,47+(hasRegionLayers()?76:0),plotWidth,48};}
    juce::Rectangle<int> inspectionBounds() const {return {plotLeft,47,plotWidth,plotBounds().getBottom()+42-47};}
    juce::String inspectionText() {return detail.getTooltip();}
    std::function<void(int64_t,int64_t,int64_t)> onSelection;
    std::function<void()> onChanged;
    void resized() override;
    void paint(juce::Graphics&) override;
private:
    StandaloneSessionModel& session;
    SectionEditPrompt sectionEdit{session};
    juce::Label legend,notice,detail,dbLabel;
    juce::TextButton previous{"<"},next{">"},apply{"Apply gains"},discard{"Discard gains"};
    juce::TextButton undoDelete{"Undo delete"};
    std::optional<proposal::SectionMap> deletedAnnotation;
    uint64_t deletedRevision=0;
    bool deleteArmed=false,rangeSelection=false;
    juce::TextEditor gain;
    juce::ComboBox selectionKind;
    juce::String programme,runtime,input,supportError,draftError,feedback;
    uint64_t audioRevision=0,mapRevision=0,editSerial=0,boundarySerial=0,appliedBoundarySerial=0;
    bool initialized=false,hasDraft=false,fresh=false,geometryMatches=true,dirty=false,usableMap=false,measured=false;
    std::vector<uint64_t> selectedIds;
    size_t row=0;
    int64_t cursor=0;
    double viewBegin=0,viewEnd=1;
    int plotLeft=0,plotWidth=0;
    std::optional<processing::SectionGainPlan> applied,draft;
    std::shared_ptr<const SectionProposalResult> automatic;
    std::optional<ProcessingPlan> proposed;
    std::vector<proposal::ProgrammeSection> lraRows;
    bool automaticStale=false;
    size_t regionAttempt=std::numeric_limits<size_t>::max();
    const RegionAttempt* inspectedAttempt() const;
    bool hasRegionLayers() const {return bool(session.regionPlanning)||(automatic&&bool(automatic->regionReview()));}
    void rebuildProposalRows();
    void paintRegionLayers(juce::Graphics&);
    void paintReadouts(juce::Graphics&);
    const std::vector<proposal::ProgrammeSection>& rows() const;
    void adjacent(int direction);
    const proposal::ProgrammeSection* selected() const;
    void showSelection();
    void updateControls();
    bool stageValue();
    bool bindingMatches() const;
    void afterAction();
};
}
