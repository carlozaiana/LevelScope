#pragma once
#include "StandaloneSessionModel.h"
#include <juce_gui_basics/juce_gui_basics.h>

namespace levelscope::standalone
{
// Shared by the main lane, waveform boundary handles and the Sections window.
// A modal answer is accepted only for the exact opening, audio revision and map.
class SectionEditPrompt final : public juce::Component
{
public:
    explicit SectionEditPrompt(StandaloneSessionModel& s):session(s){}
    juce::AlertWindow* window() const {return dialog.get();}
    bool active() const {return bool(dialog);}
    void request(juce::Component& owner,const proposal::SectionMap& before,const proposal::SectionMap& after,
                 bool deleting,int64_t begin,int64_t end,int64_t oldBoundary,int64_t newBoundary,
                 std::function<void(bool,bool)> completed)
    {
        if(dialog||session.sectionTrimDraft||session.audioEditsBlocked()||session.job.isActive()||session.persistenceLoading)
        {session.sectionTrimNotice="Apply or Discard the current draft, stop playback and finish jobs before changing section geometry.";completed(false,false);return;}
        if(!session.currentStateDocument.sectionGains)
        {completed(session.editSectionAnnotations(before,after),false);return;}
        const auto opening=session.sessionId;const auto revision=session.revision;
        auto plan=*session.currentStateDocument.sectionGains;
        if(deleting)section_edit::setRange(plan,begin,end,0);
        else
        {
            // Move the processing boundary over the swept span. Existing plateaus
            // outside that span stay fixed; the canonical processor recomputes ramps.
            const auto plateau=[&](int64_t frame){for(const auto& r:plan.intervals)if(r.begin<=frame&&frame<r.end)return r.db;return 0.;};
            const auto left=plateau(oldBoundary-1),right=plateau(oldBoundary);
            section_edit::setRange(plan,std::min(oldBoundary,newBoundary),std::max(oldBoundary,newBoundary),newBoundary>oldBoundary?left:right);
        }
        dialog=std::make_unique<juce::AlertWindow>(deleting?"Delete section":"Change section boundary",
            deleting?"Keep its applied gain, or also remove the local gain in this passage? Removing gain restores 0 dB local; global gain stays unchanged.\n\nA section + gain edit remains a draft until Apply. Discard restores both.":
                "Keep the applied gain in place, or move the gain boundary with the section? The gain on the expanding side fills the swept passage. Other plateaus stay fixed; boundary smoothing is recomputed on Apply.\n\nA section + gain edit remains a draft until Apply. Discard restores both.",
            juce::MessageBoxIconType::QuestionIcon,&owner);
        dialog->addButton(deleting?"Keep applied gain":"Keep gain in place",1);
        dialog->addButton(deleting?"Remove applied gain":"Move gain with section",2);
        dialog->addButton("Cancel",0,juce::KeyPress(juce::KeyPress::escapeKey));
        dialog->enterModalState(true,juce::ModalCallbackFunction::create(
            [safe=juce::Component::SafePointer<SectionEditPrompt>(this),opening,revision,before,after,plan,completed](int answer)
            {
                if(!safe)return;auto& self=*safe;auto closed=std::move(self.dialog);closed->setVisible(false);
                if(answer==0){completed(false,false);return;}
                if(opening!=self.session.sessionId||revision!=self.session.revision||self.session.sectionDocument.map().digest()!=before.digest())
                {self.session.sectionTrimNotice="Section edit cancelled: the programme or section changed while the question was open.";completed(false,false);return;}
                const bool accepted=self.session.editSectionAnnotations(before,after,answer==2?std::optional<processing::SectionGainPlan>(plan):std::nullopt);
                if(!accepted)self.session.sectionTrimNotice="Section edit unavailable or stale. Finish jobs and resolve any unapplied proposal first.";
                completed(accepted,answer==2);
            }),false);
    }
private:
    StandaloneSessionModel& session;
    std::unique_ptr<juce::AlertWindow> dialog;
};
}
