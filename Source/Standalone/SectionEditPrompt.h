#pragma once
#include "StandaloneSessionModel.h"
#include <juce_gui_basics/juce_gui_basics.h>

namespace levelscope::standalone
{
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
        if(deleting){requestDelete(owner,begin,end,std::move(completed));return;}
        if(!available()){completed(false,false);return;}
        auto next=session.manualState();next.map=after;
        if(!next.gains){completed(session.commitManualEdit(next,"Move section boundary"),false);return;}
        auto moved=next;const auto plateau=[&](int64_t f){for(const auto& r:moved.gains->intervals)if(r.begin<=f&&f<r.end)return r.db;return 0.;};
        section_edit::setRange(*moved.gains,std::min(oldBoundary,newBoundary),std::max(oldBoundary,newBoundary),plateau(newBoundary>oldBoundary?oldBoundary-1:oldBoundary));
        show(owner,before,"Move section boundary","Keep the gain in place, or move it with the section? Changes apply immediately. Undo: Command-Z.",
             {{"Keep gain in place",next},{"Move gain with section",moved}},std::move(completed));
    }
    void requestDelete(juce::Component& owner,int64_t begin,int64_t end,std::function<void(bool,bool)> completed)
    {
        if(!available()||begin<0||end<=begin||end>session.sectionDocument.map().frames){completed(false,false);return;}
        const auto before=session.manualState();
        if(before.map.sections.size()==1&&!before.gains){session.sectionTrimNotice="No local gain or internal section boundaries to remove.";completed(false,false);return;}
        const auto adjacent=[&](bool left)->std::optional<double>
        {
            for(const auto& r:before.map.sections)if(left?r.end==begin:r.begin==end)return section_edit::uniform(before.gains?&*before.gains:nullptr,r.begin,r.end);
            return {};
        };
        const auto left=adjacent(true),right=adjacent(false);
        const auto result=[&](double db)
        {
            auto next=before;auto p=next.gains.value_or(processing::SectionGainPlan{});
            if(!next.gains){p.inputIdentity=next.map.inputIdentity;p.sampleRate=next.map.sampleRate;p.frames=next.map.frames;p.intervals={{0,p.frames,0}};}
            section_edit::setRange(p,begin,end,db);next.gains=p;
            std::vector<proposal::ProgrammeSection> rows;
            for(auto r:next.map.sections)if(r.end<=begin)rows.push_back(r);else if(r.begin<begin){r.end=begin;rows.push_back(r);}
            const auto selected=rows.size();rows.push_back({next.map.nextId++,begin,end,begin==0&&end==next.map.frames?"Whole programme":"Section"});
            for(auto r:next.map.sections)if(r.begin>=end)rows.push_back(r);else if(r.end>end){r.begin=end;r.id=next.map.nextId++;rows.push_back(r);}
            auto i=selected;
            if(i>0&&section_edit::uniform(&p,rows[i-1].begin,rows[i-1].end)==std::optional<double>(db))
            {rows[i-1].end=rows[i].end;rows.erase(rows.begin()+ptrdiff_t(i));--i;}
            if(i+1<rows.size()&&section_edit::uniform(&p,rows[i+1].begin,rows[i+1].end)==std::optional<double>(db))
            {rows[i].end=rows[i+1].end;rows.erase(rows.begin()+ptrdiff_t(i+1));}
            next.map.sections=std::move(rows);return next;
        };
        if((!left||*left==0)&&(!right||*right==0))
        {completed(session.commitManualEdit(result(0),"Reset and merge selected sections"),true);return;}
        std::vector<Choice> choices{{"Set to 0 dB",result(0)}};
        if(left&&right&&*left==*right)choices.push_back({"Merge into neighbours ("+juce::String(*left,2)+" dB)",result(*left)});
        else
        {
            if(left&&*left!=0)choices.push_back({"Merge left ("+juce::String(*left,2)+" dB)",result(*left)});
            if(right&&*right!=0)choices.push_back({"Merge right ("+juce::String(*right,2)+" dB)",result(*right)});
        }
        show(owner,before.map,"Reset or merge sections","Set local gain in the selected passage to zero, or extend a neighbour's gain through it. Global gain stays unchanged. One action; Command-Z undoes it.",std::move(choices),std::move(completed));
    }
private:
    struct Choice {juce::String label;ManualEditState state;};
    StandaloneSessionModel& session;
    std::unique_ptr<juce::AlertWindow> dialog;
    bool available()
    {
        if(!dialog&&session.manualEditAllowed())return true;
        session.sectionTrimNotice="Stop playback and finish the current job or proposal before editing sections.";return false;
    }
    void show(juce::Component& owner,const proposal::SectionMap& before,const juce::String& title,const juce::String& message,
              std::vector<Choice> choices,std::function<void(bool,bool)> completed)
    {
        const auto opening=session.sessionId;const auto revision=session.revision;
        dialog=std::make_unique<juce::AlertWindow>(title,message,juce::MessageBoxIconType::QuestionIcon,&owner);
        for(size_t i=0;i<choices.size();++i)dialog->addButton(choices[i].label,int(i+1));
        dialog->addButton("Cancel",0,juce::KeyPress(juce::KeyPress::escapeKey));
        dialog->enterModalState(true,juce::ModalCallbackFunction::create(
            [safe=juce::Component::SafePointer<SectionEditPrompt>(this),opening,revision,before,choices=std::move(choices),completed](int answer)
            {
                if(!safe)return;auto& self=*safe;auto closed=std::move(self.dialog);closed->setVisible(false);
                if(answer<=0||size_t(answer)>choices.size()){completed(false,false);return;}
                if(opening!=self.session.sessionId||revision!=self.session.revision||self.session.sectionDocument.map().digest()!=before.digest())
                {self.session.sectionTrimNotice="Section edit cancelled: the programme or sections changed while the question was open.";completed(false,false);return;}
                const bool accepted=self.session.commitManualEdit(choices[size_t(answer-1)].state,choices[size_t(answer-1)].label);
                completed(accepted,true);
            }),false);
    }
};
}
