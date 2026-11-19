#include "RegionProposalReviewPanel.h"
#include "RegionReviewFormatting.h"

namespace levelscope::standalone
{
RegionProposalReviewPanel::RegionProposalReviewPanel(std::shared_ptr<const SectionProposalResult> r,std::function<juce::String()> s):result(std::move(r)),statusText(std::move(s))
{
    setName("Region attempts - forecasts and complete private measurements");
    for(auto* c:std::initializer_list<juce::Component*>{&attempts,&regions,&label,&authority,&inspect,&useAlternative,&coverage})addAndMakeVisible(c);
    label.setText("Browse measured attempts for inspection. Recommendation and explicit selection are separate; inspection never changes audio.",juce::dontSendNotification);
    coverage.setName("Measured alternative coverage and gains");coverage.setMultiLine(true);coverage.setReadOnly(true);coverage.setScrollbarsShown(true);coverage.setFont(13.f);
    coverage.setColour(juce::TextEditor::backgroundColourId,palette::card());useAlternative.onClick=[this]{useMeasuredAlternative();};
    authority.setFont(13.f);authority.setMinimumHorizontalScale(1.f);
    attempts.setName("Measured region attempt");regions.setName("Proposed processing root");attempts.addItem("Complete explanation - all attempts, restrictions and provenance",1);
    const auto review=result->regionReview();reports.resize(review?review->attempts.size()+1:1);
    if(review)for(size_t n=0;n<review->attempts.size();++n)
    {
        const auto& a=review->attempts[n];attempts.addItem("Attempt "+juce::String(n+1)+(a.evaluation.verdict==proposal::Verdict::pass?" - measured PASS":" - measured failure / unavailable")+(review->eventMode?" | "+eventAttemptCoverage(a):n==review->retainedAttempt?" (retained)":""),int(n+2));
    }
    attempts.setSelectedId(1,juce::dontSendNotification);attempts.onChange=[this]{chooseAttempt();inspectSelection();};
    inspect.onClick=[this]{inspectSelection();};chooseAttempt();setSize(1100,920);startTimerHz(4);
}
juce::String RegionProposalReviewPanel::attemptText(size_t n) const
{
    const auto review=result->regionReview();const auto& a=review->attempts[n];
    juce::String text="Complete private measured attempt "+juce::String(n+1)+"\nPlan identity "+a.plan.identity+"\nSnapshot "+a.snapshotIdentity+"\n";
    if(review->eventMode)
    {
        text+="EVENT TREATMENT - evidence, not a delivery rule or musical-quality verdict:\n"+eventAttemptCoverage(a)+"\n"+a.coverage.text()+"\n";
        if(review->eventEvidence)for(const auto& g:review->eventEvidence->groups)
            if(std::any_of(a.coverage.groups.begin(),a.coverage.groups.end(),[&](const auto& treatment){return treatment.groupId==g.id;}))
                text+="Event "+g.id+"; onset bracket ["+juce::String(g.onset.begin)+", "+juce::String(g.onset.end)+"), last attack ["+juce::String(g.lastAttack.begin)+", "+juce::String(g.lastAttack.end)+"), supported/uncertain decay ["+juce::String(g.decay.begin)+", "+juce::String(g.decay.end)+").\nNeed windows "+juce::String(g.firstNeedWindow)+" to "+juce::String(g.lastNeedWindow)+"; channel witness mask "+juce::String(g.witnesses)+"; "+g.reason+"\n";
        text+="\nExact permissions, requested/effective areas, pins and omitted comparisons remain in Complete explanation.\n";
    }
    text+="ACTUAL measured rule result:\n"+evaluationText(a.evaluation)+"\n";
    const char* names[]={"IL LUFS","True peak dBTP","LRA LU","Max M LUFS","Max S LUFS","Sample peak dBFS"};
    for(size_t k=0;k<a.measured.size();++k)text+=juce::String(names[k])+": "+(a.measured[k].value?regionExact(*a.measured[k].value):"unavailable (status "+juce::String(int(a.measured[k].status))+")")+"\n";
    const auto forecast=[&](const char* label,std::optional<double> value){text+=juce::String(label)+": "+(value?regionExact(*value):"unknown")+"\n";};
    text+="\nFORECASTS - approximate, never acceptance authority:\n";forecast("IL LUFS",a.prediction.integrated);forecast("LRA LU",a.prediction.lra);forecast("P10 LUFS",a.prediction.p10);forecast("P95 LUFS",a.prediction.p95);forecast("Sampled M LUFS (not sample-resolution maximum)",a.prediction.momentarySampled);forecast("Sampled S LUFS (not sample-resolution maximum)",a.prediction.shortTermSampled);forecast("TP upper bound (inconclusive if above ceiling)",a.prediction.truePeakUpperBound);
    text+=a.prediction.limitations+"\nGate-sensitive: "+juce::String(a.prediction.gateSensitive?"yes":"no")+"; EOF uncertainty: "+juce::String(a.prediction.eofUncertain?"yes":"no")+".\n";
    text+="Forecast gate entry/exit: I "+juce::String(a.prediction.integratedEntered)+" / "+juce::String(a.prediction.integratedLeft)+"; LRA "+juce::String(a.prediction.lraEntered)+" / "+juce::String(a.prediction.lraLeft)+".\nActual minus forecast (IL/LRA/P10/P95/sampled M):";
    for(auto e:a.signedPredictionError)text+=" "+(e?regionExact(*e):"unknown");text+="\nEmpirical errors are not bounds or rule tolerances.\n\n";
    text+=regionActivityDescription(a.plan,review->constraints.sampleRate,review->constraints.frames);
    for(size_t k=0;k<a.plan.assignments.size();++k)
    {
        const auto& part=a.plan.assignments[k];const auto root=part.geometry.root<review->constraints.roots.size()?review->constraints.roots[part.geometry.root].id:"unavailable";
        text+="\nRoot "+root+": "+part.geometry.reason+"\nNominal partition reductions (additional attenuation dB):";for(auto d:part.reductions)text+=" "+regionExact(d);
        text+="\nRequested / achieved reserve / uncovered deficit per part:\n";
        for(size_t j=0;j<part.requested.size()&&j<part.reserve.size()&&j<part.uncovered.size();++j)text+=juce::String(j+1)+": "+regionExact(part.requested[j])+" / "+regionExact(part.reserve[j])+" / "+regionExact(part.uncovered[j])+" dB\n";
    }
    text+="\nChoose Complete explanation for every alternative, baseline, restriction and search omission. No acceptable candidate in a bounded search does not prove no solution exists.\n";
    return text;
}
void RegionProposalReviewPanel::chooseAttempt()
{
    regions.clear(juce::dontSendNotification);regions.addItem("Whole plan",1);regions.setSelectedId(1,juce::dontSendNotification);
    const auto review=result->regionReview();const int n=attempts.getSelectedId()-2;
    const bool actualSelected=review&&n>=0&&size_t(n)<review->attempts.size();
    if(actualSelected)
    {
        const auto& a=review->attempts[size_t(n)];
        for(size_t k=0;k<a.plan.assignments.size();++k)
        {
            const auto& part=a.plan.assignments[k];
            regions.addItem("Root "+juce::String(part.geometry.root+1)+" support ["+juce::String(part.geometry.support.begin)+", "+juce::String(part.geometry.support.end)+")",int(k+2));
        }
    }
    inspect.setEnabled(actualSelected||(review&&n<0&&review->retainedAttempt<review->attempts.size()));regions.setEnabled(actualSelected);
    const size_t slot=actualSelected?size_t(n+1):0;
    if(report)report->setVisible(false);
    if(!reports[slot])
    {
        reports[slot]=std::make_unique<ProposalInspectionPanel>(actualSelected?attemptText(size_t(n)):result->text(),statusText,pageBytes);
        addChildComponent(*reports[slot]);
    }
    report=reports[slot].get();report->setVisible(true);refreshSelection();resized();
}
void RegionProposalReviewPanel::inspectSelection()
{
    const auto review=result->regionReview();if(!review||!onInspect)return;
    const int selected=attempts.getSelectedId()-2;
    const size_t n=selected>=0?size_t(selected):review->retainedAttempt;if(n>=review->attempts.size())return;
    const auto& plan=review->attempts[n].plan;proposal::RegionSpan span{0,review->constraints.frames};const int part=regions.getSelectedId()-2;
    if(part>=0&&size_t(part)<plan.assignments.size())span=plan.assignments[size_t(part)].geometry.support;
    onInspect(n,span.begin,span.end);refreshSelection();
}
void RegionProposalReviewPanel::configureMeasuredSelection(bool qualified,std::function<juce::String()> blocked,Selection select)
{selectionQualified=qualified;selectionBlocked=std::move(blocked);selectMeasured=std::move(select);refreshSelection();resized();}
void RegionProposalReviewPanel::refreshSelection()
{
    const auto review=result->regionReview();const int n=attempts.getSelectedId()-2;
    const bool event=review&&review->eventMode,actual=review&&n>=0&&size_t(n)<review->attempts.size();
    useAlternative.setVisible(event);coverage.setVisible(event);
    juce::String blocked;
    if(!selectionQualified)blocked="Measured-alternative selection awaits integrated authority qualification.";
    else if(!selectMeasured||!selectionBlocked)blocked="Selection is unavailable in this inspection snapshot.";
    else blocked=selectionBlocked();
    if(blocked.isEmpty()&&!actual)blocked="Choose an actual measured attempt to inspect its gains and coverage.";
    if(blocked.isEmpty()&&(!review->attempts[size_t(n)].receipt||review->attempts[size_t(n)].attemptIdentity.isEmpty()||review->attempts[size_t(n)].evaluation.verdict!=proposal::Verdict::pass))blocked="Only an owned actual measured pass may be selected.";
    if(blocked.isEmpty()&&review->selectedAttempt==size_t(n))blocked=result->explicitlySelectedAlternative()?"This measured alternative is selected. Apply remains explicit in the main window.":"This measured recommendation is already selected. Apply remains explicit in the main window.";
    useAlternative.setEnabled(event&&blocked.isEmpty());useAlternative.setTooltip(blocked.isNotEmpty()?blocked:"Select this measured plan without applying it. Its coverage limitations stay visible.");
    juce::String text=actual?regionAttemptAuthority(*result,size_t(n)):"Complete explanation: all attempts, baseline, restrictions and evidence.";
    if(event)text+="\nRecommended: "+(review->recommendedAttempt?"attempt "+juce::String(*review->recommendedAttempt+1):"none")+"; selected: "+(review->selectedAttempt?"attempt "+juce::String(*review->selectedAttempt+1):"none")+". Browsing does not change selection.";
    else text+=" Inspection never auditions a private candidate.";
    if(statusText)text+="\n"+statusText();
    authority.setText(text,juce::dontSendNotification);authority.setTooltip(text);
    juce::String details=selectionNotice;
    if(event&&actual)
    {
        const auto& a=review->attempts[size_t(n)];details+=eventAttemptCoverage(a)+"\n"+eventCoverageSummary(a.coverage)+"\n";
        details+="Actual measured rules: "+evaluationText(a.evaluation)+"\nGlobal gain "+regionExact(a.plan.globalDb)+" dB; whole-programme delta "+regionExact(a.plan.activity.globalDelta)+" dB. Additional local attenuation per part:";
        for(const auto& assignment:a.plan.assignments)for(auto d:assignment.reductions)details+=" "+regionExact(d)+" dB;";
        details+="\nExact interval gains, ramp support, all measured values and forecast limitations are in the attempt report below.";
    }
    if(event&&blocked.isNotEmpty())details+="\n"+blocked;
    if(coverage.getText()!=details){coverage.setText(details,false);coverage.moveCaretToTop(false);}
    if(report)report->refreshStatus();
}
void RegionProposalReviewPanel::useMeasuredAlternative()
{
    refreshSelection();if(!useAlternative.isEnabled())return;
    const auto review=result->regionReview();const auto n=size_t(attempts.getSelectedId()-2);
    const auto selected=selectMeasured(review->resultIdentity,review->attempts[n].attemptIdentity);
    if(!selected){selectionNotice="Selection rejected; audio and measurements are unchanged. Review the current binding/status below.\n";refreshSelection();return;}
    result=selected;selectionNotice.clear();if(report==reports[0].get())report=nullptr;reports[0].reset();chooseAttempt();
}
void RegionProposalReviewPanel::resized()
{
    auto r=getLocalBounds().reduced(10);label.setBounds(r.removeFromTop(28));attempts.setBounds(r.removeFromTop(32).reduced(1));
    auto line=r.removeFromTop(32);inspect.setBounds(line.removeFromRight(180).reduced(1));if(useAlternative.isVisible())useAlternative.setBounds(line.removeFromRight(245).reduced(1));regions.setBounds(line.reduced(1));authority.setBounds(r.removeFromTop(76));if(coverage.isVisible())coverage.setBounds(r.removeFromTop(150).reduced(0,3));if(report)report->setBounds(r);
}
}
