#include "LraPassagesPanel.h"
#include "StandaloneLraProposal.h"
#include "StandaloneProcessing.h"
#include "MeasurementTimeline.h"

namespace levelscope::standalone
{
LraPassagesPanel::LraPassagesPanel(StandaloneSessionModel& s,MeasurementDomain d):session(s),domain(d)
{
    setName("LRA passages and candidate estimates");
    for(auto* c:std::initializer_list<juce::Component*>{&source,&current,&refresh,&cancel,&assess,&previous,&next,&quiet,&status,&details,&list})addAndMakeVisible(c);
    source.onClick=[this]{domain=MeasurementDomain::source;};current.onClick=[this]{domain=MeasurementDomain::currentState;};
    refresh.onClick=[this]{if(control)control->cancel=true;key.clear();runningKey.clear();};cancel.onClick=[this]{if(control)control->cancel=true;result={};updateList();status.setText("View job cancelled. Audio and measurements unchanged; Refresh to retry.",juce::dontSendNotification);};
    quiet.onClick=[this]{page=0;updateList();};previous.onClick=[this]{--page;updateList();};next.onClick=[this]{++page;updateList();};
    assess.onClick=[this]{if(assess.isEnabled()&&key==selectedKey()&&onAssess)onAssess();};
    details.setMultiLine(true);details.setReadOnly(true);details.setScrollbarsShown(true);details.setFont(14.f);
    list.setRowHeight(32);status.setFont(14.f);setSize(1040,750);startTimerHz(10);timerCallback();
}
LraPassagesPanel::~LraPassagesPanel(){stopTimer();if(control)control->cancel=true;if(future.valid())future.wait();}
juce::String LraPassagesPanel::selectedKey()const
{
    const auto& measured=domain==MeasurementDomain::source?session.source:session.currentState;
    if(!measured.hasMeasurement||!measured.analysis||session.sourceContentInvalid||!session.sourceDocument.matchesFileStamps())return {};
    const auto c=session.measurementContext(domain);
    if(!measuredProgrammeForSnapshot(measured.analysis,c,domain,session.programmeApplicability).currentIdentityMatches)return {};
    return session.sessionId+";"+measured.analysis->identity+";"+juce::String(int(domain))+";"+session.sectionProposalRequest().binding();
}
void LraPassagesPanel::timerCallback()
{
    const auto nextKey=selectedKey();
    if(nextKey!=key)
    {
        if(control)control->cancel=true;key=nextKey;result={};page=0;updateList();
        status.setText(key.isEmpty()?"No current measured snapshot for this domain. Measure it first.":"Building measured passage view...",juce::dontSendNotification);
    }
    if(future.valid()&&future.wait_for(std::chrono::seconds(0))==std::future_status::ready)
    {
        auto completed=future.get();
        if(key==runningKey&&!control->isCancelled())
        {result=std::move(completed);updateList();status.setText(result.passages?(domain==MeasurementDomain::source?"SOURCE":"CURRENT")+juce::String(" inspected; audition selection unchanged. ")+result.passages->notice:result.error,juce::dontSendNotification);}
    }
    if(!future.valid()&&key.isNotEmpty()&&runningKey!=key)
    {
        runningKey=key;control=std::make_shared<JobControl>();const auto request=session.sectionProposalRequest();
        const auto snapshot=domain==MeasurementDomain::source?request.source:request.current;const auto which=domain;
        future=std::async(std::launch::async,[request,snapshot,which,token=control]
        {
            Built out;
            try
            {
                if(!snapshot)throw std::invalid_argument("Snapshot unavailable");
                auto a=proposal::analyzeProgramme(snapshot->detail,snapshot->identity,[&]{return token->isCancelled();});if(!a)return out;
                const auto& l=request.currentContext.profile.limits;
                out.passages=proposal::findLraPassages(a,l.hasMaxLoudnessRange?std::optional<double>(l.maxLoudnessRangeLu):std::nullopt,[&]{return token->isCancelled();});
                if(out.passages&&which==MeasurementDomain::currentState&&l.hasMaxLoudnessRange)
                {
                    processing::SectionGainPlan p;p.inputIdentity=sectionInputIdentity(request.currentContext);p.frames=a->frames;p.sampleRate=a->sampleRate;p.intervals={{0,p.frames,0}};
                    if(request.currentContext.processing.sectionGains)p=*request.currentContext.processing.sectionGains;
                    out.selection=proposal::selectLraTrims(*out.passages,p,request.policy,[&]{return token->isCancelled();});
                }
            }
            catch(const std::exception& e){out.error=e.what();}return out;
        });
    }
    const auto& limits=session.getSelectedTargetProfile().limits;
    assess.setEnabled(!future.valid()&&result.passages&&domain==MeasurementDomain::currentState&&limits.hasMaxLoudnessRange
        &&session.source.hasMeasurement&&session.currentState.hasMeasurement&&!session.currentStateNeedsRemeasurement()
        &&!session.job.isActive()&&!session.audioEditsBlocked()&&!session.sectionTrimDraft&&!session.persistenceLoading&&key==selectedKey());
    cancel.setEnabled(future.valid());
}
const std::vector<proposal::LraPassage>& LraPassagesPanel::rows()const
{static const std::vector<proposal::LraPassage> empty;return result.passages?(quiet.getToggleState()?result.passages->low:result.passages->high):empty;}
int LraPassagesPanel::getNumRows(){return int(std::min<size_t>(200,rows().size()-std::min(rows().size(),size_t(page)*200)));}
void LraPassagesPanel::updateList()
{
    page=std::clamp(page,0,std::max(0,int((rows().size()+199)/200)-1));list.updateContent();list.deselectAllRows();
    previous.setEnabled(page>0);next.setEnabled(size_t(page+1)*200<rows().size());
    juce::String text="Inspection never seeks, plays or changes processing. Estimates cannot satisfy measured requirements.\n";
    if(result.passages)
    {
        const auto& p=*result.passages;const auto& a=*p.analysis;
        text+="Measured LRA "+diagnosticMetricText(a.lra,"LU")+"; "+juce::String(a.lraStable?"stable":"short/unstable")+". Included "+juce::String(p.included)+", excluded "+juce::String(p.excluded)+", unavailable "+juce::String(p.unavailable)+" windows.\n";
        if(a.p10&&a.p95&&p.upperLevel)text+="P10 "+juce::String(*a.p10,9)+" / P95 "+juce::String(*a.p95,9)+" LUFS; provisional upper search level "+juce::String(*p.upperLevel,9)+" LUFS. Grouping uses full-resolution windows, not plot reduction.\n";
        text+="High/quiet rows are population support, not additive LRA contributions. "+juce::String(rows().size())+" regions; page "+juce::String(page+1)+".\n";
        if(result.selection)text+="Bounded preview summary (select a row for its exact values; complete assessment rationale is available in Review):\n"+result.selection->text().substring(0,12000)+"\nPreview only: complete measured assessment is required. Local attenuation may be followed by positive global gain.\n";
    }
    summary=text;details.setText(text,false);details.moveCaretToTop(false);repaint();
}
void LraPassagesPanel::paintListBoxItem(int row,juce::Graphics& g,int width,int height,bool selected)
{
    const auto n=size_t(page)*200+size_t(row);if(n>=rows().size()||!result.passages)return;
    if(selected)g.fillAll(juce::Colour(0xff345d65));g.setColour(juce::Colours::white);g.setFont(14.f);
    const auto& r=rows()[n];const auto rate=result.passages->analysis->sampleRate;
    g.drawText(juce::String(n+1)+"  "+MeasurementTimeline::timecode(r.begin,rate)+" – "+MeasurementTimeline::timecode(r.end,rate)
        +" | "+juce::String(r.windows)+" windows | score "+juce::String(r.score,3)+(r.paddedWindows?" | EOF padding":""),5,0,width-10,height,juce::Justification::centredLeft);
}
void LraPassagesPanel::selectedRowsChanged(int row)
{
    const auto n=size_t(page)*200+size_t(std::max(0,row));if(row<0||n>=rows().size()||!result.passages||key!=selectedKey())return;
    const auto& r=rows()[n];juce::String text="Selected measured population support ["+juce::String(r.begin)+", "+juce::String(r.end)+") audio-source frames; included windows "+juce::String(r.windows)+", EOF-padded "+juce::String(r.paddedWindows)+"; score "+juce::String(r.score,12)+" LU*windows, regional selected-window P95 "+juce::String(r.p95,12)+" LUFS.\n";
    if(result.passages->included)text+="Selected-window share of the included population "+juce::String(100.*double(r.windows)/result.passages->included,6)+"%. This is not a percentage contribution to programme LRA.\n";
    int names=0;for(const auto& a:session.sectionDocument.map().sections)if(a.begin<r.end&&r.begin<a.end&&names++<5)text+="Annotation context only: "+a.name+"\n";
    if(result.selection)for(const auto& chosen:result.selection->selected)if(chosen.passage.begin==r.begin&&chosen.passage.end==r.end)
        text+="Planning only: requested "+juce::String(chosen.requested,12)+" dB; reduction "+juce::String(chosen.reduction,12)+" dB; reserve "+juce::String(chosen.reserve,12)+"; uncovered overshoot "+juce::String(chosen.uncovered,12)+" dB.\n";
    details.setText(text+summary,false);details.moveCaretToTop(false);
    if(onInspect)onInspect(domain,r.begin,r.end);
}
void LraPassagesPanel::paint(juce::Graphics& g)
{
    g.fillAll(juce::Colour(0xff17212a));g.setColour(juce::Colours::lightgrey);g.setFont(13.f);
    g.drawText("Local gain preview only — grey applied / cyan estimated; global gain not yet selected",preview.withHeight(20),juce::Justification::centredLeft);
    if(!result.selection||!result.selection->available)return;
    auto bounds=preview.withTrimmedTop(22).toFloat();const auto& s=*result.selection;
    const auto draw=[&](const processing::SectionGainPlan& p,juce::Colour c)
    {
        const auto knots=gainOutline(&p,p.frames,0,double(p.frames),int(bounds.getWidth()));juce::Path path;bool first=true;
        for(auto k:knots){const float x=bounds.getX()+float(double(k.frame)/p.frames)*bounds.getWidth(),y=bounds.getY()+float((12-k.db)/24)*bounds.getHeight();if(first){path.startNewSubPath(x,y);first=false;}else path.lineTo(x,y);}
        g.setColour(c);g.strokePath(path,juce::PathStrokeType(1.5f));
    };
    draw(s.before,juce::Colours::grey);draw(s.plan,juce::Colours::cyan);
}
void LraPassagesPanel::resized()
{
    auto r=getLocalBounds().reduced(10),buttons=r.removeFromTop(32);
    for(auto* b:{&source,&current,&refresh,&cancel,&assess})b->setBounds(buttons.removeFromLeft(160).reduced(2));
    status.setBounds(r.removeFromTop(64));preview=r.removeFromTop(95);details.setBounds(r.removeFromBottom(185));
    auto nav=r.removeFromBottom(32);quiet.setBounds(nav.removeFromLeft(210));previous.setBounds(nav.removeFromLeft(150));next.setBounds(nav.removeFromLeft(150));list.setBounds(r.reduced(0,4));
}
}
