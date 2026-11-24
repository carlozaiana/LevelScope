#pragma once
#include "Core/Proposal/LraTimelineDiagnostics.h"
#include "TimelineLayout.h"
#include "TimelineViewPreferences.h"
#include "TimelineAppearance.h"

namespace levelscope::standalone
{
// Presentation of Logic's immutable measured index. No gate, percentile, target
// classification, proposal prediction or PCM computation takes place here.
class LraTimelineLanes final : public juce::Component, public juce::TooltipClient
{
public:
    LraTimelineLanes()
    {
        setName("Measured LRA diagnostic lanes");
        for(size_t n=0;n<dividers.size();++n)
        {
            dividers[n]=std::make_unique<TimelineDivider>(titles[int(n)]+" height");addAndMakeVisible(*dividers[n]);
            dividers[n]->move=[this,n](int delta){if(!heights)return;heights->lra[n]+=delta;heights->clamp();if(onHeightChanged)onHeightChanged();};
            dividers[n]->reset=[this,n]{if(!heights)return;heights->lra[n]=TimelineHeights{}.lra[n];if(onHeightChanged)onHeightChanged();};
            dividers[n]->finished=[this]{if(onSaved)onSaved();};
        }
    }
    void setHeights(TimelineHeights* h){heights=h;resized();}
    std::function<void()> onHeightChanged,onSaved,onInvalidated;
    std::function<void(const juce::String&)> onInspect;
    std::function<void(const juce::MouseEvent&,const juce::MouseWheelDetails&)> onWheel;
    int preferredHeight()const{int h=0;for(auto n:heights?heights->lra:TimelineHeights{}.lra)h+=n+TimelineLayout::dividerHeight;return h;}
    TimelineDivider& divider(size_t n){return *dividers.at(n);}
    juce::Rectangle<int> laneBounds(size_t n)const{return boxes.at(n);}
    void unavailable(const juce::String& domain,const juce::String& reason)
    {
        if(!index&&domainText==domain&&unavailableText==reason)return;
        invalidateInspection();index.reset();columns.clear();identity.clear();target={};domainText=domain;unavailableText=reason;selected.reset();repaint();
    }
    void measured(std::shared_ptr<const proposal::lra_timeline::Index> value,const proposal::TargetProfile& profile,proposal::Applicability applicability,const juce::String& domain)
    {
        if(!value){unavailable(domain,"Measured diagnostics are being prepared.");return;}
        const auto key=profile.digest()+";"+juce::String(int(applicability));
        if(index==value&&identity==key&&domainText==domain)return;
        invalidateInspection();index=std::move(value);identity=key;domainText=domain;unavailableText.clear();target=index->target(profile,applicability);selected.reset();query();
    }
    void viewport(double first,double last,int x,int width)
    {
        if(begin==first&&end==last&&left==x&&plotWidth==width)return;
        begin=first;end=last;left=x;plotWidth=width;query();
    }
    void markers(int64_t cursor,std::optional<std::pair<double,double>> range)
    {if(playhead==cursor&&selection==range)return;playhead=cursor;selection=range;repaint();}
    const auto& displayedColumns()const{return columns;}
    const auto& measuredIndex()const{return index;}
    const auto& measuredTarget()const{return target;}
    size_t queryCount()const{return queries;}
    std::optional<int64_t> inspectedFrame()const{return selected;}
    juce::String status(size_t n)const
    {
        using namespace proposal::lra_timeline;
        if(!index)return domainText+" - "+unavailableText;
        const auto& m=index->measured();const auto unstable=m.p10&&m.p95&&m.lra.value?(m.lraStable?juce::String():" Available but unstable LRA population (under 60 s)."):
            m.lra.status==measurement::Status::insufficientDuration?juce::String(" Insufficient programme duration; no percentile population."):
            juce::String(" No eligible percentile population; gate evidence only.");
        const auto level=[](std::optional<double> v){return v?juce::String(*v,1)+" LUFS":juce::String("unavailable");};
        if(n==0)
        {
            const auto prefix=domainText+" measured LRA "+(m.lra.value?juce::String(*m.lra.value,1)+" LU":"unavailable")+". ";
            switch(target.state)
            {
                case TargetState::noCeiling:return prefix+"Set an LRA target to highlight potential passages.";
                case TargetState::inapplicable:return prefix+"Selected target is not applicable to this programme.";
                case TargetState::unavailable:return prefix+"Insufficient eligible data for target guidance.";
                case TargetState::unstable:return prefix+"Population is unstable; no target guidance.";
                case TargetState::within:return prefix+"Within the selected "+juce::String(target.ceilingLu.value_or(0),1)+" LU ceiling.";
                case TargetState::exceeds:return prefix+"Target "+juce::String(target.ceilingLu.value_or(0),1)+" LU; highlighted windows above "+level(target.upperLevelLufs)+". Inspection aid only.";
            }
        }
        if(n==1)return domainText+" | P95 "+level(m.p95)+" | Amber: above P95; hatched: tied upper ranks. Not 5% of audio duration."+unstable;
        return domainText+" | P10 "+level(m.p10)+" | Blue: below P10; light grey: absolute gate; dark grey: relative gate; dotted: zero / too short."+unstable;
    }
    juce::String inspectionAt(int64_t frame,size_t lane)const
    {
        using namespace proposal::lra_timeline;
        juce::String text=titles[int(lane)]+"\n"+status(lane)+"\n";
        if(!index)return text;
        const auto& m=index->measured();
        text+="Measured 3 s LRA windows touching "+tc(frame,m.sampleRate)+"; overlapping support, not isolated causes.\n";
        text+="Relative LRA gate: "+(m.lraGate?juce::String(*m.lraGate,2)+" LUFS":"unavailable")+". Gate equality is included.\n";
        for(const auto& w:index->inspect(frame,target))
        {
            juce::StringArray categories;
            if(w.flags&eligible)categories.add("eligible");if(w.flags&targetPassage)categories.add("potential target passage");
            if(w.flags&aboveP95)categories.add("above P95");if(w.flags&upperRankTie)categories.add("P95 tie in upper ranks");
            if(w.flags&belowP10)categories.add("below P10");if(w.flags&atP10)categories.add("at P10");if(w.flags&atP95)categories.add("at P95");
            if(w.flags&absoluteRejected)categories.add("below absolute gate");if(w.flags&relativeRejected)categories.add("below relative gate");
            if(w.flags&zeroEnergy)categories.add("zero energy");if(w.flags&insufficientDuration)categories.add("insufficient duration");
            text+=tc(w.begin,m.sampleRate)+" - "+tc(w.end,m.sampleRate)+" | samples ["+juce::String(w.begin)+", "+juce::String(w.end)+") | "+(w.lufs?juce::String(*w.lufs,2)+" LUFS":"no finite loudness")+" | "+categories.joinIntoString(", ");
            if(w.flags&padded)text+=" | terminal padding to "+tc(w.measuredEnd,m.sampleRate)+"; display clipped to programme";
            text+="\n";
        }
        text+="Inspection only: these are measured windows, not per-section compliance verdicts. Low percentile rank alone does not establish inaudibility; global gain can change absolute level without changing percentile rank.";
        return text;
    }
    juce::String getTooltip()override
    {
        const auto p=getMouseXYRelative();const auto n=laneAt(p.y);return n?status(*n)+"\nClick to inspect the measured overlapping windows. Drag the divider below this lane to resize; double-click resets its height.":juce::String();
    }
    void mouseDown(const juce::MouseEvent& e)override
    {
        const auto n=laneAt(e.y);if(!n||!index||plotWidth<2||end<=begin||e.x<left||e.x>=left+plotWidth)return;
        selected=std::min(index->measured().frames-1,timelineFrameAt(e.position.x,left,plotWidth,begin,end,index->measured().frames));
        if(onInspect)onInspect(inspectionAt(*selected,*n));repaint();
    }
    void mouseWheelMove(const juce::MouseEvent& e,const juce::MouseWheelDetails& w)override{if(onWheel)onWheel(e,w);}
    void resized()override
    {
        int y=0;for(size_t n=0;n<3;++n){const auto h=heights?heights->lra[n]:TimelineHeights{}.lra[n];boxes[n]={0,y,getWidth(),h};y+=h;dividers[n]->setBounds(0,y,getWidth(),TimelineLayout::dividerHeight);y+=TimelineLayout::dividerHeight;}
    }
    void paint(juce::Graphics& g)override
    {
        using namespace proposal::lra_timeline;
        g.fillAll(palette::panel());
        const auto x=[&](double f){return float(left+(f-begin)*std::max(1,plotWidth-1)/(end-begin));};
        for(size_t n=0;n<3;++n)
        {
            const auto b=boxes[n];g.setColour(juce::Colours::white);g.setFont(12.f);g.drawText(titles[int(n)],8,b.getY(),getWidth()-16,18,juce::Justification::left);
            g.setColour(juce::Colours::lightgrey);g.setFont(11.f);g.drawText(status(n),8,b.getY()+18,getWidth()-16,16,juce::Justification::left);
            const juce::Rectangle<int> plot{left,b.getY()+35,plotWidth,std::max(5,b.getHeight()-36)};g.setColour(palette::card());g.fillRect(plot);
            if(!index||columns.empty()||end<=begin)continue;
            juce::Graphics::ScopedSaveState save(g);g.reduceClipRegion(plot);
            for(size_t i=0;i<columns.size();++i)
            {
                const auto f=columns[i].flags;const float px=left+float(i)*plotWidth/columns.size(),w=std::max(1.f,float(plotWidth)/columns.size());
                const auto bar=[&](juce::Colour colour,int row,int count){g.setColour(colour);const float h=float(plot.getHeight())/count;g.fillRect(px,plot.getY()+row*h,w,std::max(1.f,h-1));};
                if(n==0&& (f&targetPassage))bar(juce::Colour(0xffffb65c),0,1);
                if(n==1){if(f&aboveP95)bar(juce::Colour(0xffffb65c),0,2);if(f&upperRankTie){bar(juce::Colour(0xffb5a5d7),1,2);if(i%4==0){g.setColour(palette::card());g.drawVerticalLine(int(px),float(plot.getCentreY()),float(plot.getBottom()));}}}
                if(n==2){if(f&belowP10)bar(juce::Colour(0xff62c9ff),0,4);if(f&absoluteRejected)bar(juce::Colour(0xffb4bcc6),1,4);if(f&relativeRejected)bar(juce::Colour(0xff758392),2,4);if((f&(zeroEnergy|insufficientDuration))&&i%3==0)bar(juce::Colour(0xff89909b),3,4);}
            }
            if(selection&&selection->second>begin&&selection->first<end){g.setColour(juce::Colours::yellow.withAlpha(.15f));g.fillRect(x(std::max(begin,selection->first)),float(plot.getY()),std::max(1.f,x(std::min(end,selection->second))-x(std::max(begin,selection->first))),float(plot.getHeight()));}
            for(auto f:{std::optional<int64_t>(playhead),selected})if(f&&*f>=begin&&*f<end){g.setColour(f==selected?juce::Colours::yellow:juce::Colours::white.withAlpha(.6f));g.drawVerticalLine(int(x(double(*f))),float(plot.getY()),float(plot.getBottom()));}
        }
    }
private:
    const juce::StringArray titles{"LRA - potential passages to inspect","LRA - upper percentile tail","LRA - quiet and gate-excluded material"};
    TimelineHeights* heights=nullptr;
    std::array<juce::Rectangle<int>,3> boxes;
    std::array<std::unique_ptr<TimelineDivider>,3> dividers;
    std::shared_ptr<const proposal::lra_timeline::Index> index;
    proposal::lra_timeline::Target target;
    std::vector<proposal::lra_timeline::Column> columns;
    juce::String identity,domainText="Original",unavailableText="No completed measurement.";
    double begin=0,end=1;int left=55,plotWidth=1;int64_t playhead=0;size_t queries=0;
    std::optional<int64_t> selected;std::optional<std::pair<double,double>> selection;
    void invalidateInspection(){if(selected&&onInvalidated)onInvalidated();selected.reset();}
    void query(){columns=index?index->columns(int64_t(std::floor(begin)),int64_t(std::ceil(end)),plotWidth,target):std::vector<proposal::lra_timeline::Column>{};++queries;repaint();}
    std::optional<size_t> laneAt(int y)const{for(size_t n=0;n<3;++n)if(y>=boxes[n].getY()&&y<boxes[n].getBottom())return n;return {};}
    static juce::String tc(int64_t frame,int rate)
    {const auto ms=juce::int64(1000.*double(frame)/std::max(1,rate));return juce::String::formatted("%02lld:%02lld:%02lld.%03lld",ms/3600000,(ms/60000)%60,(ms/1000)%60,ms%1000);}
};
}
