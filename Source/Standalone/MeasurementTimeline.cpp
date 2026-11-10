#include "MeasurementTimeline.h"
#include "TimelineGainLane.h"
#include <cmath>
#include <cstdlib>
#include <juce_cryptography/juce_cryptography.h>
#include "VideoTiming.h"
#include "Core/Measurement/GatePopulation.h"

namespace levelscope::standalone
{
juce::String diagnosticMetricText (const measurement::Metric& m, const char* unit)
{
    if (m.value) return juce::String (*m.value, 6) + " " + unit;
    if (m.status == measurement::Status::silence) return "unavailable (digital silence)";
    if (m.status == measurement::Status::belowGate) return "unavailable (no eligible blocks)";
    return "unavailable (insufficient duration)";
}
juce::String diagnosticLraInterpretation (const proposal::ProgramAnalysis& a)
{
    const bool shortProgramme = a.frames < static_cast<int64_t> (a.sampleRate) * 60;
    juce::String text = a.lra.value
        ? (a.lraStable ? "LRA stability: programme >=60 s." : "LRA stability: unstable.")
        : "LRA stability: unavailable with no eligible LRA result; unavailable is not zero.";
    if (shortProgramme) text += " Under 60 s: any available LRA may be unstable or misleading.";
    return text;
}
juce::String diagnosticWindowText (const proposal::DiagnosticWindow& w, int rate, int64_t frames, bool ungated)
{
    juce::String s;
    s << "[" << juce::String (static_cast<double> (w.begin) / rate, 6) << ", " << juce::String (static_cast<double> (w.end) / rate, 6)
      << ") s | frames [" << juce::String (w.begin) << ", " << juce::String (w.end) << ") | "
      << (w.lufs ? juce::String (*w.lufs, 6) + " LUFS" : "no finite loudness (zero window energy)")
      << (ungated ? " | ungated S" : w.included ? " | included" : " | excluded");
    if (w.padding) s << " | EOF padding; source intersection ends at " << juce::String (static_cast<double> (frames) / rate, 6) << " s";
    return s;
}
MeasurementTimeline::MeasurementTimeline()
{
    setName ("Measured audio-source timeline");
    addAndMakeVisible(blend);blend.onValueChange=[this]{repaint();};
    addAndMakeVisible(overTargetPeaks);overTargetPeaks.setName("Show peaks over target");
    overTargetPeaks.setToggleState(true,juce::dontSendNotification);overTargetPeaks.setTooltip("Sample-peak warnings use the target ceiling, or 0 dBFS with No Limit. TP! marks measured per-channel maxima.");
    overTargetPeaks.onClick=[this]{repaint();};
    for (auto* b : { &zoomIn, &zoomOut, &fit, &fitSelection,&inspect, &seek }) addAndMakeVisible (*b);
    addAndMakeVisible(follow);seek.setVisible(false);follow.setVisible(false);
    addAndMakeVisible(navigationFeedback);navigationFeedback.setFont(12.f);
    follow.setTooltip("Keep the audition playhead in view. Manual panning turns Follow off.");
    follow.onClick=[this]{if(follow.getToggleState()&&programmeFrames()>0&&(playhead<viewBegin||playhead>viewEnd))
        {const auto width=viewEnd-viewBegin;setViewRange(playhead-.2*width,playhead+.8*width);}};
    addAndMakeVisible (seconds); seconds.setName ("Audio-source inspection seconds"); seconds.setText ("0");
    addAndMakeVisible (scroll); scroll.addListener (this); scroll.setName ("Audio-source horizontal scroll");
    zoomIn.onClick = [this] { zoom (.5); }; zoomOut.onClick = [this] { zoom (2); };
    fit.onClick = [this] { if (programmeFrames()>0) setViewRange (0, static_cast<double> (programmeFrames())); };
    fitSelection.onClick=[this]{if(inspectedRange)setViewRange(inspectedRange->first,inspectedRange->second);};
    seconds.setTooltip("Audio-source seconds or an exact fraction. Inspect never moves audition; Seek is stopped-only.");
    inspect.onClick = [this]
    {
        try {if(!inspectFrame(double(enteredFrame())))throw std::invalid_argument("outside");}
        catch(const std::exception&){navigationMessage("Enter finite audio-source seconds or a fraction within the programme.");}
    };
    seek.onClick=[this]
    {
        try {requestSeek(enteredFrame(),transportIdentity);}
        catch(const std::exception&){navigationMessage("Enter valid audio-source seconds or an exact fraction for Seek.");}
    };
    seconds.onReturnKey = [this] { inspect.onClick(); };
}
void MeasurementTimeline::setAnalysis (std::shared_ptr<const proposal::ProgramAnalysis> value)
{
    if(value==analysis)return;
    analysis=std::move(value);
    if(!waveformData&&!preview&&knownFrames==0){selection.reset();inspectedRange.reset();inspectedFrame.reset();inspection.clear();viewBegin=0;viewEnd=analysis?double(analysis->frames):1;}
    rebuild();
}
void MeasurementTimeline::setPeakTarget(std::optional<double> dbtp,bool measured)
{
    if(dbtp&&!std::isfinite(*dbtp))dbtp.reset();
    if(peakTarget==dbtp&&peakMeasurementCurrent==measured)return;
    peakTarget=dbtp;peakMeasurementCurrent=measured;
    peakCeiling=dbtp?std::pow(10.,*dbtp/20.):1.;
    overTargetPeaks.setEnabled(true);
    overTargetPeaks.setTooltip("Orange circles: retained sample peaks above "+juce::String(peakTarget.value_or(0.),1)+" dBFS. No Limit uses 0 dBFS. TP! marks measured per-channel true-peak maxima, not every inter-sample over. Preview warnings are provisional; monitor gain is excluded.");
    repaint();
}
bool MeasurementTimeline::sampleExceedsTarget(double value,int64_t frame) const
{
    return peakOverlayActive()&&(waveformAccepted||isGainPreview())&&std::isfinite(value)&&frame>=viewBegin&&frame<viewEnd
        &&std::abs(value)>peakCeiling;
}
void MeasurementTimeline::setProgramme(int64_t frames,int rate,const juce::String& identity)
{
    if(identity==programmeIdentity&&frames==knownFrames&&rate==knownRate)return;
    programmeIdentity=identity;knownFrames=std::max<int64_t>(0,frames);knownRate=std::max(1,rate);
    analysis.reset();preview.reset();waveformData.reset();originalMeasurement.reset();selection.reset();inspectedRange.reset();inspectedFrame.reset();inspection.clear();seconds.setText("0",false);
    clearGainPreview();waveformAccepted=false;waveformPublication=0;
    viewBegin=0;viewEnd=std::max<double>(1,knownFrames);++navigationEpoch;rebuild();
}
void MeasurementTimeline::setWaveform(std::shared_ptr<const waveform::Overview> data,bool accepted,const juce::StringArray& channels)
{
    const auto publication=data?data->publication():0;
    const bool clearingPreview=!data&&previewBase.has_value();
    if(waveformEnabled && data==waveformData && publication==waveformPublication && accepted==waveformAccepted && channels==channelNames&&!clearingPreview)return;
    const auto previousFrames=programmeFrames();const auto previousRate=programmeRate();
    // Measurement preparation and programme replacement can paint before the
    // next timer refresh. A retired overview must also retire its gain binding.
    if(!data)clearGainPreview();
    waveformEnabled=true;waveformData=std::move(data);waveformAccepted=accepted&&bool(waveformData);waveformPublication=publication;channelNames=channels;
    if(programmeFrames()!=previousFrames || programmeRate()!=previousRate || viewEnd<=1)
    {viewBegin=0;viewEnd=std::max<double>(1,programmeFrames());selection.reset();}
    inspectedChannel=waveformData?juce::jlimit(0,waveformData->spec.channels-1,inspectedChannel):0;
    resized();
}
void MeasurementTimeline::setGainPreview(std::optional<ProcessingPlan> base,const ProcessingPlan& current)
{
    const auto key=[](const ProcessingPlan& p){const auto value=juce::String(p.staticGainDb,17)+";"+(p.sectionGains?p.sectionGains->audioKey():juce::String());return juce::SHA256(value.toRawUTF8(),value.getNumBytesAsUTF8()).toHexString();};
    const auto next=base?key(*base)+"=>"+key(current):juce::String();if(next==previewGainKey)return;
    previewGainKey=next;previewBase=std::move(base);previewCurrent=current;rebuild();
}
void MeasurementTimeline::setPreview(std::shared_ptr<const LoudnessPreview> value)
{
    const auto coverage=value?value->coverage():0;
    if(value==preview&&coverage==previewCoverage)return;
    const auto frames=programmeFrames();const auto rate=programmeRate();preview=std::move(value);previewCoverage=coverage;
    if(frames!=programmeFrames()||rate!=programmeRate()||viewEnd<=1){viewBegin=0;viewEnd=std::max<double>(1,programmeFrames());}
    rebuild();
}
juce::String MeasurementTimeline::timecode(int64_t frame,int rate)
{
    if(rate<=0||frame<0)return "--:--:--.---";
    const auto ms=video::Time(frame,rate).floorTicks(1000);
    const auto two=[](int64_t n){return juce::String(n).paddedLeft('0',2);};
    return two(ms/3600000)+":"+two(ms/60000%60)+":"+two(ms/1000%60)+"."+juce::String(ms%1000).paddedLeft('0',3);
}
void MeasurementTimeline::setTransport(int64_t frame,bool allowed,const juce::String& reason,const juce::String& identity)
{
    const bool enabled=bool(onRequestSeek);
    if(identity!=transportIdentity||allowed!=canSeek||reason!=seekReason)
    {++navigationEpoch;navigationMessage(allowed?(sectionSession?"Waveform: drag section or boundary | Alt-drag to pan":"Ruler: stopped seek | Graph: inspect"):reason);}
    transportIdentity=identity;canSeek=allowed;seekReason=reason;
    seek.setVisible(enabled);follow.setVisible(enabled);seek.setEnabled(enabled&&canSeek&&programmeFrames()>0);seek.setTooltip(reason.isEmpty()?"Move audition to these seconds while stopped; no playback.":reason);
    if(playhead==frame)return;playhead=frame;
    if(follow.getToggleState()&&programmeFrames()>0&&(frame<viewBegin||frame>viewEnd))
        setViewRange(double(frame)-(viewEnd-viewBegin)*.2,double(frame)+(viewEnd-viewBegin)*.8);
    repaint();
}
void MeasurementTimeline::navigationMessage(const juce::String& text)
{navigationFeedback.setText(text,juce::dontSendNotification);navigationFeedback.setTooltip(text);}
void MeasurementTimeline::requestSeek(int64_t frame,const juce::String& identity)
{
    if(!canSeek){navigationMessage(seekReason.isEmpty()?"Stopped seeking is unavailable.":seekReason);return;}
    if(frame<0||frame>programmeFrames()){navigationMessage("Seek is outside the programme.");return;}
    if(onRequestSeek)
    {
        focusForTransport();
        const auto error=onRequestSeek(frame,identity);
        navigationMessage(error.isEmpty()?"Audition positioned at "+timecode(frame,programmeRate())+" | frame "+juce::String(frame)+". Stopped; inspection is independent.":error);
    }
}
void MeasurementTimeline::focusForTransport()
{
    // Seek is disabled as soon as playback starts. Hold focus on a stable
    // surface instead of letting JUCE choose a neighbouring text editor.
    setWantsKeyboardFocus(true);
    if(isShowing()&&getPeer()!=nullptr)grabKeyboardFocus();
}
bool MeasurementTimeline::keyPressed(const juce::KeyPress& key)
{
    if(juce::Component::getCurrentlyModalComponent())return false;
    for(auto* c=juce::Component::getCurrentlyFocusedComponent();c;c=c->getParentComponent())
        if(dynamic_cast<juce::TextEditor*>(c))return false;
    if((key==juce::KeyPress::backspaceKey||key==juce::KeyPress::deleteKey)&&!key.getModifiers().isAnyModifierKeyDown()&&gainLane
        &&juce::Component::getCurrentlyFocusedComponent()==this)
        return gainLane->deleteSelected();
    if(key.getModifiers().isCommandDown())
    {
        if((key.getKeyCode()=='Z'||key.getKeyCode()=='z'))
        {
            const bool redo=key.getModifiers().isShiftDown();
            if(onUndo)return onUndo(redo);
            if(sectionSession&&sectionSession->undoManual(redo)){refreshSectionEditing();if(sectionEdited)sectionEdited(0);return true;}
            return false;
        }
        if((key.getKeyCode()=='R'||key.getKeyCode()=='r')&&onRemeasure){onRemeasure();return true;}
    }
    if(!key.getModifiers().isAnyModifierKeyDown()
        &&(key.getKeyCode()=='A'||key.getKeyCode()=='a')&&onToggleInspectionAudition)
    {onToggleInspectionAudition();return true;}
    if(key!=juce::KeyPress::spaceKey||!onTogglePlayback)return false;
    focusForTransport();onTogglePlayback();return true;
}
void MeasurementTimeline::attachGainLane(TimelineGainLane& lane)
{
    gainLane=&lane;addAndMakeVisible(lane);
    lane.onSelection=[this](int64_t begin,int64_t end,int64_t frame)
    {
        inspectedFrame=frame;inspectedEntry=juce::String(double(frame)/programmeRate(),9);seconds.setText(inspectedEntry,false);
        highlight(double(begin),double(end),false);
        if(frame<viewBegin||frame>viewEnd){const auto width=viewEnd-viewBegin;setViewRange(frame-width*.2,frame+width*.8);}
        if(onInspect)onInspect("Current section gain settings  |  "+gainLane->inspectionText()+"\nInspected audio: "+rangePeakText(begin,end));
    };
    resized();
}
int64_t MeasurementTimeline::frameAt(float x) const
{
    if(graph.getWidth()<=0)return 0;
    return std::clamp<int64_t>(int64_t(std::floor(viewBegin+std::clamp((double(x)-graph.getX())/std::max(1,graph.getWidth()-1),0.,1.)*(viewEnd-viewBegin))),0,programmeFrames());
}
int64_t MeasurementTimeline::enteredFrame() const
{
    const auto entry=seconds.getText().trim();
    // Inspect displays decimal seconds. Retain the exact cursor instead of
    // quantizing its rounded presentation again when the entry hasn't changed.
    return inspectedFrame&&entry==inspectedEntry?*inspectedFrame:video::Time::parse(entry.toStdString()).floorTicks(programmeRate());
}
void MeasurementTimeline::setViewRange (double begin, double end)
{
    const auto frames=double(programmeFrames());
    if(frames<=0 || !std::isfinite(begin) || !std::isfinite(end) || end<=begin)return;
    const auto finest=waveformData?double(waveformData->binFrames):programmeRate()*.1;
    const auto width=std::clamp(end-begin,std::min(frames,finest),frames);
    viewBegin=std::clamp(begin,0.,frames-width);viewEnd=viewBegin+width;++navigationEpoch;rebuild();
}
void MeasurementTimeline::zoom (double factor)
{
    const auto centre = selection?std::clamp(selection->first,viewBegin,viewEnd):(viewBegin + viewEnd) / 2, half = (viewEnd - viewBegin) * factor / 2;
    setViewRange (centre - half, centre + half);
}
void MeasurementTimeline::zoomAt(double factor,float px)
{
    if(graph.getWidth()<2||programmeFrames()<=0)return;
    const auto fraction=std::clamp((double(px)-graph.getX())/(graph.getWidth()-1),0.,1.);
    const auto anchor=viewBegin+fraction*(viewEnd-viewBegin);
    const auto finest=waveformData?double(waveformData->binFrames):programmeRate()*.1;
    const auto width=std::clamp((viewEnd-viewBegin)*factor,std::min(double(programmeFrames()),finest),double(programmeFrames()));
    setViewRange(anchor-fraction*width,anchor+(1-fraction)*width);
}
void MeasurementTimeline::setReferenceMeasurement(std::shared_ptr<const measurement::Result> value,bool current)
{
    if(value==originalMeasurement&&current==currentInspection)return;
    if(value!=originalMeasurement)
        originalIntegratedGateEnergy=value&&value->integratedRelativeGate
            ?measurement::relativeGateEnergy(value->momentary,false):std::nullopt;
    originalMeasurement=std::move(value);currentInspection=current;rebuild();
}
juce::String MeasurementTimeline::rangePeakText(int64_t begin,int64_t end) const
{
    if(previewBase)return "Current gain preview only. Peak values and loudness need re-measurement.";
    juce::String text="Selected frames ["+juce::String(begin)+", "+juce::String(end)+")\n";
    bool foundTrue=false;
    if(waveformData)
    {
        const auto bins=waveformData->reduce(begin,end,1);
        if(!bins.empty()&&bins[0].available())
        {
            const auto& bin=bins[0];
            text+="PCM overview covers ["+juce::String(bin.begin)+", "+juce::String(bin.coveredEnd)+") frames";
            if(!waveformAccepted||bin.coveredEnd<end)text+=" — INCOMPLETE";
            text+=".\n";
            for(int c=0;c<waveformData->spec.channels;++c)
            {
                const auto& e=bin.channels[size_t(c)];
                const bool negative=std::abs(e.minimum)>std::abs(e.maximum)||(std::abs(e.minimum)==std::abs(e.maximum)&&e.minimumFrame<e.maximumFrame);
                const auto magnitude=std::abs(negative?e.minimum:e.maximum);const auto frame=negative?e.minimumFrame:e.maximumFrame;
                const bool exact=waveformAccepted&&bin.coveredEnd>=end&&frame>=begin&&frame<end;
                const auto name=channelNames[c].isEmpty()?juce::String(c+1):channelNames[c];
                text+=name+": "+(exact?"selected sample peak ":"covering overview peak ")+(magnitude>0?juce::String(20*std::log10(magnitude),3)+" dBFS":"-inf dBFS (silence)")+" at frame "+juce::String(frame);
                if(!exact)text+="; boundary bins/coverage prevent an exact selected-interval maximum";
                text+="\n";
            }
        }
        else text+="PCM overview coverage pending.\n";
    }
    else text+="Sample peak unavailable for this range: Build waveform for the inspected state.\n";
    if(analysis&&!preview)for(const auto& peak:analysis->peaks)
        if(peak.trueDbtp&&peak.truePeakEighthFrame>=begin*8&&peak.truePeakEighthFrame<end*8)
        {foundTrue=true;text+="Channel "+juce::String(peak.channel+1)+" measured programme true-peak maximum lies here: "+juce::String(*peak.trueDbtp,3)+" dBTP at "+juce::String(peak.truePeakEighthFrame)+"/8 frames.\n";}
    if(!foundTrue)text+="Regional true peak is not available from the PCM overview; no peak-safety claim.\n";
    return text;
}
void MeasurementTimeline::highlight (double begin, double end, bool reveal)
{
    if (programmeFrames()<=0) return;
    selection = std::make_pair (begin, end);
    inspectedRange=std::make_pair(begin,std::max(begin+1,end));fitSelection.setEnabled(true);
    const auto duration = static_cast<double> (programmeFrames());
    if (reveal && (begin < viewBegin || end > viewEnd))
    {
        const auto width = std::max (viewEnd - viewBegin, std::min (duration, programmeRate() * 5.));
        setViewRange (std::clamp ((begin + end) / 2, 0., duration) - width / 2,
                      std::clamp ((begin + end) / 2, 0., duration) + width / 2);
    }
    repaint();
}
bool MeasurementTimeline::inspectSeconds (double value)
{return inspectFrame(value*programmeRate());}
bool MeasurementTimeline::inspectFrame(double frame)
{
    const auto rate=programmeRate();const auto frames=programmeFrames();
    if(frames<=0||!std::isfinite(frame)||frame<0||frame>frames)return false;
    const auto value=frame/rate;
    if(gainLane)gainLane->selectFrame(int64_t(frame),false);
    inspectedFrame=int64_t(std::floor(frame));inspectedEntry=juce::String(double(*inspectedFrame)/rate,9);seconds.setText(inspectedEntry,false);
    inspectedRange=std::make_pair(std::min(double(frames-1),frame),std::min(double(frames),frame+1));
    inspection="Audio-source "+juce::String(value,6)+" s; cursor frame "+juce::String(int64_t(std::floor(frame)))+"\n";
    if(waveformData&&!previewBase)
    {
        const auto bin=waveformData->inspect(std::min<int64_t>(frames-1,int64_t(std::floor(frame))));
        if(bin)
        {
            inspectedRange=std::make_pair(double(bin->begin),double(bin->end));
            const auto name=channelNames[inspectedChannel].isEmpty()?juce::String(inspectedChannel+1):channelNames[inspectedChannel];
            inspection<<"PCM overview | channel "<<name<<" | represented frames ["<<juce::String(bin->begin)<<", "<<juce::String(bin->end)
                <<") | "<<juce::String(bin->end-bin->begin)<<" frames / "<<juce::String(1000.*(bin->end-bin->begin)/rate,6)<<" ms\n";
            if(bin->available())
            {
                const auto& e=bin->channels[size_t(inspectedChannel)];
                inspection<<"Minimum "<<juce::String(e.minimum,9)<<" at frame "<<juce::String(e.minimumFrame)
                    <<"; maximum "<<juce::String(e.maximum,9)<<" at frame "<<juce::String(e.maximumFrame)<<"\n";
            }
            else inspection<<"Pending coverage\n";
            if(!waveformAccepted)inspection<<"ANALYSING - incomplete waveform overview\n";
        }
        selection=std::make_pair(frame,frame);repaint();
    }
    if(preview)
    {
        inspection<<"ANALYSING - incomplete, not a programme result. M gate membership pending.\n";
        for(bool s:{false,true})
        {
            inspection<<(s?"S: ":"M: ");const auto w=preview->inspect(s,frame,previewCoverage);
            if(w&&!s&&!waveformData)inspectedRange=std::make_pair(double(w->endFrame-rate*4/10),double(w->endFrame));
            if(w)inspection<<"frames ["<<juce::String(w->endFrame-(s?rate*3:rate*4/10))<<", "<<juce::String(w->endFrame)<<") | "
                <<(w->energy>0?juce::String(measurement::OfflineMeasurement::loudness(w->energy),6)+" LUFS":"no finite loudness (zero window energy)")<<(s?" | ungated S\n":" | gate pending\n");
            else inspection<<"Pending coverage or no complete window yet\n";
        }
    }
    else if(analysis)
    {
        inspection<<"Measured loudness: nearest sampled window ends (earlier wins ties).\n";
        for(bool shortTerm:{false,true})
        {
            const auto& windows=shortTerm?analysis->shortTerm:analysis->momentary;const auto n=proposal::inspectWindow(windows,frame);
            inspection<<(shortTerm?"S: ":"M: ");
            if(n)
            {
                inspection<<diagnosticWindowText(windows[*n],analysis->sampleRate,analysis->frames,shortTerm)<<"\n";
                if(!waveformData && !shortTerm)highlight(double(windows[*n].begin),double(windows[*n].end));
            }
            else inspection<<"no complete windows\n";
        }
    }
    else inspection<<"Momentary and short-term measurements pending.\n";
    fitSelection.setEnabled(true);if(onInspect)onInspect(inspection);return true;
}
void MeasurementTimeline::scrollBarMoved (juce::ScrollBar*, double start)
{ if (! updatingScroll) {follow.setToggleState(false,juce::dontSendNotification);setViewRange (start, start + viewEnd - viewBegin);} }
void MeasurementTimeline::rebuild()
{
    if(gainLane)gainLane->setViewport(viewBegin,viewEnd,graph.getX(),graph.getWidth());
    plot = {};
    if(preview){previewCoverage=preview->coverage();plot=preview->reduce(viewBegin,viewEnd,std::max(1,graph.getWidth()),previewCoverage);}
    else if (analysis) plot = proposal::prepareTimeline (*analysis, viewBegin, viewEnd, std::max (1, graph.getWidth()));
    // UI-only classification. Short-term loudness remains the original ungated
    // series; its colour describes membership in the independently measured LRA population.
    gatedShortTerm.clear();changedMomentary.clear();
    if(analysis&&!preview&&viewEnd>viewBegin)
    {
        const auto count=std::max<size_t>(1,plot.shortTerm.size());gatedShortTerm.resize(count);
        const auto column=[&](int64_t end){return size_t(std::clamp(int((end-viewBegin)*count/(viewEnd-viewBegin)),0,int(count)-1));};
        auto w=std::lower_bound(analysis->lraWindows.begin(),analysis->lraWindows.end(),viewBegin,[](const auto& v,double f){return v.end<f;});
        for(;w!=analysis->lraWindows.end()&&w->end<=viewEnd&&w->end<=analysis->frames;++w)
            if(w->lufs)(w->included?gatedShortTerm[column(w->end)].included:gatedShortTerm[column(w->end)].excluded).add(*w->lufs,w->end);
        if(currentInspection&&originalMeasurement&&originalIntegratedGateEnergy&&analysis->integratedGate
            &&originalMeasurement->frames==analysis->frames&&originalMeasurement->sampleRate==analysis->sampleRate)
        {
            changedMomentary.resize(count,false);
            const auto& source=originalMeasurement->momentary;
            auto m=std::lower_bound(analysis->momentary.begin(),analysis->momentary.end(),viewBegin,[](const auto& v,double f){return v.end<f;});
            for(;m!=analysis->momentary.end()&&m->end-analysis->sampleRate/10<viewEnd;++m)
            {
                const auto n=size_t(m-analysis->momentary.begin());
                // Compare the exact measured IL populations, each against its own
                // gate. A gain/value difference alone does not change membership.
                if(m->end>viewBegin&&n<source.size()&&source[n].endFrame==m->end
                    &&measurement::includedInPopulation(source[n].energy,*originalIntegratedGateEnergy,false)!=m->included)
                {
                    const auto first=column(std::max<int64_t>(0,m->end-analysis->sampleRate/10)),last=column(m->end);
                    for(auto c=first;c<=last;++c)changedMomentary[c]=true;
                }
            }
        }
    }
    if(waveformData)waveformPublication=waveformData->publication();
    waveformColumns=waveformData?waveformData->reduce(int64_t(std::floor(viewBegin)),int64_t(std::ceil(viewEnd)),std::max(1,graph.getWidth()),waveformPublication):std::vector<waveform::Column>{};
    if(isGainPreview())
    {
        const auto gain=[](const ProcessingPlan& p,int64_t f){return p.staticGainDb+(p.sectionGains?p.sectionGains->gainAt(f):0.);};
        for(auto& col:waveformColumns)if(col.available())for(int c=0;c<waveformData->spec.channels;++c)
        {
            auto& e=col.channels[size_t(c)];
            e.minimum=float(e.minimum*std::pow(10.,(gain(previewCurrent,e.minimumFrame)-gain(*previewBase,e.minimumFrame))/20.));
            e.maximum=float(e.maximum*std::pow(10.,(gain(previewCurrent,e.maximumFrame)-gain(*previewBase,e.maximumFrame))/20.));
            if(e.minimum>e.maximum){std::swap(e.minimum,e.maximum);std::swap(e.minimumFrame,e.maximumFrame);}
        }
    }
    updatingScroll = true;
    scroll.setRangeLimits (0, std::max<double>(1,programmeFrames()));
    scroll.setCurrentRange (viewBegin, viewEnd - viewBegin, juce::dontSendNotification);
    updatingScroll = false;
    for (auto* b : { &zoomIn, &zoomOut, &fit, &inspect }) b->setEnabled (programmeFrames()>0);
    fitSelection.setEnabled(programmeFrames()>0&&inspectedRange.has_value());
    seconds.setEnabled (programmeFrames()>0); scroll.setEnabled (programmeFrames()>0);
    repaint();
}
void MeasurementTimeline::resized()
{
    const auto previousGraph=graph;
    auto b = getLocalBounds(); auto row = b.removeFromTop (28);
    zoomIn.setBounds (row.removeFromLeft (80)); zoomOut.setBounds (row.removeFromLeft (80));
    fit.setBounds(row.removeFromLeft(125));fitSelection.setBounds(row.removeFromLeft(100));row.removeFromLeft(8);seconds.setBounds(row.removeFromLeft(110));inspect.setBounds(row.removeFromLeft(118));seek.setBounds(row.removeFromLeft(65));follow.setBounds(row.removeFromLeft(76));
    navigationFeedback.setBounds(row);
    overTargetPeaks.setBounds(std::max(10,getWidth()-185),30,180,20);
    scroll.setBounds (b.removeFromBottom (16).withTrimmedLeft (55));
    if(gainLane)gainLane->setBounds(b.removeFromBottom(gainLane->preferredHeight()));
    b.removeFromTop (44); ruler=b.removeFromTop(waveformEnabled?24:0).withTrimmedLeft(waveformEnabled?85:55).withTrimmedRight(waveformEnabled?60:12); b.removeFromBottom (25); graph = b.withTrimmedLeft(waveformEnabled?85:55).withTrimmedRight(waveformEnabled?60:12);
    blend.setBounds(2,graph.getY(),18,std::max(0,graph.getHeight()));
    if(graph!=previousGraph)++navigationEpoch;
    rebuild();
}
void MeasurementTimeline::bindSectionEditing(StandaloneSessionModel& session,std::function<void(uint64_t)> changed)
{
    if(sectionSession!=&session)sectionEdit=std::make_unique<SectionEditPrompt>(session);
    sectionSession=&session;
    sectionEdited=std::move(changed);
    refreshSectionEditing();
}
void MeasurementTimeline::refreshSectionEditing()
{
    if(!sectionSession)return;
    auto& session=*sectionSession;
    const auto& map=session.sectionDocument.map();
    const bool editable=knownFrames>0 && map.frames==knownFrames && !map.sections.empty() && session.sourceDocument.layout.confirmed
        && session.manualEditAllowed();
    if(sectionRuntime==session.sessionId && audioRevision==session.revision && sectionInput==map.inputIdentity
        && sectionRevision==map.revision && sectionEditable==editable)return;
    sectionRuntime=session.sessionId;audioRevision=session.revision;sectionInput=map.inputIdentity;
    sectionRevision=map.revision;sectionEditable=editable;
    sectionBoundaries.clear();
    for(const auto& row:map.sections)sectionBoundaries.push_back(row.end);
    if(sectionGesture)gesture=false;
    sectionGesture=false;
    ++navigationEpoch;
    repaint();
}
std::optional<size_t> MeasurementTimeline::boundaryHandle(float px) const
{
    if(!sectionEditable || sectionBoundaries.size()<2 || graph.getWidth()<2)return {};
    const auto last=sectionBoundaries.end()-1; // Programme start and end are fixed.
    const auto near=std::lower_bound(sectionBoundaries.begin(),last,frameAt(px));
    std::optional<size_t> best;
    double distance=5.0;
    const auto consider=[&](auto it)
    {
        if(it==last || *it<viewBegin || *it>viewEnd)return;
        const auto x=graph.getX()+(*it-viewBegin)*(graph.getWidth()-1)/(viewEnd-viewBegin);
        if(std::abs(px-x)<=distance){distance=std::abs(px-x);best=size_t(it-sectionBoundaries.begin());}
    };
    if(near!=sectionBoundaries.begin())consider(near-1);
    consider(near);
    return best;
}
juce::String MeasurementTimeline::peakHoverText(juce::Point<float> pointer) const
{
    if(!peakOverlayActive()||blend.waveformAlpha()<=0||!getLocalBounds().toFloat().contains(pointer)||viewEnd<=viewBegin)return {};
    const auto x=[&](double frame){return float(graph.getX()+(frame-viewBegin)*(graph.getWidth()-1)/(viewEnd-viewBegin));};
    float distance=8.f;juce::String text;
    const auto describe=[&](int channel,double db,double frame,bool truePeak)
    {
        const auto threshold=truePeak&&!currentInspection?std::min(peakTarget.value_or(0.),0.):peakTarget.value_or(0.);
        const auto name=channelNames[channel].isNotEmpty()?channelNames[channel]:"Channel "+juce::String(channel+1);
        return juce::String(currentInspection?"Current":"Original")+" | "+name+" | "+timecode(int64_t(frame),programmeRate())
            +"\n"+(truePeak?"Measured programme true peak: ":previewBase?"Preview sample peak: ":"Sample peak: ")+juce::String(db,3)+(truePeak?" dBTP":" dBFS")
            +"\n"+(peakTarget?"Target: "+juce::String(*peakTarget,2)+" dBTP":juce::String(truePeak?"Warning threshold: 0 dBTP (No Limit)":"Warning threshold: 0 dBFS (No Limit)"))
            +"\nApproximate additional reduction: "+juce::String(threshold-db,3)+" dB"
            +"\n"+(truePeak?"Guidance for this measured maximum. Re-measure Current after editing.":"Sample-based guidance; inter-sample peaks may be higher. Re-measure Current after editing.");
    };
    if(waveformData&&(waveformAccepted||previewBase)&&graph.toFloat().contains(pointer))
    {
        const int channels=waveformData->spec.channels;const float height=float(graph.getHeight())/channels;
        double scale=1;for(const auto& col:waveformColumns)if(col.available())for(int c=0;c<channels;++c)
            scale=std::max(scale,double(std::max(std::abs(col.channels[size_t(c)].minimum),std::abs(col.channels[size_t(c)].maximum))));
        const int channel=std::clamp(int((pointer.y-graph.getY())/height),0,channels-1);
        const float centre=graph.getY()+(channel+.5f)*height;
        const auto visit=[&](float value,int64_t frame)
        {
            if(!sampleExceedsTarget(value,frame))return;
            const auto d=pointer.getDistanceFrom({x(double(frame)),centre-float(value/scale)*height*.44f});
            if(d<distance){distance=d;text=describe(channel,20*std::log10(std::abs(value)),double(frame),false);}
        };
        for(const auto& col:waveformColumns)if(col.available()){const auto& e=col.channels[size_t(channel)];visit(e.minimum,e.minimumFrame);visit(e.maximum,e.maximumFrame);}
    }
    if(text.isNotEmpty())return text;
    if(analysis&&peakMeasurementCurrent&&!previewBase)for(const auto& peak:analysis->peaks)if(peak.trueDbtp&&*peak.trueDbtp>std::min(peakTarget.value_or(0.),currentInspection?peakTarget.value_or(0.):0.))
    {
        const double frame=peak.truePeakEighthFrame/8.;if(frame<viewBegin||frame>viewEnd)continue;
        const auto top=waveformData?graph.getY()+peak.channel*graph.getHeight()/waveformData->spec.channels:graph.getY();
        if(juce::Rectangle<float>(x(frame)-3,float(top),51,15).contains(pointer))return describe(peak.channel,*peak.trueDbtp,frame,true);
    }
    return {};
}
void MeasurementTimeline::mouseMove(const juce::MouseEvent& e)
{
    setMouseCursor(graph.contains(e.getPosition()) && boundaryHandle(e.position.x)
        ? juce::MouseCursor::LeftRightResizeCursor : juce::MouseCursor::NormalCursor);
}
void MeasurementTimeline::changeSectionEnd(uint64_t id,int64_t frame)
{
    if(!sectionSession||!sectionEdit)return;
    const auto before=sectionSession->sectionDocument.map();auto edited=sectionSession->sectionDocument;
    const auto it=std::find_if(before.sections.begin(),before.sections.end(),[&](const auto& r){return r.id==id;});
    if(it==before.sections.end()||it->end==frame)return;
    juce::String error;if(!edited.moveEnd(id,frame,error)){navigationMessage(error);return;}
    sectionEdit->request(*this,before,edited.map(),false,it->begin,it->end,it->end,frame,
        [safe=juce::Component::SafePointer<MeasurementTimeline>(this),id](bool accepted,bool)
        {
            if(!safe)return;auto& self=*safe;self.refreshSectionEditing();
            if(self.sectionEdited)self.sectionEdited(id);
            self.navigationMessage(accepted?("Section changed. Command-Z undoes it."):self.sectionSession->sectionTrimNotice);
        });
}
void MeasurementTimeline::finishSectionGesture(int64_t end)
{
    refreshSectionEditing();
    if(!sectionEditable || gestureEpoch!=navigationEpoch || !sectionSession->sourceDocument.matchesFileStamps())return;
    auto document=sectionSession->sectionDocument;
    const auto& map=document.map();
    juce::String error;
    uint64_t selected=0;
    const auto before=map.revision;
    bool accepted=false;
    if(movingBoundary)
    {
        const auto n=*movingBoundary;
        changeSectionEnd(map.sections[n].id,std::clamp(end,map.sections[n].begin+1,map.sections[n+1].end-1));
        return;
    }
    else
    {
        const auto start=frameAt(downPoint.x);
        if(start==end)return;
        accepted=document.defineInterval(std::min(start,end),std::max(start,end),selected,error);
    }
    if(!accepted){navigationMessage(error);return;}
    if(map.revision!=before){auto next=sectionSession->manualState();next.map=map;if(!sectionSession->commitManualEdit(std::move(next),"Draw section")){navigationMessage(sectionSession->sectionTrimNotice);return;}}
    refreshSectionEditing();
    if(sectionEdited)sectionEdited(selected);
    const auto it=std::find_if(map.sections.begin(),map.sections.end(),[selected](const auto& row){return row.id==selected;});
    if(it!=map.sections.end())highlight(double(it->begin),double(it->end),false);
    navigationMessage("Section saved. Drag a boundary to adjust it.");
}
void MeasurementTimeline::mouseDown(const juce::MouseEvent& e)
{
    if(gainLane&&!gainLane->commitPending())return;
    refreshSectionEditing();
    gestureGain=gainLane&&gainLane->inspectionBounds().contains(e.getPosition()-gainLane->getPosition());
    if(gainLane&&!gestureGain)gainLane->clearDeleteSelection();
    selectGesture=gestureGain||(e.mods.isRightButtonDown()&&graph.contains(e.getPosition()));
    gesture=(e.mods.isLeftButtonDown()||selectGesture)&&programmeFrames()>0&&(graph.contains(e.getPosition())||ruler.contains(e.getPosition())||gestureGain);
    sectionGesture=gesture && !selectGesture && sectionEditable && graph.contains(e.getPosition()) && !e.mods.isAltDown();
    movingBoundary=sectionGesture?boundaryHandle(e.position.x):std::nullopt;
    dragFrame=frameAt(e.position.x);
    dragging=false;gestureRuler=ruler.contains(e.getPosition());downPoint=e.position;dragBegin=viewBegin;dragEnd=viewEnd;gestureEpoch=navigationEpoch;gestureIdentity=transportIdentity;
}
void MeasurementTimeline::mouseDrag(const juce::MouseEvent& e)
{
    if(!gesture)return;
    if(selectGesture)
    {
        if(gestureEpoch!=navigationEpoch||!gainLane)return;
        if(e.position.getDistanceFrom(downPoint)>=4){dragging=true;gainLane->selectRange(frameAt(downPoint.x),frameAt(e.position.x));}
        return;
    }
    if(sectionGesture)
    {
        if(gestureEpoch!=navigationEpoch)return;
        dragFrame=frameAt(e.position.x);
        if(movingBoundary || std::abs(e.position.x-downPoint.x)>=4)dragging=true;
        if(movingBoundary)
        {
            const auto& rows=sectionSession->sectionDocument.map().sections;
            if(*movingBoundary+1<rows.size())dragFrame=std::clamp(dragFrame,rows[*movingBoundary].begin+1,rows[*movingBoundary+1].end-1);
        }
        repaint();
        return;
    }
    if(e.position.getDistanceFrom(downPoint)<4)return;
    dragging=true;follow.setToggleState(false,juce::dontSendNotification);
    const auto shift=-(e.position.x-downPoint.x)/std::max(1,graph.getWidth())*(dragEnd-dragBegin);
    setViewRange(dragBegin+shift,dragEnd+shift);
}
void MeasurementTimeline::mouseUp(const juce::MouseEvent& e)
{
    if(!gesture)return;
    gesture=false;
    if(sectionGesture)
    {
        sectionGesture=false;
        repaint();
        const auto delta=std::abs(e.position.x-downPoint.x);
        if((movingBoundary && delta>0) || delta>=4){finishSectionGesture(frameAt(e.position.x));return;}
    }
    if(selectGesture&&dragging){focusForTransport();return;}
    if(dragging||e.position.getDistanceFrom(downPoint)>=4||gestureEpoch!=navigationEpoch)return;
    if(gestureRuler&&ruler.contains(e.getPosition()))requestSeek(frameAt(e.position.x),gestureIdentity);
    else if(gestureGain&&gainLane&&gainLane->inspectionBounds().contains(e.getPosition()-gainLane->getPosition()))
    {focusForTransport();const auto y=e.y-gainLane->getY();
        if(y>=gainLane->plotBounds().getBottom())gainLane->editReadout(frameAt(e.position.x));
        else gainLane->selectAt(frameAt(e.position.x),y);}
    else if(!gestureRuler&&graph.contains(e.getPosition()))
    {
        if(waveformData)inspectedChannel=juce::jlimit(0,waveformData->spec.channels-1,(e.y-graph.getY())*waveformData->spec.channels/std::max(1,graph.getHeight()));
        inspectFrame(double(frameAt(e.position.x)));
    }
}
void MeasurementTimeline::mouseWheelMove(const juce::MouseEvent& e,const juce::MouseWheelDetails& w)
{
    follow.setToggleState(false,juce::dontSendNotification);
    if((e.mods.isCommandDown()||e.mods.isCtrlDown())&&!e.mods.isShiftDown())zoomAt(w.deltaY>0?.8:1.25,e.position.x);
    else {const auto shift=-(w.deltaX!=0?w.deltaX:w.deltaY)*(viewEnd-viewBegin);setViewRange(viewBegin+shift,viewEnd+shift);}
}
void MeasurementTimeline::paint (juce::Graphics& g)
{
    g.fillAll (palette::graph()); g.setFont (12.0f); g.setColour (juce::Colours::white);
    if(waveformEnabled)
    {
        juce::String heading=waveformData?"PCM waveform overview | "+juce::String(waveformData->binFrames)+" frames/bin ("+juce::String(1000.*waveformData->binFrames/programmeRate(),3)+" ms)":
            "PCM waveform unavailable | use Build waveform for an existing measurement";
        if(waveformData && !waveformAccepted&&!previewBase)heading+=" | ANALYSING - incomplete: "+juce::String(waveformData->coveredFrames(waveformPublication))+" / "+juce::String(programmeFrames())+" frames";
        if(isGainPreview())heading="CURRENT GAIN PREVIEW | "+juce::String(1000.*waveformData->binFrames/programmeRate(),3)+" ms/bin | Not measured";
        if(overTargetPeaks.getToggleState())
        {
            bool overs=false;for(const auto& col:waveformColumns)if(col.available())for(int c=0;waveformData&&c<waveformData->spec.channels;++c)
                overs=overs||sampleExceedsTarget(col.channels[size_t(c)].minimum,col.channels[size_t(c)].minimumFrame)||sampleExceedsTarget(col.channels[size_t(c)].maximum,col.channels[size_t(c)].maximumFrame);
            if(analysis&&peakMeasurementCurrent)for(const auto& p:analysis->peaks)overs=overs||(p.trueDbtp&&*p.trueDbtp>std::min(peakTarget.value_or(0.),currentInspection?peakTarget.value_or(0.):0.)&&p.truePeakEighthFrame/8.>=viewBegin&&p.truePeakEighthFrame/8.<viewEnd);
            heading+=!(waveformAccepted||isGainPreview())?" | Peak data unavailable":overs?" | Peak warnings":" | "+juce::String(isGainPreview()?"Preview: ":"")+"No overs in this view.";
        }
        g.drawText(heading,10,31,std::max(0,getWidth()-205),18,juce::Justification::left);
        g.drawText(preview?"ANALYSING - incomplete | M amber, S blue: gates pending | LUFS at right; PCM at left.":"M amber / excluded grey; S blue / outside LRA gate blue-grey. Current: red IL gate = M gate membership differs from Original.",10,49,getWidth()-20,18,juce::Justification::left);
    }
    else
    {
        g.drawText("Audio-source time | M included: orange; M excluded: grey; S: blue (ungated)",55,31,std::max(0,getWidth()-250),18,juce::Justification::left);
        g.drawText("Points at window ends, every 100 ms. Per-column min/max retains extrema and mixed membership; no interpolation across gaps.",55,49,getWidth()-60,18,juce::Justification::left);
    }
    if(!ruler.isEmpty())
    {
        g.setColour(juce::Colour(0xff28323e));g.fillRect(ruler);
        g.setColour(juce::Colours::lightgrey);g.drawText("Timecode",0,ruler.getY(),80,ruler.getHeight(),juce::Justification::right);
        for(int n=0;n<=4;++n)
        {
            const auto frame=int64_t(std::floor(viewBegin+(viewEnd-viewBegin)*n/4));
            const int px=ruler.getX()+(ruler.getWidth()-1)*n/4;
            g.drawVerticalLine(px,float(ruler.getY()),float(ruler.getY()+4));
            g.drawText(timecode(frame,programmeRate()),std::clamp(px-(n==4?100:0),ruler.getX(),std::max(ruler.getX(),ruler.getRight()-100)),ruler.getY()+4,100,20,juce::Justification::left);
        }
        if(onRequestSeek&&playhead>=viewBegin&&playhead<=viewEnd)
        {
            const auto px=float(ruler.getX()+(playhead-viewBegin)*(ruler.getWidth()-1)/(viewEnd-viewBegin));
            const auto bottom=float(ruler.getBottom());juce::Path head;head.addTriangle(px-4,bottom-5,px+4,bottom-5,px,bottom);g.setColour(juce::Colours::white);g.fillPath(head);
        }
    }
    if(programmeFrames()==0||graph.isEmpty()){g.drawText("No measured timeline",getLocalBounds(),juce::Justification::centred);return;}
    if(waveformData){g.beginTransparencyLayer(blend.waveformAlpha());paintWaveform(g);g.endTransparencyLayer();}
    double low = -75, high = -10;
    for (const auto* series : { &plot.momentary, &plot.shortTerm }) for (const auto& b : *series)
        for (const auto* v : { &b.included, &b.excluded, &b.ungated }) if (v->count)
        { low = std::min (low, v->minimum - 2); high = std::max (high, v->maximum + 2); }
    const auto y = [&] (double value) { return static_cast<float> (graph.getBottom() - (value - low) * graph.getHeight() / (high - low)); };
    const auto x = [&] (double frame) { return static_cast<float> (graph.getX() + (frame - viewBegin) * (graph.getWidth()-1) / (viewEnd - viewBegin)); };
    g.beginTransparencyLayer(blend.loudnessAlpha());
    for (double v = std::ceil (low / 10) * 10; v <= high; v += 10)
    {
        g.setColour (juce::Colour (0xff303744).withAlpha(waveformData?.3f:1.f)); g.drawHorizontalLine (static_cast<int> (y(v)), static_cast<float> (graph.getX()), static_cast<float> (graph.getRight()));
        g.setColour (juce::Colours::lightgrey); g.drawText(juce::String(v,0)+(waveformEnabled?" LUFS":""),waveformEnabled?graph.getRight()+3:0,static_cast<int>(y(v))-7,waveformEnabled?57:48,15,juce::Justification::right);
    }
    g.endTransparencyLayer();g.saveState(); g.reduceClipRegion (graph);
    // Intervals outside the viewport must not become fabricated edge markers.
    // Zero-width channel-peak markers are visible only at their actual position.
    if (selection && ((selection->second > viewBegin && selection->first < viewEnd)
        || (selection->first == selection->second && selection->first >= viewBegin && selection->first <= viewEnd)))
    {
        const auto left = x (std::clamp (selection->first, viewBegin, viewEnd));
        const auto right = x (std::clamp (selection->second, viewBegin, viewEnd));
        g.setColour (juce::Colours::yellow.withAlpha (.16f)); g.fillRect (left, static_cast<float> (graph.getY()), std::max (2.f, right - left), static_cast<float> (graph.getHeight()));
    }
    g.beginTransparencyLayer(blend.loudnessAlpha());
    for (const auto* series : { &plot.momentary, gatedShortTerm.empty()?&plot.shortTerm:&gatedShortTerm })
    {
        std::array<const proposal::ValueExtent*,3> previous{};
        for (size_t n = 0; n < series->size(); ++n)
        {
            const auto px = static_cast<float> (graph.getX()) + (static_cast<float> (n) + .5f) * graph.getWidth() / static_cast<float> (series->size());
            const auto& b = (*series)[n]; int category = 0;
            for (const auto* v : { &b.excluded, &b.included, &b.ungated })
            {
                g.setColour (series==&gatedShortTerm?(category==0?juce::Colour(0xff75899e):juce::Colour(0xff62c9ff)):
                    preview&&series==&plot.momentary ? juce::Colour(0xffffce80) : category == 0 ? juce::Colour (0xff89909b) : category == 1 ? juce::Colour (0xffffb65c) : juce::Colour (0xff62c9ff));
                auto& last=previous[size_t(category++)];if(!v->count)continue;
                const auto point=v->count==1?x(double(v->firstEnd)):px;
                g.drawLine(point,y(v->maximum),point,y(v->minimum)+1.f,1.3f);
                // Only join adjacent actual 100 ms windows in the same category.
                // A zero-energy window, missing time or gate change breaks the line.
                if(last&&v->follows(*last,programmeRate()/10))g.drawLine(x(double(last->lastEnd)),y(last->lastValue),x(double(v->firstEnd)),y(v->firstValue),1.2f);
                last=v;
            }
        }
    }
    g.endTransparencyLayer();g.beginTransparencyLayer(blend.waveformAlpha());
    if(analysis&&!preview&&!previewBase) for (const auto& peak : analysis->peaks) for (bool truePeak : { false, true }) if (truePeak ? peak.trueDbtp.has_value() : peak.sampleDbfs.has_value())
    {
        const auto frame = truePeak ? peak.truePeakEighthFrame / 8. : static_cast<double> (peak.sampleFrame);
        if (frame >= viewBegin && frame <= viewEnd)
        {
            const bool exceeds=truePeak&&peakOverlayActive()&&*peak.trueDbtp>std::min(peakTarget.value_or(0.),currentInspection?peakTarget.value_or(0.):0.);
            const auto top = exceeds&&waveformData?graph.getY()+peak.channel*graph.getHeight()/waveformData->spec.channels:graph.getY()+(truePeak?0:16);
            g.setColour (exceeds?juce::Colour(0xffff9d42):truePeak ? juce::Colours::pink : juce::Colours::lightgrey);
            g.drawVerticalLine (static_cast<int> (x(frame)), static_cast<float> (top), static_cast<float> (top + 13));
            g.drawText (juce::String (exceeds?"TP!":truePeak ? "TP" : "SP") + juce::String (peak.channel + 1), static_cast<int> (x(frame)) + 2, top, 48, 14, juce::Justification::left);
        }
    }
    g.endTransparencyLayer();g.beginTransparencyLayer(blend.loudnessAlpha());
    g.setColour (juce::Colours::lightgreen.withAlpha (preview?0.f:.75f)); g.drawHorizontalLine (static_cast<int> (y(-70)), static_cast<float> (graph.getX()), static_cast<float> (graph.getRight()));
    if (!preview && analysis && analysis->integratedGate)
    {
        const auto gateY=y(*analysis->integratedGate);
        g.drawHorizontalLine(int(gateY),float(graph.getX()),float(graph.getRight()));
        g.setColour(juce::Colour(0xffff3232));
        for(size_t n=0;n<changedMomentary.size();)
        {
            if(!changedMomentary[n]){++n;continue;}const auto first=n;
            while(n<changedMomentary.size()&&changedMomentary[n])++n;
            g.drawLine(graph.getX()+float(first)*graph.getWidth()/changedMomentary.size(),gateY,
                       graph.getX()+float(n)*graph.getWidth()/changedMomentary.size(),gateY,2.f);
        }
    }
    g.endTransparencyLayer();
    // Source-frame annotations are independent of video offset and plot reduction.
    // Do not clamp offscreen boundaries into fabricated edge markers.
    g.setColour (juce::Colours::violet.withAlpha (.7f));
    int previousBoundaryPixel=graph.getX()-1;
    for(auto it=std::lower_bound(sectionBoundaries.begin(),sectionBoundaries.end(),viewBegin);
        it!=sectionBoundaries.end() && *it<=viewEnd;++it)
    {
        if(*it<=0 || *it>=programmeFrames())continue;
        const int pixel=int(x(double(*it)));
        if(pixel==previousBoundaryPixel)continue;
        g.drawVerticalLine(pixel,float(graph.getY()),float(graph.getBottom()));
        previousBoundaryPixel=pixel;
    }
    if(gesture && sectionGesture && dragging && gestureEpoch==navigationEpoch)
    {
        const auto left=x(double(movingBoundary?dragFrame:frameAt(downPoint.x))),right=x(double(dragFrame));
        g.setColour(juce::Colours::aqua.withAlpha(.20f));
        g.fillRect(std::min(left,right),float(graph.getY()),std::max(1.f,std::abs(right-left)),float(graph.getHeight()));
        g.setColour(juce::Colours::aqua);
        g.drawVerticalLine(int(right),float(graph.getY()),float(graph.getBottom()));
        g.drawText("Frame "+juce::String(dragFrame),graph.reduced(5),juce::Justification::topLeft);
    }
    if(preview&&previewCoverage<viewEnd)
    {
        const auto left=std::clamp(x(double(previewCoverage)),float(graph.getX()),float(graph.getRight()));
        g.setColour(juce::Colours::lightgrey);g.drawText("Loudness pending",int(left),graph.getY()+5,std::max(0,graph.getRight()-int(left)),18,juce::Justification::centred);
    }
    if(onRequestSeek&&playhead>=viewBegin&&playhead<=viewEnd)
    {g.setColour(juce::Colours::white.withAlpha(.8f));g.drawVerticalLine(int(x(double(playhead))),float(graph.getY()),float(graph.getBottom()));}
    g.restoreState(); g.setColour (juce::Colours::lightgrey);
    for (int n = 0; n <= 4; ++n)
    {
        const auto f = viewBegin + (viewEnd - viewBegin) * n / 4;
        const auto px = juce::jlimit (graph.getX(), graph.getRight() - 85, static_cast<int> (x(f)) - (n == 4 ? 85 : 0));
        g.drawText (juce::String (f / programmeRate(), (viewEnd-viewBegin)/programmeRate()<.01?6:3) + " s", px, graph.getBottom() + 2, 90, 20, juce::Justification::left);
    }
}
void MeasurementTimeline::paintWaveform(juce::Graphics& g)
{
    const int channels=waveformData->spec.channels;
    double scale=1;for(const auto& col:waveformColumns)if(col.available())for(int c=0;c<channels;++c)
        scale=std::max(scale,double(std::max(std::abs(col.channels[size_t(c)].minimum),std::abs(col.channels[size_t(c)].maximum))));
    const auto x=[&](double f){return float(graph.getX()+(f-viewBegin)*(graph.getWidth()-1)/(viewEnd-viewBegin));};
    const float height=float(graph.getHeight())/channels;
    for(int c=0;c<channels;++c)
    {
        const float top=graph.getY()+c*height,centre=top+height/2;
        const auto name=channelNames[c].isEmpty()?juce::String(c+1):channelNames[c];
        g.setColour(juce::Colours::lightgrey);
        if(height>=30){g.drawText(name,2,int(centre)-14,78,16,juce::Justification::right);g.drawText("+/-"+juce::String(scale,2),2,int(centre)+1,78,15,juce::Justification::right);}
        else g.drawText(name+" +/-"+juce::String(scale,2),2,int(centre)-7,78,15,juce::Justification::right);
        g.setColour(juce::Colour(0xff394351));g.drawHorizontalLine(int(centre),float(graph.getX()),float(graph.getRight()));
        g.drawHorizontalLine(int(top+height),float(graph.getX()),float(graph.getRight()));
        const auto y=[&](float v){return centre-float(v/scale)*(height*.44f);};
        g.saveState();g.reduceClipRegion(graph);
        if(peakOverlayActive()&&(waveformAccepted||previewBase))
        {
            const auto level=peakCeiling;g.setColour(juce::Colour(0xffff9d42).withAlpha(.45f));
            if(level<=scale)for(auto sign:{-1.,1.})for(int px=graph.getX();px<graph.getRight();px+=9)
                g.drawLine(float(px),y(float(sign*level)),float(std::min(px+4,graph.getRight())),y(float(sign*level)),1.f);
        }
        for(const auto& col:waveformColumns)if(col.available())
        {
            const auto& e=col.channels[size_t(c)];const float midpoint=x((double(col.begin)+col.coveredEnd)/2);
            g.setColour(juce::Colour(0xff697b8c));g.drawLine(midpoint,y(e.maximum),midpoint,y(e.minimum)+.8f,1.3f);
            // Endpoints retain their actual sample coordinates even after reduction.
            if(x(double(col.coveredEnd))-x(double(col.begin))>3)
            {
                g.setColour(juce::Colour(0xffa7b7c5));
                g.fillEllipse(x(double(e.minimumFrame))-1,y(e.minimum)-1,2,2);
                g.fillEllipse(x(double(e.maximumFrame))-1,y(e.maximum)-1,2,2);
            }
            g.setColour(juce::Colour(0xffff9d42));
            if(sampleExceedsTarget(e.minimum,e.minimumFrame))g.drawEllipse(x(double(e.minimumFrame))-3,y(e.minimum)-3,6,6,1.8f);
            if(e.maximumFrame!=e.minimumFrame&&sampleExceedsTarget(e.maximum,e.maximumFrame))g.drawEllipse(x(double(e.maximumFrame))-3,y(e.maximum)-3,6,6,1.8f);
        }
        g.restoreState();
    }
    const auto coverage=waveformData->coveredFrames(waveformPublication);
    if(coverage<programmeFrames() && coverage<viewEnd)
    {
        const int begin=int(std::clamp(x(double(coverage)),float(graph.getX()),float(graph.getRight())));
        auto pending=graph.withLeft(begin);g.setColour(juce::Colour(0xff20252c).withAlpha(.85f));g.fillRect(pending);
        g.setColour(juce::Colours::lightgrey);g.drawText("Pending",pending.reduced(3),juce::Justification::centred);
    }
}

}
