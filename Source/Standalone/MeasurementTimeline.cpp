#include "MeasurementTimeline.h"
#include "TimelineGainLane.h"
#include <cmath>
#include <cstdlib>
#include "VideoTiming.h"

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
void MeasurementTimeline::setProgramme(int64_t frames,int rate,const juce::String& identity)
{
    if(identity==programmeIdentity&&frames==knownFrames&&rate==knownRate)return;
    programmeIdentity=identity;knownFrames=std::max<int64_t>(0,frames);knownRate=std::max(1,rate);
    analysis.reset();preview.reset();waveformData.reset();selection.reset();inspectedRange.reset();inspectedFrame.reset();inspection.clear();seconds.setText("0",false);
    viewBegin=0;viewEnd=std::max<double>(1,knownFrames);++navigationEpoch;rebuild();
}
void MeasurementTimeline::setWaveform(std::shared_ptr<const waveform::Overview> data,bool accepted,const juce::StringArray& channels)
{
    const auto publication=data?data->publication():0;
    if(waveformEnabled && data==waveformData && publication==waveformPublication && accepted==waveformAccepted && channels==channelNames)return;
    const auto previousFrames=programmeFrames();const auto previousRate=programmeRate();
    waveformEnabled=true;waveformData=std::move(data);waveformAccepted=accepted;waveformPublication=publication;channelNames=channels;
    if(programmeFrames()!=previousFrames || programmeRate()!=previousRate || viewEnd<=1)
    {viewBegin=0;viewEnd=std::max<double>(1,programmeFrames());selection.reset();}
    inspectedChannel=waveformData?juce::jlimit(0,waveformData->spec.channels-1,inspectedChannel):0;
    resized();
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
    if(key!=juce::KeyPress::spaceKey||!onTogglePlayback||juce::Component::getCurrentlyModalComponent())return false;
    for(auto* c=juce::Component::getCurrentlyFocusedComponent();c;c=c->getParentComponent())
        if(dynamic_cast<juce::TextEditor*>(c))return false;
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
        if(onInspect)onInspect("Current section gain settings  |  "+gainLane->inspectionText());
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
    if(waveformData)
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
    if(waveformData)waveformPublication=waveformData->publication();
    waveformColumns=waveformData?waveformData->reduce(int64_t(std::floor(viewBegin)),int64_t(std::ceil(viewEnd)),std::max(1,graph.getWidth()),waveformPublication):std::vector<waveform::Column>{};
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
    scroll.setBounds (b.removeFromBottom (16).withTrimmedLeft (55));
    if(gainLane)gainLane->setBounds(b.removeFromBottom(TimelineGainLane::height));
    b.removeFromTop (44); ruler=b.removeFromTop(waveformEnabled?24:0).withTrimmedLeft(waveformEnabled?85:55).withTrimmedRight(waveformEnabled?60:12); b.removeFromBottom (25); graph = b.withTrimmedLeft(waveformEnabled?85:55).withTrimmedRight(waveformEnabled?60:12);
    if(graph!=previousGraph)++navigationEpoch;
    rebuild();
}
void MeasurementTimeline::bindSectionEditing(StandaloneSessionModel& session,std::function<void(uint64_t)> changed)
{
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
        && !session.sourceContentInvalid && !session.audioEditsBlocked() && !session.job.isActive();
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
void MeasurementTimeline::mouseMove(const juce::MouseEvent& e)
{
    setMouseCursor(graph.contains(e.getPosition()) && boundaryHandle(e.position.x)
        ? juce::MouseCursor::LeftRightResizeCursor : juce::MouseCursor::NormalCursor);
}
void MeasurementTimeline::finishSectionGesture(int64_t end)
{
    refreshSectionEditing();
    if(!sectionEditable || gestureEpoch!=navigationEpoch || !sectionSession->sourceDocument.matchesFileStamps())return;
    auto& document=sectionSession->sectionDocument;
    const auto& map=document.map();
    juce::String error;
    uint64_t selected=0;
    const auto before=map.revision;
    bool accepted=false;
    if(movingBoundary)
    {
        const auto n=*movingBoundary;
        selected=map.sections[n].id;
        accepted=document.moveEnd(selected,std::clamp(end,map.sections[n].begin+1,map.sections[n+1].end-1),error);
    }
    else
    {
        const auto start=frameAt(downPoint.x);
        if(start==end)return;
        accepted=document.defineInterval(std::min(start,end),std::max(start,end),selected,error);
    }
    if(!accepted){navigationMessage(error);return;}
    if(map.revision!=before)sectionSession->documentChanged();
    refreshSectionEditing();
    if(sectionEdited)sectionEdited(selected);
    const auto it=std::find_if(map.sections.begin(),map.sections.end(),[selected](const auto& row){return row.id==selected;});
    if(it!=map.sections.end())highlight(double(it->begin),double(it->end),false);
    navigationMessage("Section saved. Drag a boundary to adjust it.");
}
void MeasurementTimeline::mouseDown(const juce::MouseEvent& e)
{
    refreshSectionEditing();
    gestureGain=gainLane&&gainLane->plotBounds().contains(e.getPosition()-gainLane->getPosition());
    gesture=e.mods.isLeftButtonDown()&&programmeFrames()>0&&(graph.contains(e.getPosition())||ruler.contains(e.getPosition())||gestureGain);
    sectionGesture=gesture && sectionEditable && graph.contains(e.getPosition()) && !e.mods.isAltDown();
    movingBoundary=sectionGesture?boundaryHandle(e.position.x):std::nullopt;
    dragFrame=frameAt(e.position.x);
    dragging=false;gestureRuler=ruler.contains(e.getPosition());downPoint=e.position;dragBegin=viewBegin;dragEnd=viewEnd;gestureEpoch=navigationEpoch;gestureIdentity=transportIdentity;
}
void MeasurementTimeline::mouseDrag(const juce::MouseEvent& e)
{
    if(!gesture)return;
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
    if(dragging||e.position.getDistanceFrom(downPoint)>=4||gestureEpoch!=navigationEpoch)return;
    if(gestureRuler&&ruler.contains(e.getPosition()))requestSeek(frameAt(e.position.x),gestureIdentity);
    else if(gestureGain&&gainLane&&gainLane->plotBounds().contains(e.getPosition()-gainLane->getPosition()))
        gainLane->selectFrame(frameAt(e.position.x),true);
    else if(!gestureRuler&&graph.contains(e.getPosition()))
    {
        if(waveformData)inspectedChannel=juce::jlimit(0,waveformData->spec.channels-1,(e.y-graph.getY())*waveformData->spec.channels/std::max(1,graph.getHeight()));
        inspectFrame(double(frameAt(e.position.x)));
    }
}
void MeasurementTimeline::mouseWheelMove(const juce::MouseEvent& e,const juce::MouseWheelDetails& w)
{
    follow.setToggleState(false,juce::dontSendNotification);
    if(e.mods.isCommandDown()||e.mods.isCtrlDown())zoom(w.deltaY>0?.8:1.25);
    else {const auto shift=-(w.deltaX!=0?w.deltaX:w.deltaY)*(viewEnd-viewBegin);setViewRange(viewBegin+shift,viewEnd+shift);}
}
void MeasurementTimeline::paint (juce::Graphics& g)
{
    g.fillAll (juce::Colour (0xff171c24)); g.setFont (12.0f); g.setColour (juce::Colours::white);
    if(waveformEnabled)
    {
        juce::String heading=waveformData?"PCM waveform overview | "+juce::String(waveformData->binFrames)+" frames/bin ("+juce::String(1000.*waveformData->binFrames/programmeRate(),3)+" ms)":
            "PCM waveform unavailable | use Build waveform for an existing measurement";
        if(waveformData && !waveformAccepted)heading+=" | ANALYSING - incomplete: "+juce::String(waveformData->coveredFrames(waveformPublication))+" / "+juce::String(programmeFrames())+" frames";
        g.drawText(heading,10,31,getWidth()-20,18,juce::Justification::left);
        g.drawText(preview?"ANALYSING - incomplete, not a programme result | M amber: gate pending; S blue: ungated | LUFS at right, PCM at left.":"Loudness overlay: measured M orange (excluded grey), S blue; window-end times. LUFS scale at right; PCM amplitude at left.",10,49,getWidth()-20,18,juce::Justification::left);
    }
    else
    {
        g.drawText("Audio-source time | M included: orange; M excluded: grey; S: blue (ungated)",55,31,getWidth()-60,18,juce::Justification::left);
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
    if(waveformData)paintWaveform(g);
    double low = -75, high = -10;
    for (const auto* series : { &plot.momentary, &plot.shortTerm }) for (const auto& b : *series)
        for (const auto* v : { &b.included, &b.excluded, &b.ungated }) if (v->count)
        { low = std::min (low, v->minimum - 2); high = std::max (high, v->maximum + 2); }
    const auto y = [&] (double value) { return static_cast<float> (graph.getBottom() - (value - low) * graph.getHeight() / (high - low)); };
    const auto x = [&] (double frame) { return static_cast<float> (graph.getX() + (frame - viewBegin) * (graph.getWidth()-1) / (viewEnd - viewBegin)); };
    for (double v = std::ceil (low / 10) * 10; v <= high; v += 10)
    {
        g.setColour (juce::Colour (0xff303744).withAlpha(waveformData?.3f:1.f)); g.drawHorizontalLine (static_cast<int> (y(v)), static_cast<float> (graph.getX()), static_cast<float> (graph.getRight()));
        g.setColour (juce::Colours::lightgrey); g.drawText(juce::String(v,0)+(waveformEnabled?" LUFS":""),waveformEnabled?graph.getRight()+3:0,static_cast<int>(y(v))-7,waveformEnabled?57:48,15,juce::Justification::right);
    }
    g.saveState(); g.reduceClipRegion (graph);
    // Intervals outside the viewport must not become fabricated edge markers.
    // Zero-width channel-peak markers are visible only at their actual position.
    if (selection && ((selection->second > viewBegin && selection->first < viewEnd)
        || (selection->first == selection->second && selection->first >= viewBegin && selection->first <= viewEnd)))
    {
        const auto left = x (std::clamp (selection->first, viewBegin, viewEnd));
        const auto right = x (std::clamp (selection->second, viewBegin, viewEnd));
        g.setColour (juce::Colours::yellow.withAlpha (.16f)); g.fillRect (left, static_cast<float> (graph.getY()), std::max (2.f, right - left), static_cast<float> (graph.getHeight()));
    }
    for (const auto* series : { &plot.momentary, &plot.shortTerm })
    {
        std::array<const proposal::ValueExtent*,3> previous{};
        for (size_t n = 0; n < series->size(); ++n)
        {
            const auto px = static_cast<float> (graph.getX()) + (static_cast<float> (n) + .5f) * graph.getWidth() / static_cast<float> (series->size());
            const auto& b = (*series)[n]; int category = 0;
            for (const auto* v : { &b.excluded, &b.included, &b.ungated })
            {
                g.setColour (preview&&series==&plot.momentary ? juce::Colour(0xffffce80) : category == 0 ? juce::Colour (0xff89909b) : category == 1 ? juce::Colour (0xffffb65c) : juce::Colour (0xff62c9ff));
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
    if(analysis&&!preview) for (const auto& peak : analysis->peaks) for (bool truePeak : { false, true }) if (truePeak ? peak.trueDbtp.has_value() : peak.sampleDbfs.has_value())
    {
        const auto frame = truePeak ? peak.truePeakEighthFrame / 8. : static_cast<double> (peak.sampleFrame);
        if (frame >= viewBegin && frame <= viewEnd)
        {
            const auto top = graph.getY() + (truePeak ? 0 : 16);
            g.setColour (truePeak ? juce::Colours::pink : juce::Colours::lightgrey);
            g.drawVerticalLine (static_cast<int> (x(frame)), static_cast<float> (top), static_cast<float> (top + 13));
            g.drawText (juce::String (truePeak ? "TP" : "SP") + juce::String (peak.channel + 1), static_cast<int> (x(frame)) + 2, top, 40, 14, juce::Justification::left);
        }
    }
    g.setColour (juce::Colours::lightgreen.withAlpha (preview?0.f:.75f)); g.drawHorizontalLine (static_cast<int> (y(-70)), static_cast<float> (graph.getX()), static_cast<float> (graph.getRight()));
    if (!preview && analysis && analysis->integratedGate) g.drawHorizontalLine (static_cast<int> (y(*analysis->integratedGate)), static_cast<float> (graph.getX()), static_cast<float> (graph.getRight()));
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
