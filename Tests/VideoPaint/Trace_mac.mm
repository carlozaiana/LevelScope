#include "Trace.h"
#include "Standalone/StandaloneAudition.h"
#import <Cocoa/Cocoa.h>
#import <CoreVideo/CoreVideo.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstdio>

namespace levelscope::standalone::video
{
// The sole Player change is this test friend's declaration. It neither changes
// object layout nor production selection/decoding. Extra locking is diagnostic.
struct DiagnosticPlayerProbe
{
    static void sample(Player& p, paintdiagnostic::Event& e)
    {
        std::lock_guard<std::mutex> lock(p.mutex);
        e.generation=p.generation; e.epoch=p.viewEpoch.load();
        e.desired=p.desired ? static_cast<int64_t>(*p.desired) : -1;
        e.queueCount=p.queue.size();
        e.queueFirst=p.queue.empty() ? -1 : static_cast<int64_t>(p.queue.front().index);
        e.queueLast=p.queue.empty() ? -1 : static_cast<int64_t>(p.queue.back().index);
    }
};
}
namespace levelscope::paintdiagnostic
{
using standalone::StandaloneAudition;
struct Trace::Native
{
    Trace& owner;
    NSWindow* window=nil; // not owned; Trace stops before its Component is destroyed
    NSMutableArray* observers=nil;
    CVDisplayLinkRef link=nullptr;
    CGDirectDisplayID displayId=0;
    uint32_t previousFlags=UINT32_MAX;
    uint64_t tickCount=0;
    double nominalHz=0;
    int pixelWidth=0,pixelHeight=0;
    Native(Trace& t, juce::Component& c) : owner(t)
    {
        if (auto* peer=c.getPeer()) window=[(NSView*)peer->getNativeHandle() window];
        if (window==nil) { owner.failure="Diagnostic window has no native peer"; return; }
        observers=[[NSMutableArray alloc] init];
        NSArray* names=@[NSWindowDidChangeOcclusionStateNotification,NSWindowDidBecomeKeyNotification,
            NSWindowDidResignKeyNotification,NSWindowDidChangeScreenNotification,NSWindowDidMiniaturizeNotification,
            NSWindowDidDeminiaturizeNotification,NSWindowWillStartLiveResizeNotification,NSWindowDidEndLiveResizeNotification,
            NSApplicationDidBecomeActiveNotification,NSApplicationDidResignActiveNotification];
        for (NSString* name in names)
        {
            id token=[[NSNotificationCenter defaultCenter] addObserverForName:name object:nil queue:nil
                usingBlock:^(NSNotification* note) {
                    if ([note object]==window || [note object]==NSApp) capture();
                }];
            [observers addObject:token];
        }
        capture();
    }
    ~Native()
    {
        for (id token in observers) [[NSNotificationCenter defaultCenter] removeObserver:token];
        [observers release];
        if (link) { CVDisplayLinkStop(link); CVDisplayLinkRelease(link); }
    }
    void capture()
    {
        if (owner.stopped) return;
        auto screen=window.screen;
        auto id=static_cast<CGDirectDisplayID>([[screen.deviceDescription objectForKey:@"NSScreenNumber"] unsignedIntValue]);
        uint32_t flags=([window isVisible]?1u:0u) | (([window occlusionState]&NSWindowOcclusionStateVisible)?2u:0u)
            | ([NSApp isActive]?4u:0u) | ([window isKeyWindow]?8u:0u)
            | ([window isMiniaturized]?16u:0u) | ([window inLiveResize]?32u:0u);
        if (id!=displayId)
        {
            if (link) { CVDisplayLinkStop(link); CVDisplayLinkRelease(link); link=nullptr; }
            displayId=id;
            if (auto mode=CGDisplayCopyDisplayMode(id))
            {
                nominalHz=CGDisplayModeGetRefreshRate(mode);
                pixelWidth=static_cast<int>(CGDisplayModeGetPixelWidth(mode));
                pixelHeight=static_cast<int>(CGDisplayModeGetPixelHeight(mode));
                CGDisplayModeRelease(mode);
            }
            auto result=CVDisplayLinkCreateWithCGDisplay(id,&link);
            if (result==kCVReturnSuccess) result=CVDisplayLinkSetOutputCallback(link,tick,this);
            if (result==kCVReturnSuccess) result=CVDisplayLinkStart(link);
            if (result!=kCVReturnSuccess) owner.failure="Independent display link failed: "+juce::String(result);
        }
        if (id!=owner.last.display || flags!=previousFlags)
        {
            previousFlags=flags; owner.last.display=id; owner.last.flags=flags;
            owner.mark(Kind::windowState);
        }
    }
    static CVReturn tick(CVDisplayLinkRef,const CVTimeStamp*,const CVTimeStamp* output,
                         CVOptionFlags,CVOptionFlags*,void* context)
    {
        auto& self=*static_cast<Native*>(context);
        Event event; event.hostNs=StandaloneAudition::hostTimeNs(); event.kind=Kind::displayTick;
        event.display=self.displayId; event.generation=++self.tickCount;
        // Raw predicted display host ticks are supplemental, not actual scanout.
        event.auxiliary=output ? output->hostTime : 0;
        self.owner.display.push(event);
        return kCVReturnSuccess;
    }
};
Trace::Trace()=default;
Trace::~Trace() { stop(); }
bool Trace::start(juce::Component& c)
{
    if (started) return failure.isEmpty();
    started=true;
    native=std::make_unique<Native>(*this,c);
    vblank=juce::VBlankAttachment(&c,[this] { ++mainVblanks; mark(Kind::mainVblank); });
    return failure.isEmpty();
}
void Trace::mark(Kind kind) noexcept
{
    if (stopped) return;
    auto event=last; event.kind=kind; event.hostNs=StandaloneAudition::hostTimeNs();
    if (kind==Kind::mainVblank) event.auxiliary=mainVblanks;
    ui.push(event);
}
void Trace::selection(standalone::video::Player& player,const standalone::video::Player::View& view,
                      standalone::video::Time time)
{
    standalone::video::DiagnosticPlayerProbe::sample(player,last);
    last.timeNumerator=time.numerator; last.timeDenominator=time.denominator;
    last.picture=view.picture ? static_cast<int64_t>(view.picture->index) : -1;
    mark(Kind::selection);
}
void Trace::stop()
{
    if (stopped) return;
    mark(Kind::stopping); stopped=true;
    vblank=juce::VBlankAttachment();
    // Stop the producer before reading its buffer, but retain metadata for save.
    if (native && native->link) { CVDisplayLinkStop(native->link); }
}
bool Trace::valid() const noexcept
{
    return stopped && failure.isEmpty() && ui.lost()==0 && display.lost()==0
        && ui.size()>0 && display.size()>0 && mainVblanks>0;
}
juce::String Trace::error() const { return failure; }
juce::String Trace::save(const juce::File& path,double rate)
{
    stop();
    int fd=::open(path.getFullPathName().toRawUTF8(),O_WRONLY|O_CREAT|O_EXCL,0600);
    if (fd<0) return "Could not exclusively create diagnostic trace";
    FILE* stream=fdopen(fd,"w");
    if (!stream) { ::close(fd); return "Could not open diagnostic trace stream"; }
    fprintf(stream,"{\"schema\":\"levelscope-paint-trace-v1\",\"scope\":\"Diagnostic callbacks; not physical presentation\",\"valid\":%s,\"sample_rate\":%.0f,\"ui_dropped\":%llu,\"display_dropped\":%llu,\"allocated_buffer_bytes\":%zu,\"main_vblanks\":%llu,\"display_id\":%u,\"nominal_hz\":%.6f,\"pixel_width\":%d,\"pixel_height\":%d,",
        valid()?"true":"false",rate,(unsigned long long)ui.lost(),(unsigned long long)display.lost(),
        Buffer<131072>::bytes+Buffer<32768>::bytes,(unsigned long long)mainVblanks,
        native?native->displayId:0,native?native->nominalHz:0,native?native->pixelWidth:0,native?native->pixelHeight:0);
    fputs("\"kinds\":{\"1\":\"timer_enter\",\"2\":\"timer_exit\",\"3\":\"selection\",\"4\":\"repaint_request\",\"5\":\"paint_enter\",\"6\":\"paint_exit\",\"7\":\"main_vblank\",\"8\":\"display_tick\",\"9\":\"window_state\",\"10\":\"stopping\"},\"window_flags\":{\"1\":\"visible\",\"2\":\"unoccluded\",\"4\":\"app_active\",\"8\":\"key_window\",\"16\":\"minimized\",\"32\":\"live_resize\"},\"columns\":[\"host_ns\",\"kind\",\"window_flags\",\"time_numerator\",\"time_denominator\",\"generation_or_tick_sequence\",\"view_epoch\",\"desired\",\"picture\",\"queue_first\",\"queue_last\",\"queue_count\",\"display_id\",\"auxiliary\"],",stream);
    auto write=[&](const char* name,const auto& buffer)
    {
        fprintf(stream,"\"%s\":[",name);
        for(size_t i=0;i<buffer.size();++i)
        {
            const auto& e=buffer.at(i);
            fprintf(stream,"%s[%llu,%u,%u,%lld,%lld,%llu,%llu,%lld,%lld,%lld,%lld,%llu,%llu,%llu]",i?",":"",
                (unsigned long long)e.hostNs,(unsigned)e.kind,e.flags,(long long)e.timeNumerator,(long long)e.timeDenominator,
                (unsigned long long)e.generation,(unsigned long long)e.epoch,(long long)e.desired,(long long)e.picture,
                (long long)e.queueFirst,(long long)e.queueLast,(unsigned long long)e.queueCount,
                (unsigned long long)e.display,(unsigned long long)e.auxiliary);
        }
        fputs("]",stream);
    };
    write("ui_events",ui); fputs(",",stream); write("display_events",display); fputs("}\n",stream);
    const bool failed=ferror(stream)!=0;
    const bool closed=fclose(stream)==0;
    return failed || !closed ? "Failed to write complete diagnostic trace" : juce::String();
}
}
