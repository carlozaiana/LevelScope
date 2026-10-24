#pragma once
#include "TraceBuffer.h"
#include "Standalone/StandaloneVideo.h"
#include <juce_gui_basics/juce_gui_basics.h>

namespace levelscope::paintdiagnostic
{
class Trace
{
public:
    Trace();
    ~Trace();
    bool start(juce::Component&);
    void mark(Kind) noexcept;
    void selection(standalone::video::Player&, const standalone::video::Player::View&,
                   standalone::video::Time);
    void stop();
    juce::String save(const juce::File&, double sampleRate);
    bool valid() const noexcept;
    juce::String error() const;
    struct Span
    {
        Trace& trace; Kind end;
        Span(Trace& t, Kind begin, Kind finish) : trace(t),end(finish) { trace.mark(begin); }
        ~Span() { trace.mark(end); }
    };
private:
    struct Native;
    friend struct Native;
    Buffer<131072> ui;
    Buffer<32768> display;
    Event last;
    bool stopped=false, started=false;
    uint64_t mainVblanks=0;
    juce::String failure;
    std::unique_ptr<Native> native;
    juce::VBlankAttachment vblank;
};
}
