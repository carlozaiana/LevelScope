#include "MeasurementTimeline.h"
#include <cmath>
#include <cstdlib>

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
    for (auto* b : { &zoomIn, &zoomOut, &fit, &inspect }) addAndMakeVisible (*b);
    addAndMakeVisible (seconds); seconds.setName ("Audio-source inspection seconds"); seconds.setText ("0");
    addAndMakeVisible (scroll); scroll.addListener (this); scroll.setName ("Audio-source horizontal scroll");
    zoomIn.onClick = [this] { zoom (.5); }; zoomOut.onClick = [this] { zoom (2); };
    fit.onClick = [this] { if (analysis) setViewRange (0, static_cast<double> (analysis->frames)); };
    inspect.onClick = [this]
    {
        const auto value = seconds.getText().trim().toStdString(); char* end = nullptr;
        const auto n = std::strtod (value.c_str(), &end);
        if (value.empty() || end != value.c_str() + value.size() || ! inspectSeconds (n))
        { inspection = "Enter finite audio-source seconds within the programme."; if (onInspect) onInspect (inspection); }
    };
    seconds.onReturnKey = [this] { inspect.onClick(); };
}
void MeasurementTimeline::setAnalysis (std::shared_ptr<const proposal::ProgramAnalysis> value)
{
    if (value == analysis) return;
    analysis = std::move (value); selection.reset(); inspection.clear();
    viewBegin = 0; viewEnd = analysis ? static_cast<double> (analysis->frames) : 1;
    rebuild();
}
void MeasurementTimeline::setViewRange (double begin, double end)
{
    if (! analysis || ! std::isfinite (begin) || ! std::isfinite (end) || end <= begin) return;
    const auto width = std::clamp (end - begin, std::min (static_cast<double> (analysis->frames), analysis->sampleRate * .1), static_cast<double> (analysis->frames));
    viewBegin = std::clamp (begin, 0., static_cast<double> (analysis->frames) - width); viewEnd = viewBegin + width;
    rebuild();
}
void MeasurementTimeline::zoom (double factor)
{
    const auto centre = (viewBegin + viewEnd) / 2, half = (viewEnd - viewBegin) * factor / 2;
    setViewRange (centre - half, centre + half);
}
void MeasurementTimeline::highlight (double begin, double end, bool reveal)
{
    if (! analysis) return;
    selection = std::make_pair (begin, end);
    const auto duration = static_cast<double> (analysis->frames);
    if (reveal && (begin < viewBegin || end > viewEnd))
    {
        const auto width = std::max (viewEnd - viewBegin, std::min (duration, analysis->sampleRate * 5.));
        setViewRange (std::clamp ((begin + end) / 2, 0., duration) - width / 2,
                      std::clamp ((begin + end) / 2, 0., duration) + width / 2);
    }
    repaint();
}
bool MeasurementTimeline::inspectSeconds (double value)
{
    if (! analysis || ! std::isfinite (value) || value < 0 || value > static_cast<double> (analysis->frames) / analysis->sampleRate) return false;
    const auto frame = value * analysis->sampleRate;
    inspection = "Audio-source " + juce::String (value, 6) + " s; nearest sampled window ends (earlier wins ties).\n";
    for (bool shortTerm : { false, true })
    {
        const auto& windows = shortTerm ? analysis->shortTerm : analysis->momentary;
        const auto n = proposal::inspectWindow (windows, frame);
        inspection << (shortTerm ? "S: " : "M: ");
        if (n)
        {
            inspection << diagnosticWindowText (windows[*n], analysis->sampleRate, analysis->frames, shortTerm) << "\n";
            if (! shortTerm) highlight (static_cast<double> (windows[*n].begin), static_cast<double> (windows[*n].end));
        }
        else inspection << "no complete windows\n";
    }
    if (onInspect) onInspect (inspection);
    return true;
}
void MeasurementTimeline::scrollBarMoved (juce::ScrollBar*, double start)
{ if (! updatingScroll) setViewRange (start, start + viewEnd - viewBegin); }
void MeasurementTimeline::rebuild()
{
    plot = {};
    if (analysis) plot = proposal::prepareTimeline (*analysis, viewBegin, viewEnd, std::max (1, graph.getWidth()));
    updatingScroll = true;
    scroll.setRangeLimits (0, analysis ? static_cast<double> (analysis->frames) : 1);
    scroll.setCurrentRange (viewBegin, viewEnd - viewBegin, juce::dontSendNotification);
    updatingScroll = false;
    for (auto* b : { &zoomIn, &zoomOut, &fit, &inspect }) b->setEnabled (analysis != nullptr);
    seconds.setEnabled (analysis != nullptr); scroll.setEnabled (analysis != nullptr);
    repaint();
}
void MeasurementTimeline::resized()
{
    auto b = getLocalBounds(); auto row = b.removeFromTop (28);
    zoomIn.setBounds (row.removeFromLeft (80)); zoomOut.setBounds (row.removeFromLeft (80));
    fit.setBounds (row.removeFromLeft (155)); row.removeFromLeft (12); seconds.setBounds (row.removeFromLeft (110)); inspect.setBounds (row.removeFromLeft (145));
    scroll.setBounds (b.removeFromBottom (16).withTrimmedLeft (55));
    b.removeFromTop (44); b.removeFromBottom (25); graph = b.withTrimmedLeft (55).withTrimmedRight (12);
    rebuild();
}
void MeasurementTimeline::mouseDown (const juce::MouseEvent& e)
{
    if (analysis && graph.contains (e.getPosition())) inspectSeconds ((viewBegin + (e.position.x - graph.getX()) / graph.getWidth() * (viewEnd - viewBegin)) / analysis->sampleRate);
}
void MeasurementTimeline::mouseWheelMove (const juce::MouseEvent&, const juce::MouseWheelDetails& w)
{ const auto shift = -(w.deltaX != 0 ? w.deltaX : w.deltaY) * (viewEnd - viewBegin); setViewRange (viewBegin + shift, viewEnd + shift); }
void MeasurementTimeline::paint (juce::Graphics& g)
{
    g.fillAll (juce::Colour (0xff171c24)); g.setFont (12.0f); g.setColour (juce::Colours::white);
    g.drawText ("Audio-source time | M included: orange; M excluded: grey; S: blue (ungated)", 55, 31, getWidth() - 60, 18, juce::Justification::left);
    g.drawText ("Points at window ends, every 100 ms. Per-column min/max retains extrema and mixed membership; no interpolation across gaps.", 55, 49, getWidth() - 60, 18, juce::Justification::left);
    if (! analysis || graph.isEmpty()) { g.drawText ("No measured timeline", getLocalBounds(), juce::Justification::centred); return; }
    double low = -75, high = -10;
    for (const auto* series : { &plot.momentary, &plot.shortTerm }) for (const auto& b : *series)
        for (const auto* v : { &b.included, &b.excluded, &b.ungated }) if (v->count)
        { low = std::min (low, v->minimum - 2); high = std::max (high, v->maximum + 2); }
    const auto y = [&] (double value) { return static_cast<float> (graph.getBottom() - (value - low) * graph.getHeight() / (high - low)); };
    const auto x = [&] (double frame) { return static_cast<float> (graph.getX() + (frame - viewBegin) * graph.getWidth() / (viewEnd - viewBegin)); };
    for (double v = std::ceil (low / 10) * 10; v <= high; v += 10)
    {
        g.setColour (juce::Colour (0xff303744)); g.drawHorizontalLine (static_cast<int> (y(v)), static_cast<float> (graph.getX()), static_cast<float> (graph.getRight()));
        g.setColour (juce::Colours::lightgrey); g.drawText (juce::String (v, 0), 0, static_cast<int> (y(v)) - 7, 48, 15, juce::Justification::right);
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
        for (size_t n = 0; n < series->size(); ++n)
        {
            const auto px = static_cast<float> (graph.getX()) + (static_cast<float> (n) + .5f) * graph.getWidth() / static_cast<float> (series->size());
            const auto& b = (*series)[n]; int category = 0;
            for (const auto* v : { &b.excluded, &b.included, &b.ungated })
            {
                g.setColour (category == 0 ? juce::Colour (0xff89909b) : category == 1 ? juce::Colour (0xffffb65c) : juce::Colour (0xff62c9ff)); ++category;
                if (v->count) g.drawLine (px, y(v->maximum), px, y(v->minimum) + 1.f, 1.3f);
            }
        }
    for (const auto& peak : analysis->peaks) for (bool truePeak : { false, true }) if (truePeak ? peak.trueDbtp.has_value() : peak.sampleDbfs.has_value())
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
    g.setColour (juce::Colours::lightgreen.withAlpha (.75f)); g.drawHorizontalLine (static_cast<int> (y(-70)), static_cast<float> (graph.getX()), static_cast<float> (graph.getRight()));
    if (analysis->integratedGate) g.drawHorizontalLine (static_cast<int> (y(*analysis->integratedGate)), static_cast<float> (graph.getX()), static_cast<float> (graph.getRight()));
    g.restoreState(); g.setColour (juce::Colours::lightgrey);
    for (int n = 0; n <= 4; ++n)
    {
        const auto f = viewBegin + (viewEnd - viewBegin) * n / 4;
        const auto px = juce::jlimit (graph.getX(), graph.getRight() - 85, static_cast<int> (x(f)) - (n == 4 ? 85 : 0));
        g.drawText (juce::String (f / analysis->sampleRate, 3) + " s", px, graph.getBottom() + 2, 90, 20, juce::Justification::left);
    }
}
}
