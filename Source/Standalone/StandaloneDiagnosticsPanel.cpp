#include "StandaloneDiagnosticsPanel.h"

namespace levelscope::standalone
{
StandaloneDiagnosticsPanel::StandaloneDiagnosticsPanel (const StandaloneSessionModel& s) : session (s)
{
    domainBox.addItem ("Measured Source", 1); domainBox.addItem ("Measured Current State", 2);
    domainBox.setSelectedId (1, juce::dontSendNotification); domainBox.setName ("Diagnostic measured state");
    domainBox.onChange = [this] { timerCallback(); };
    rankingBox.addItem ("Integrated contributors (energy share)", 1);
    rankingBox.addItem ("LRA lower tail (at/below P10)", 2);
    rankingBox.addItem ("LRA upper tail (at/above P95)", 3);
    rankingBox.addItem ("Excluded integrated windows", 4);
    rankingBox.addItem ("Per-channel maximum peak locations", 5);
    rankingBox.setSelectedId (1, juce::dontSendNotification); rankingBox.setName ("Diagnostic ranking");
    rankingBox.onChange = [this] { ranking.deselectAllRows(); timeline.clearHighlight(); details.clear(); ranking.updateContent(); ranking.repaint(); };
    rebuildButton.onClick = [this] { controller.retry(); timerCallback(); };
    cancelButton.onClick = [this] { controller.cancel(); timerCallback(); };
    timeline.onInspect = [this] (const juce::String& text) { ranking.deselectAllRows(); details.setText (text, false); };
    for (auto* box : { &summary, &details })
    {
        box->setMultiLine (true); box->setReadOnly (true); box->setFont (13.0f);
        box->setColour (juce::TextEditor::backgroundColourId, juce::Colour (0xff202631));
        box->setColour (juce::TextEditor::textColourId, juce::Colours::white); addAndMakeVisible (*box);
    }
    summary.setName ("Measured diagnostics availability and interpretation"); details.setName ("Full-resolution diagnostic inspection");
    for (auto* c : std::array<juce::Component*, 7> { &domainBox, &rankingBox, &rebuildButton, &cancelButton, &status, &timeline, &ranking }) addAndMakeVisible (*c);
    ranking.setRowHeight (25); ranking.setColour (juce::ListBox::backgroundColourId, juce::Colour (0xff202631));
    setSize (1040, 850); startTimerHz (10);
}
StandaloneDiagnosticsPanel::~StandaloneDiagnosticsPanel() { stopTimer(); controller.cancel(); }
void StandaloneDiagnosticsPanel::suspend() { suspended = true; controller.cancel(); }
void StandaloneDiagnosticsPanel::resume() { suspended = false; controller.retry(); timerCallback(); }
void StandaloneDiagnosticsPanel::timerCallback()
{
    const auto which = domainBox.getSelectedId() == 2 ? MeasurementDomain::currentState : MeasurementDomain::source;
    const auto& measured = which == MeasurementDomain::source ? session.source : session.currentState;
    auto snapshot = ! session.sourceContentInvalid && measured.hasMeasurement
        && (which != MeasurementDomain::currentState || ! session.currentStateNeedsRemeasurement()) ? measured.analysis : nullptr;
    if (! suspended) controller.select (snapshot, session.measurementContext (which), which);
    controller.poll();
    const auto next = controller.result();
    if (next != displayed)
    {
        displayed = next; timeline.setAnalysis (displayed ? displayed->analysis : nullptr);
        ranking.deselectAllRows(); ranking.updateContent(); details.clear(); annotationKey.clear(); refreshSummary(); repaint();
    }
    if (displayed)
    {
        const auto profile = session.getSelectedTargetProfile();
        const auto key = profile.digest() + ";" + juce::String (static_cast<int> (session.programmeApplicability));
        if (key != annotationKey)
        {
            annotationKey = key; annotations = annotateDiagnostics (*displayed, profile, session.programmeApplicability);
            refreshSummary(); ranking.repaint();
        }
    }
    status.setText (controller.status, juce::dontSendNotification);
    cancelButton.setEnabled (controller.busy()); rebuildButton.setEnabled (snapshot != nullptr && ! controller.busy());
}
void StandaloneDiagnosticsPanel::refreshSummary()
{
    if (! displayed) { summary.setText ("No current measured diagnostic result. Select Source or Current State and measure it in the main window. Predictions cannot populate this view.", false); return; }
    const auto& a = *displayed->analysis;
    juce::String s;
    s << (displayed->snapshot->domain == MeasurementDomain::source ? "SOURCE" : "CURRENT STATE") << " | "
      << juce::String (static_cast<double> (a.frames) / a.sampleRate, 6) << " s | I " << diagnosticMetricText (a.integrated, "LUFS")
      << " | LRA " << diagnosticMetricText (a.lra, "LU") << "\n";
    s << diagnosticLraInterpretation (a) << "\n";
    s << "Measured maxima: M " << diagnosticMetricText (a.momentaryMax, "LUFS") << "; S " << diagnosticMetricText (a.shortTermMax, "LUFS")
      << ". Maxima can occur between 100 ms plot points.\n";
    s << "I gates: absolute -70 LUFS; relative " << (a.integratedGate ? juce::String (*a.integratedGate, 6) + " LUFS" : "unavailable")
      << ". Green lines apply to M inclusion only. S timeline is ungated.\n" << annotations.text;
    summary.setText (s, false);
}
int StandaloneDiagnosticsPanel::getNumRows()
{
    if (! displayed) return 0;
    const auto& a = *displayed->analysis;
    switch (rankingBox.getSelectedId())
    {
        case 1: return static_cast<int> (a.integratedRank.size()); case 2: return static_cast<int> (a.lowerTail.size());
        case 3: return static_cast<int> (a.upperTail.size()); case 4: return static_cast<int> (a.excludedIntegrated.size());
        default: return static_cast<int> (a.peaks.size());
    }
}
const proposal::DiagnosticWindow* StandaloneDiagnosticsPanel::windowForRow (int row) const
{
    if (! displayed || row < 0 || rankingBox.getSelectedId() == 5) return nullptr;
    const auto& a = *displayed->analysis;
    const auto kind = rankingBox.getSelectedId();
    const auto& indices = kind == 1 ? a.integratedRank : kind == 2 ? a.lowerTail : kind == 3 ? a.upperTail : a.excludedIntegrated;
    if (static_cast<size_t> (row) >= indices.size()) return nullptr;
    return &(kind == 2 || kind == 3 ? a.lraWindows : a.momentary)[indices[static_cast<size_t> (row)]];
}
juce::String StandaloneDiagnosticsPanel::rowText (int row) const
{
    if (! displayed || row < 0) return {};
    const auto& a = *displayed->analysis; juce::String s;
    if (const auto* w = windowForRow (row))
    {
        s << juce::String (row + 1) << ". " << diagnosticWindowText (*w, a.sampleRate, a.frames);
        if (rankingBox.getSelectedId() == 1) s << " | block-energy share " << juce::String (100 * w->energyShare, 6) << "%";
    }
    else if (rankingBox.getSelectedId() == 5 && static_cast<size_t> (row) < a.peaks.size())
    {
        const auto& p = a.peaks[static_cast<size_t> (row)];
        s << "File channel " << juce::String (row + 1) << " " << juce::AudioChannelSet::getAbbreviatedChannelTypeName (displayed->snapshot->context.fileLayout.roles[row]);
        if (p.sampleDbfs) s << " | SP " << juce::String (*p.sampleDbfs, 6) << " dBFS @ " << juce::String (static_cast<double> (p.sampleFrame) / a.sampleRate, 9) << " s (frame " << juce::String (p.sampleFrame) << ")";
        else s << " | SP unavailable (silence; no location)";
        if (p.trueDbtp)
        {
            s << " | TP " << juce::String (*p.trueDbtp, 6) << " dBTP @ " << juce::String (p.truePeakEighthFrame / (8. * a.sampleRate), 9)
              << " s (frame " << juce::String (p.truePeakEighthFrame / 8., 3) << ")";
            if (p.truePeakEighthFrame < 0 || p.truePeakEighthFrame >= a.frames * 8) s << " [boundary interpolation]";
        }
        else s << " | TP unavailable (silence; no location)";
        if (static_cast<size_t> (row) < annotations.peakHeadroom.size() && annotations.peakHeadroom[static_cast<size_t> (row)])
            s << " | TP headroom " << juce::String (*annotations.peakHeadroom[static_cast<size_t> (row)], 6) << " dB";
    }
    return s;
}
void StandaloneDiagnosticsPanel::paintListBoxItem (int row, juce::Graphics& g, int width, int height, bool selected)
{
    if (selected) g.fillAll (juce::Colour (0xff384c66)); g.setColour (juce::Colours::white); g.setFont (12.0f);
    g.drawText (rowText (row), 5, 0, width - 10, height, juce::Justification::centredLeft, true);
}
void StandaloneDiagnosticsPanel::selectedRowsChanged (int row)
{
    if (! displayed || row < 0) return;
    details.setText (rowText (row) + "\nIntervals are [start,end) in audio-source time; plotted values sit at window ends. Overlapping-block shares are not independent passage energy or removal effects. Tail ties are retained; no LRA ceiling or processing recommendation.", false);
    if (const auto* w = windowForRow (row)) timeline.highlight (static_cast<double> (w->begin), static_cast<double> (w->end));
    else if (rankingBox.getSelectedId() == 5 && static_cast<size_t> (row) < displayed->analysis->peaks.size())
    {
        const auto& p = displayed->analysis->peaks[static_cast<size_t> (row)];
        if (p.trueDbtp) timeline.highlight (p.truePeakEighthFrame / 8., p.truePeakEighthFrame / 8.);
    }
}
void StandaloneDiagnosticsPanel::resized()
{
    auto b = getLocalBounds().reduced (12); auto row = b.removeFromTop (28);
    domainBox.setBounds (row.removeFromLeft (225)); row.removeFromLeft (12); rebuildButton.setBounds (row.removeFromLeft (160)); cancelButton.setBounds (row.removeFromLeft (160));
    status.setBounds (b.removeFromTop (25)); summary.setBounds (b.removeFromTop (105)); b.removeFromTop (6);
    timeline.setBounds (b.removeFromTop (std::max (210, b.getHeight() * 46 / 100))); b.removeFromTop (6);
    distributionBounds = b.removeFromTop (135); b.removeFromTop (6);
    rankingBox.setBounds (b.removeFromTop (28).removeFromLeft (400)); b.removeFromTop (4);
    details.setBounds (b.removeFromBottom (72)); b.removeFromBottom (4); ranking.setBounds (b);
}
void StandaloneDiagnosticsPanel::paintDistribution (juce::Graphics& g)
{
    g.setColour (juce::Colour (0xff202631)); g.fillRect (distributionBounds); g.setColour (juce::Colours::white); g.setFont (12.f);
    auto b = distributionBounds.reduced (8);
    g.drawText ("LRA DISTRIBUTION (loudness axis, not programme time). Green included / grey excluded; white CDF at 1 LU bin edges.", b.removeFromTop (18), juce::Justification::left);
    if (! displayed) return;
    const auto& a = *displayed->analysis;
    juce::String text = "Gate: ";
    text << (a.lraGate ? juce::String (std::max (-70., *a.lraGate), 6) + " LUFS" : "unavailable")
         << " | P10 " << (a.p10 ? juce::String (*a.p10, 6) : "unavailable") << " | P95 " << (a.p95 ? juce::String (*a.p95, 6) : "unavailable")
         << " | includes marked EOF drain windows; no extension of programme duration";
    g.drawText (text, b.removeFromTop (18), juce::Justification::left);
    if (a.histogram.empty()) { g.drawText ("No finite short-term population; silence/unavailable is not 0 LUFS.", b, juce::Justification::centred); return; }
    auto labels = b.removeFromBottom (18); auto graph = b;
    size_t maximum = 1; for (const auto& bin : a.histogram) maximum = std::max (maximum, bin.included + bin.excluded);
    const auto low = a.histogram.front().lower, high = a.histogram.back().lower + 1;
    const auto x = [&] (double value) { return static_cast<float> (graph.getX() + (value - low) * graph.getWidth() / (high - low)); };
    juce::Path cdf; cdf.startNewSubPath (x(low), static_cast<float> (graph.getBottom()));
    for (const auto& bin : a.histogram)
    {
        const auto width = std::max (1.f, x(bin.lower + 1) - x(bin.lower));
        const auto height = static_cast<float> (graph.getHeight()) * static_cast<float> (bin.included + bin.excluded) / static_cast<float> (maximum);
        g.setColour (juce::Colours::grey); g.fillRect (x(bin.lower), graph.getBottom() - height, width, height);
        const auto includedHeight = static_cast<float> (graph.getHeight()) * static_cast<float> (bin.included) / static_cast<float> (maximum);
        g.setColour (juce::Colours::lightgreen); g.fillRect (x(bin.lower), graph.getBottom() - includedHeight, width, includedHeight);
        if (! a.lraOrder.empty()) cdf.lineTo (x(bin.lower + 1), static_cast<float> (graph.getBottom()) - static_cast<float> (graph.getHeight()) * static_cast<float> (bin.cumulativeIncluded) / static_cast<float> (a.lraOrder.size()));
    }
    if (! a.lraOrder.empty()) { g.setColour (juce::Colours::white); g.strokePath (cdf, juce::PathStrokeType (1.2f)); }
    g.setColour (juce::Colours::yellow);
    for (const auto value : { a.p10, a.p95, a.lraGate ? std::optional<double> (std::max (-70., *a.lraGate)) : std::nullopt })
        if (value && *value >= low && *value <= high) g.drawVerticalLine (static_cast<int> (x(*value)), static_cast<float> (graph.getY()), static_cast<float> (graph.getBottom()));
    g.setColour (juce::Colours::white); g.drawText (juce::String (low, 0) + " LUFS", labels, juce::Justification::left);
    g.drawText (juce::String (high, 0) + " LUFS", labels, juce::Justification::right);
}
void StandaloneDiagnosticsPanel::paint (juce::Graphics& g) { g.fillAll (juce::Colour (0xff151b24)); paintDistribution (g); }
}
