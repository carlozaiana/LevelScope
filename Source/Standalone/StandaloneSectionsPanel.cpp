#include "StandaloneSectionsPanel.h"
#include <charconv>

namespace levelscope::standalone
{
namespace
{
juce::String intervalText (const proposal::ProgrammeSection& s, int rate)
{
    return "[" + juce::String (s.begin) + ", " + juce::String (s.end) + ") frames | ["
        + juce::String (static_cast<double> (s.begin) / rate, 6) + ", " + juce::String (static_cast<double> (s.end) / rate, 6) + ") s";
}
juce::String summaryText (const char* name, const proposal::SectionWindowSummary& s, bool gated)
{
    juce::String text; text << name << " sampled windows: " << juce::String (static_cast<int64_t> (s.contained.size())) << " contained, "
        << juce::String (static_cast<int64_t> (s.crossing.size())) << " crossing, " << juce::String (static_cast<int64_t> (s.padding.size())) << " EOF-padded. ";
    if (s.median) text << "Min / median / max " << juce::String (*s.minimum, 6) << " / " << juce::String (*s.median, 6) << " / " << juce::String (*s.maximum, 6) << " LUFS.";
    else text << (s.contained.size() == 0 ? "Unavailable: no wholly contained windows." : "Unavailable: no finite loudness (zero window energy).");
    text << " Missing " << juce::String (static_cast<int64_t> (s.missingCount)) << ".";
    if (gated) text << " Programme gate: " << juce::String (static_cast<int64_t> (s.included)) << " included / " << juce::String (static_cast<int64_t> (s.excluded)) << " excluded.";
    else text << " Ungated.";
    return text + "\n";
}
}
StandaloneSectionsPanel::StandaloneSectionsPanel (const StandaloneSessionModel& s) : session (s)
{
    domainBox.addItem ("Measured Source", 1); domainBox.addItem ("Measured Current State", 2); domainBox.setSelectedId (1, juce::dontSendNotification);
    domainBox.setName ("Section measured state"); domainBox.onChange = [this] { synchronize(); };
    windowKind.addItem ("Momentary windows", 1); windowKind.addItem ("Short-term windows", 2); windowKind.addItem ("LRA population windows", 3);
    windowScope.addItem ("Fully contained", 1); windowScope.addItem ("Boundary crossing", 2); windowScope.addItem ("EOF padding", 3);
    windowKind.setSelectedId (1, juce::dontSendNotification); windowScope.setSelectedId (1, juce::dontSendNotification);
    windowKind.setName ("Section window type"); windowScope.setName ("Section window membership");
    for (auto* b : { &windowKind, &windowScope }) b->onChange = [this] { windowList.deselectAllRows(); windowList.updateContent(); windowList.repaint(); inspection.clear(); };
    boundary.setName ("Section boundary source frame"); name.setName ("Section name");
    boundaryLabel.setText ("Boundary frame", juce::dontSendNotification); nameLabel.setText ("Name", juce::dontSendNotification);
    for (auto* e : { &summary, &inspection }) { e->setMultiLine (true); e->setReadOnly (true); e->setFont (13.f); }
    summary.setName ("Measured section summary"); inspection.setName ("Section window inspection"); status.setName ("Section work status");
    splitButton.onClick = [this] { performEdit (0); }; moveButton.onClick = [this] { performEdit (1); };
    mergeButton.onClick = [this] { performEdit (2); }; renameButton.onClick = [this] { performEdit (3); };
    rebuildButton.onClick = [this] { editNotice.clear(); diagnostics.retry(); sections.retry(); synchronize(); };
    cancelButton.onClick = [this] { editNotice.clear(); diagnostics.cancel(); sections.cancel(); synchronize(); };
    timeline.onInspect = [this] (const juce::String& text) { inspection.setText (text, false); };
    for (auto* c : std::array<juce::Component*, 20> { &domainBox, &windowKind, &windowScope, &splitButton, &moveButton, &mergeButton, &renameButton,
         &rebuildButton, &cancelButton, &boundary, &name, &summary, &inspection, &boundaryLabel, &nameLabel, &status, &timeline, &sectionList, &windowList, nullptr })
        if (c) addAndMakeVisible (*c);
    sectionList.setRowHeight (27); windowList.setRowHeight (25); setSize (1080, 850); startTimerHz (10);
}
StandaloneSectionsPanel::~StandaloneSectionsPanel() { stopTimer(); diagnostics.cancel(); sections.cancel(); }
void StandaloneSectionsPanel::suspend() { suspended = true; diagnostics.cancel(); sections.cancel(); synchronize(); }
void StandaloneSectionsPanel::resume() { suspended = false; diagnostics.retry(); sections.retry(); synchronize(); }
void StandaloneSectionsPanel::timerCallback()
{
    const auto context = session.measurementContext (MeasurementDomain::source);
    const auto input = session.sourceContentInvalid ? juce::String() : sectionInputIdentity (context);
    if (document.bind (input, input.isNotEmpty() ? static_cast<int> (context.media.sampleRate) : 0, input.isNotEmpty() ? context.media.frames : 0))
    { selectedId = 0; editNotice.clear(); }
    const auto key = document.map().valid() ? document.map().digest() : juce::String();
    const bool mapChanged = key != mapKey;
    if (mapChanged)
    {
        mapKey = key; refreshing = true; sectionList.updateContent(); sectionList.deselectAllRows();
        int selectedRow = 0;
        for (size_t n = 0; n < document.map().sections.size(); ++n) if (document.map().sections[n].id == selectedId) selectedRow = static_cast<int> (n);
        if (! document.map().sections.empty()) { selectedId = document.map().sections[static_cast<size_t> (selectedRow)].id; sectionList.selectRow (selectedRow); }
        refreshing = false; sectionList.repaint();
        std::vector<int64_t> boundaries; for (const auto& row : document.map().sections) boundaries.push_back (row.end);
        timeline.setSectionBoundaries (std::move (boundaries));
    }
    const auto domain = domainBox.getSelectedId() == 2 ? MeasurementDomain::currentState : MeasurementDomain::source;
    const auto& m = domain == MeasurementDomain::source ? session.source : session.currentState;
    auto snapshot = input.isNotEmpty() && m.hasMeasurement && (domain != MeasurementDomain::currentState || ! session.currentStateNeedsRemeasurement()) ? m.analysis : nullptr;
    if (! suspended) diagnostics.select (snapshot, session.measurementContext (domain), domain);
    diagnostics.poll(); const auto nextMeasured = diagnostics.result();
    const bool measuredChanged = nextMeasured != measured;
    if (measuredChanged) { measured = nextMeasured; timeline.setAnalysis (measured ? measured->analysis : nullptr); annotationKey.clear(); }
    if (! suspended) sections.select (measured, document.map());
    sections.poll(); const auto next = sections.result(); const bool resultChanged = next != displayed; displayed = next;
    const auto newAnnotation = session.getSelectedTargetProfile().digest() + ";" + juce::String (static_cast<int> (session.programmeApplicability));
    if (mapChanged || measuredChanged || resultChanged || newAnnotation != annotationKey) { annotationKey = newAnnotation; refreshDetails(); }
    const bool busy = diagnostics.busy() || sections.busy(); cancelButton.setEnabled (busy);
    rebuildButton.setEnabled (snapshot != nullptr && ! busy);
    status.setText (editNotice.isNotEmpty() ? editNotice : diagnostics.busy() ? diagnostics.status : sections.status, juce::dontSendNotification);
}
const proposal::ProgrammeSection* StandaloneSectionsPanel::selectedSection() const
{ for (const auto& s : document.map().sections) if (s.id == selectedId) return &s; return nullptr; }
const proposal::SectionSummary* StandaloneSectionsPanel::selectedSummary() const
{ if (displayed) for (const auto& s : displayed->sections) if (s.section.id == selectedId) return &s; return nullptr; }
const proposal::SectionWindowSummary* StandaloneSectionsPanel::windowSummary() const
{ if (const auto* s = selectedSummary()) return windowKind.getSelectedId() == 1 ? &s->momentary : windowKind.getSelectedId() == 2 ? &s->shortTerm : &s->lraWindows; return nullptr; }
const proposal::WindowSubset* StandaloneSectionsPanel::subset() const
{ if (const auto* s = windowSummary()) return windowScope.getSelectedId() == 1 ? &s->contained : windowScope.getSelectedId() == 2 ? &s->crossing : &s->padding; return nullptr; }
const proposal::DiagnosticWindow* StandaloneSectionsPanel::windowAt (int row) const
{
    const auto* selected = subset(); if (! measured || ! selected || row < 0) return nullptr;
    const auto n = selected->index (static_cast<size_t> (row)); if (! n) return nullptr;
    const auto& a = *measured->analysis; const auto& windows = windowKind.getSelectedId() == 1 ? a.momentary : windowKind.getSelectedId() == 2 ? a.shortTerm : a.lraWindows;
    return *n < windows.size() ? &windows[*n] : nullptr;
}
void StandaloneSectionsPanel::performEdit (int operation)
{
    const auto* row = selectedSection(); if (! row) return;
    int64_t frame = 0;
    if (operation < 2)
    {
        const auto bytes = boundary.getText().trim().toStdString();
        const auto parsed = std::from_chars (bytes.data(), bytes.data() + bytes.size(), frame);
        if (bytes.empty() || parsed.ec != std::errc() || parsed.ptr != bytes.data() + bytes.size())
        { editNotice = "Enter a whole audio-source frame number; fractional or rounded seconds are not accepted."; synchronize(); return; }
    }
    juce::String error;
    const bool ok = operation == 0 ? document.split (row->id, frame, error) : operation == 1 ? document.moveEnd (row->id, frame, error)
        : operation == 2 ? document.mergeNext (row->id, error) : document.rename (row->id, name.getText(), error);
    editNotice = ok ? "Section edit accepted. Audio, measurements, verdicts and proposals unchanged." : error;
    synchronize();
}
void StandaloneSectionsPanel::refreshDetails()
{
    windowList.deselectAllRows(); windowList.updateContent(); windowList.repaint(); inspection.clear();
    const auto* row = selectedSection();
    for (auto* b : { &splitButton, &moveButton, &mergeButton, &renameButton }) b->setEnabled (row != nullptr);
    boundary.setEnabled (row != nullptr); name.setEnabled (row != nullptr);
    if (! row) { boundary.clear(); name.clear(); summary.setText ("Import a source and confirm its channel roles. Sections are in-memory annotations; they do not change audio.", false); timeline.clearHighlight(); return; }
    boundary.setText (juce::String (row->end), false); name.setText (row->name, false);
    timeline.highlight (static_cast<double> (row->begin), static_cast<double> (row->end), false);
    juce::String text = row->name + " | " + intervalText (*row, document.map().sampleRate) + "\n";
    const auto* a = selectedSummary();
    if (! a || ! measured) text += "No current measured section result. Select and measure this state; predictions cannot populate this view.\n";
    else
    {
        text << (measured->snapshot->domain == MeasurementDomain::source ? "SOURCE" : "CURRENT STATE") << " | Full-resolution, wholly contained windows only.\n";
        text << summaryText ("M", a->momentary, true) << summaryText ("S", a->shortTerm, false) << summaryText ("LRA population", a->lraWindows, true);
        text << "Programme I " << diagnosticMetricText (measured->analysis->integrated, "LUFS") << " | LRA " << diagnosticMetricText (measured->analysis->lra, "LU") << ".\n";
        text << "Programme included block-energy share: " << (a->momentary.programmeBlockEnergyShare ? juce::String (100 * *a->momentary.programmeBlockEnergyShare, 6) + "%" : "unavailable (no eligible programme blocks)")
             << " (overlapping blocks, not independent section energy).\n";
        text << "Programme LRA-tail membership: " << juce::String (static_cast<int64_t> (a->lraWindows.lowerTail)) << " lower / " << juce::String (static_cast<int64_t> (a->lraWindows.upperTail)) << " upper. " << diagnosticLraInterpretation (*measured->analysis) << "\n";
        if (row->end - row->begin < static_cast<int64_t> (document.map().sampleRate) * 3) text << "Short section: S/LRA summaries may have no complete windows.\n";
        for (const auto& located : a->programmeMaxima)
        {
            const auto& peak = measured->analysis->peaks[located.peakIndex];
            text << "Programme maximum located here: file channel " << juce::String (peak.channel + 1);
            if (located.sample) text << " SP " << juce::String (*peak.sampleDbfs, 6) << " dBFS @ frame " << juce::String (peak.sampleFrame);
            if (located.truePeak) text << " TP " << juce::String (*peak.trueDbtp, 6) << " dBTP @ eighth-frame " << juce::String (peak.truePeakEighthFrame);
            text << "\n";
        }
        text << "Local peak maxima are not measured; absence of a programme maximum does not establish peak safety.\n";
        text << annotateDiagnostics (*measured, session.getSelectedTargetProfile(), session.programmeApplicability).text << "\n";
    }
    text << "No section I/LRA or section compliance verdict. Crossing/padded windows are excluded from summaries and inspectable below. Section edits never seek, play or process audio.";
    summary.setText (text, false);
}
int StandaloneSectionsPanel::getNumRows() { return static_cast<int> (document.map().sections.size()); }
void StandaloneSectionsPanel::paintListBoxItem (int n, juce::Graphics& g, int w, int h, bool selected)
{
    if (n < 0 || n >= getNumRows()) return;
    if (selected) g.fillAll (juce::Colour (0xff384c66)); g.setColour (juce::Colours::white); g.setFont (12.f);
    const auto& s = document.map().sections[static_cast<size_t> (n)];
    g.drawText (juce::String (n + 1) + ". " + s.name + " | " + intervalText (s, document.map().sampleRate), 5, 0, w - 10, h, juce::Justification::centredLeft, true);
}
void StandaloneSectionsPanel::selectedRowsChanged (int row)
{
    if (refreshing || row < 0 || row >= getNumRows()) return;
    selectedId = document.map().sections[static_cast<size_t> (row)].id; refreshDetails();
    if (const auto* s = selectedSection()) timeline.highlight (static_cast<double> (s->begin), static_cast<double> (s->end));
}
int StandaloneSectionsPanel::Windows::getNumRows() { const auto* s = panel.subset(); return s ? static_cast<int> (s->size()) : 0; }
void StandaloneSectionsPanel::Windows::paintListBoxItem (int row, juce::Graphics& g, int w, int h, bool selected)
{
    if (const auto* window = panel.windowAt (row))
    {
        if (selected) g.fillAll (juce::Colour (0xff384c66)); g.setColour (juce::Colours::white); g.setFont (12.f);
        g.drawText (diagnosticWindowText (*window, panel.document.map().sampleRate, panel.document.map().frames, panel.windowKind.getSelectedId() == 2), 5, 0, w - 10, h, juce::Justification::centredLeft, true);
    }
}
void StandaloneSectionsPanel::Windows::selectedRowsChanged (int row)
{
    if (const auto* w = panel.windowAt (row))
    {
        panel.inspection.setText (diagnosticWindowText (*w, panel.document.map().sampleRate, panel.document.map().frames, panel.windowKind.getSelectedId() == 2)
            + "\nOriginal programme gate membership; window values sit at their ends. Overlap counts are not independent durations.", false);
        panel.timeline.highlight (static_cast<double> (w->begin), static_cast<double> (w->end));
    }
}
void StandaloneSectionsPanel::resized()
{
    auto b = getLocalBounds().reduced (12); auto top = b.removeFromTop (28);
    domainBox.setBounds (top.removeFromLeft (225)); rebuildButton.setBounds (top.removeFromLeft (145)); cancelButton.setBounds (top.removeFromLeft (145));
    status.setBounds (b.removeFromTop (27)); auto edit = b.removeFromTop (28);
    boundaryLabel.setBounds (edit.removeFromLeft (105)); boundary.setBounds (edit.removeFromLeft (140));
    splitButton.setBounds (edit.removeFromLeft (120)); moveButton.setBounds (edit.removeFromLeft (100)); mergeButton.setBounds (edit.removeFromLeft (100));
    nameLabel.setBounds (edit.removeFromLeft (45)); renameButton.setBounds (edit.removeFromRight (90)); name.setBounds (edit);
    b.removeFromTop (6); sectionList.setBounds (b.removeFromTop (100)); b.removeFromTop (6);
    timeline.setBounds (b.removeFromTop (240)); b.removeFromTop (6);
    summary.setBounds (b.removeFromTop (std::max (160, b.getHeight() * 54 / 100))); b.removeFromTop (6);
    auto choose = b.removeFromTop (28); windowKind.setBounds (choose.removeFromLeft (240)); windowScope.setBounds (choose.removeFromLeft (200));
    inspection.setBounds (b.removeFromBottom (55)); b.removeFromBottom (4); windowList.setBounds (b);
}
void StandaloneSectionsPanel::paint (juce::Graphics& g) { g.fillAll (juce::Colour (0xff151b24)); }
}
