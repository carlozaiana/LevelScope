#include "StandaloneSections.h"
#include <juce_cryptography/juce_cryptography.h>
#include <sstream>
#include <locale>

namespace levelscope::standalone
{
juce::String sectionInputIdentity (const ServiceContext& c)
{
    if (! c.isValid() || ! c.isFullProgram() || (c.media.sampleRate != 44100 && c.media.sampleRate != 48000
        && c.media.sampleRate != 96000 && c.media.sampleRate != 192000)) return {};
    std::ostringstream s; s.imbue (std::locale::classic()); const auto hash = c.media.sha256.toStdString();
    s << "section-input-v1;" << hash.size() << ':' << hash << ';' << static_cast<int> (c.media.sampleRate)
      << ';' << c.media.frames << ';' << c.media.channels << ';' << c.range.begin << ';' << c.range.end;
    for (int n = 0; n < c.fileLayout.count; ++n) s << ';' << static_cast<int> (c.fileLayout.roles[n]);
    const auto bytes = s.str(); return juce::SHA256 (bytes.data(), bytes.size()).toHexString();
}
SectionsController::~SectionsController() { cancel(); if (future.valid()) future.wait(); }
void SectionsController::select (std::shared_ptr<const MeasuredDiagnostics> next, const proposal::SectionMap& map)
{
    const auto key = next && next->analysis && next->snapshot && map.valid() ? next->analysis->identity + ";" + map.digest() : juce::String();
    if (next == selected && key == selectedKey) return;
    cancel(); selected = std::move (next); selectedMap = map; selectedKey = key; attempted = false;
    status = key.isNotEmpty() ? "Preparing section diagnostics..." : "No valid measured section diagnostics; measurements are unchanged.";
}
void SectionsController::cancel()
{
    if (control) control->cancel.store (true);
    ++generation; accepted.reset(); attempted = true; status = "Section diagnostics cancelled; section edits and audio are unchanged.";
}
void SectionsController::retry() { if (! busy()) { cancel(); attempted = false; } }
void SectionsController::poll()
{
    if (future.valid() && future.wait_for (std::chrono::seconds (0)) == std::future_status::ready)
    {
        try
        {
            auto done = future.get();
            if (runningGeneration == generation && control && ! control->isCancelled())
            { accepted = std::move (done); status = accepted ? "Measured section diagnostics ready. No audio change." : "Section diagnostics cancelled."; }
        }
        catch (const std::exception& e) { if (runningGeneration == generation) status = e.what(); }
        control.reset();
    }
    if (future.valid() || attempted || selectedKey.isEmpty() || ! selected) return;
    attempted = true; ++buildsStarted; runningGeneration = generation; control = std::make_shared<JobControl>();
    status = "Building section diagnostics...";
    future = std::async (std::launch::async, [d = selected, map = selectedMap, token = control]
    {
        if (! d->snapshot || ! d->analysis || ! d->measuredBinding.currentIdentityMatches
            || d->analysis->snapshotIdentity != d->measuredBinding.snapshotIdentity
            || d->snapshot->measurementVersion != measurement::version
            || (d->snapshot->domain != MeasurementDomain::source && d->snapshot->domain != MeasurementDomain::currentState)
            || d->snapshot->measuredContentHash != d->snapshot->context.media.sha256)
            throw std::invalid_argument ("Valid measured diagnostics required; predictions cannot populate sections");
        return proposal::analyzeSections (*d->analysis, map, sectionInputIdentity (d->snapshot->context), [token] { return token->isCancelled(); });
    });
}
}
