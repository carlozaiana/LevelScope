#include "StandaloneDiagnostics.h"
#include <juce_cryptography/juce_cryptography.h>

namespace levelscope::standalone
{
DiagnosticAnnotations annotateDiagnostics (const MeasuredDiagnostics& d, const TargetProfile& profile, proposal::Applicability applicability)
{
    auto input = d.measuredBinding; input.applicability = applicability;
    const auto evaluation = proposal::evaluateMeasuredProgramme (profile, input);
    DiagnosticAnnotations a;
    const auto key = d.analysis->identity + ";" + evaluation.identity;
    a.identity = juce::SHA256 (key.toRawUTF8(), key.getNumBytesAsUTF8()).toHexString();
    if (evaluation.verdict != proposal::Verdict::notEvaluated)
    {
        a.ceiling = profile.limits.truePeakCeilingDbtp;
        a.text = "Selected profile: " + profile.displayName + " | TP ceiling " + juce::String (*a.ceiling, 6) + " dBTP. Headroom = ceiling - measured channel maximum.";
    }
    else a.text = "Profile annotation unavailable: " + evaluation.message;
    for (const auto& peak : d.analysis->peaks)
        a.peakHeadroom.push_back (a.ceiling && peak.trueDbtp ? std::optional<double> (*a.ceiling - *peak.trueDbtp) : std::nullopt);
    a.text += " LRA descriptions add no compliance rule.";
    return a;
}
DiagnosticsController::~DiagnosticsController()
{
    cancel();
    if (future.valid()) future.wait();
}
void DiagnosticsController::select (std::shared_ptr<const AnalysisSnapshot> next, const ServiceContext& context, MeasurementDomain which)
{
    const auto identity = next ? OfflineAnalysisService::identityFor (context, which, context.media.sha256) : juce::String();
    if (next == selected && identity == expectedIdentity && which == domain) return;
    cancel(); ++generation; selected = std::move (next); expected = context; domain = which;
    expectedIdentity = identity; accepted.reset(); attempted = false;
    status = selected ? "Preparing measured diagnostics..." : "No valid measured snapshot selected; measure this state first.";
}
void DiagnosticsController::cancel()
{
    if (control) control->cancel.store (true);
    ++generation; attempted = true;
    status = "Diagnostics cancelled; measurements and processing unchanged. Use Build diagnostics to retry.";
}
void DiagnosticsController::retry() { if (! busy()) { attempted = false; accepted.reset(); } }
void DiagnosticsController::poll()
{
    if (future.valid() && future.wait_for (std::chrono::seconds (0)) == std::future_status::ready)
    {
        auto done = future.get();
        if (runningGeneration == generation && selected == running && control && ! control->isCancelled())
        {
            if (done.value) { accepted = std::make_shared<const MeasuredDiagnostics> (std::move (*done.value)); status = "Measured diagnostics ready. Inspection only; no audio changes."; }
            else status = done.message;
        }
        running.reset(); control.reset();
    }
    if (future.valid() || ! selected || attempted) return;
    attempted = true; ++buildsStarted; running = selected; runningGeneration = generation;
    control = std::make_shared<JobControl>(); status = "Building measured diagnostics...";
    future = std::async (std::launch::async, [snapshot = selected, context = expected, which = domain, token = control]
    {
        ServiceResult<MeasuredDiagnostics> result;
        try
        {
            if (token->isCancelled()) { result.error = ErrorCode::cancelled; return result; }
            if (which != MeasurementDomain::source && which != MeasurementDomain::currentState)
                throw std::invalid_argument ("Diagnostics require measured Source or Current State");
            auto binding = measuredProgrammeForSnapshot (snapshot.get(), context, which, proposal::Applicability::unconfirmed);
            if (! binding.currentIdentityMatches || binding.measurementVersion != measurement::version)
                throw std::invalid_argument ("Diagnostic snapshot is stale or unsupported; re-measure");
            auto analysis = proposal::analyzeProgramme (snapshot->detail, binding.snapshotIdentity, [token] { return token->isCancelled(); });
            if (! analysis) { result.error = ErrorCode::cancelled; result.message = "Diagnostics cancelled"; return result; }
            result.value = MeasuredDiagnostics { snapshot, analysis, binding }; result.error = ErrorCode::none;
        }
        catch (const std::exception& e) { result.error = ErrorCode::invalidRequest; result.message = e.what(); }
        return result;
    });
}
}
