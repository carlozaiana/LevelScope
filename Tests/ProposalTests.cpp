#include "Standalone/StandaloneSessionModel.h"
#include "Standalone/StandaloneAudition.h"
#include <cmath>
#include <iostream>
#include <limits>

using namespace levelscope;
using namespace levelscope::standalone;
using namespace levelscope::proposal;
namespace
{
juce::Array<juce::var> checks;
void expect (bool value, const char* name)
{
    auto* c = new juce::DynamicObject(); c->setProperty ("name", name); c->setProperty ("passed", value);
    checks.add (juce::var (c));
    if (! value) throw std::runtime_error (name);
}
TrimRequest input (double i = -26, double tp = -5, double gain = 0)
{
    TrimRequest r; r.profile = ebuR128ProgrammeQc(); r.currentGainDb = gain;
    r.source.domain = SnapshotDomain::source; r.current.domain = SnapshotDomain::currentState;
    for (auto* m : { &r.source, &r.current })
    {
        m->snapshotIdentity = m == &r.source ? "source-boundary-fixture" : "current-boundary-fixture";
        m->measurementVersion = measurement::version; m->currentIdentityMatches = true;
        m->applicability = Applicability::preProducedLinearProgramme;
        m->integrated = { i, measurement::Status::available };
        m->truePeak = { tp, measurement::Status::available };
    }
    r.sourceInputIdentity = r.currentInputIdentity = "common-input-fixture";
    return r;
}
void scalar()
{
    auto r = input(); auto p = draftGlobalTrim (r);
    expect (p.status == TrimStatus::draft && p.totalGainDb == 3 && p.deltaDb == 3, "increase reaches intended target");
    expect (p.predictedIntegrated == -23 && p.predictedTruePeak == -2, "scalar predictions explicitly distinct");
    expect (p.identity == draftGlobalTrim (r).identity && p.explanation == draftGlobalTrim (r).explanation, "deterministic repeated plan");
    expect (trimProposalText (p).contains ("prediction only"), "draft has no measured-pass wording");
    p = draftGlobalTrim (input (-20, 1));
    expect (p.status == TrimStatus::draft && p.totalGainDb == -3, "attenuation solves both selected rules");
    p = draftGlobalTrim (input (-29, -8, -6));
    expect (p.totalGainDb == 0 && p.deltaDb == 6, "delta is added once to measured Current State gain");
    p = draftGlobalTrim (input (-26, -3.9));
    expect (p.status == TrimStatus::draft && std::abs (*p.deltaDb - 2.9) < 1e-14
            && *p.predictedTruePeak <= -1 && *p.predictedIntegrated >= -23.2, "headroom constrained point inside QC allowance");
    p = draftGlobalTrim (input (-26, -3));
    expect (p.status == TrimStatus::noSolution && ! p.totalGainDb, "infeasible intersection gives no applicable plan");
    expect (draftGlobalTrim (input (-23.3, -9, 24)).status == TrimStatus::noSolution, "total upper gain bound");
    expect (draftGlobalTrim (input (-22.7, -9, -60)).status == TrimStatus::noSolution, "total lower gain bound");
    for (double i : { -23.2, -23.0, -22.8 })
    {
        p = draftGlobalTrim (input (i, -1, -6.125));
        expect (p.status == TrimStatus::noChange && p.totalGainDb == -6.125 && p.deltaDb == 0
                && ! p.predictedIntegrated, "already passing retains existing gain without renormalizing");
    }
    const auto inf = std::numeric_limits<double>::infinity();
    for (double edge : { -23.2, -22.8 })
    {
        p = draftGlobalTrim (input (std::nextafter (edge, -23), -9));
        expect (p.status == TrimStatus::noChange, "one ULP inside loudness allowance is no-op");
        p = draftGlobalTrim (input (std::nextafter (edge, edge < -23 ? -inf : inf), -9));
        expect (p.status == TrimStatus::draft, "one ULP outside loudness allowance requires proposal");
    }
    expect (draftGlobalTrim (input (-23, std::nextafter (-1., -inf))).status == TrimStatus::noChange, "one ULP below TP ceiling no-op");
    p = draftGlobalTrim (input (-23, std::nextafter (-1., inf)));
    expect (p.status == TrimStatus::draft && *p.deltaDb < 0 && *p.predictedTruePeak <= -1, "one ULP above TP ceiling requires strict correction");
    r = input (-23, -.7); p = draftGlobalTrim (r);
    expect (p.status == TrimStatus::noSolution, "meter tolerance is not extra peak headroom");
    r = input (-26.123456789, -10, -.123456789); p = draftGlobalTrim (r);
    expect (p.status == TrimStatus::draft && *p.totalGainDb == r.currentGainDb + *p.deltaDb
            && std::abs (*p.deltaDb - 3.123456789) < 1e-14, "gain retains full precision");
    r = input(); r.profile.limits.integratedLufs = -24; r.profile.definitionVersion = "test-target";
    expect (draftGlobalTrim (r).deltaDb == 2, "target is read from versioned profile");
    r = input(); r.current.lra = { {}, measurement::Status::insufficientDuration };
    expect (draftGlobalTrim (r).status == TrimStatus::draft, "unavailable informative LRA never blocks");
    r.profile.limits.hasMaxLoudnessRange = true; r.profile.limits.maxLoudnessRangeLu = 10;
    expect (draftGlobalTrim (r).status == TrimStatus::unavailable, "LRA steering profiles unsupported explicitly");
    for (int invalid = 0; invalid < 14; ++invalid)
    {
        r = input();
        switch (invalid)
        {
            case 0: r.currentInputIdentity = "different-content-or-roles-or-range"; break;
            case 1: r.sourceInputIdentity.clear(); break;
            case 2: r.current.currentIdentityMatches = false; break;
            case 3: r.source.currentIdentityMatches = false; break;
            case 4: r.current.domain = SnapshotDomain::predicted; break;
            case 5: r.source.domain = SnapshotDomain::currentState; break;
            case 6: r.current.applicability = Applicability::unconfirmed; break;
            case 7: r.profile = {}; break;
            case 8: r.processingVersion = "dynamics"; break;
            case 9: r.currentGainDb = inf; break;
            case 10: r.currentGainDb = 24.01; break;
            case 11: r.source.integrated = {}; break;
            case 12: r.current.truePeak = {}; break;
            case 13: r.current.measurementVersion = "unsupported"; break;
        }
        expect (draftGlobalTrim (r).status == TrimStatus::unavailable, "invalid missing stale unsupported or mismatched inputs rejected");
    }
    for (bool source : { false, true })
        for (bool peak : { false, true })
            for (double bad : { inf, -inf, std::numeric_limits<double>::quiet_NaN() })
            {
                r = input(); auto& m = source ? r.source : r.current;
                (peak ? m.truePeak : m.integrated).value = bad;
                expect (draftGlobalTrim (r).status == TrimStatus::unavailable, "non-finite metric cannot draft");
            }
    r = input(); const auto original = draftGlobalTrim (r).identity;
    r.profile.definitionVersion = "changed";
    expect (draftGlobalTrim (r).identity != original, "profile version binds draft identity");
    r = input(); r.current.snapshotIdentity = "new-measurement";
    expect (draftGlobalTrim (r).identity != original, "exact measured snapshot binds draft identity");
}

struct Fixture
{
    juce::File file = juce::File::getSpecialLocation (juce::File::tempDirectory).getNonexistentChildFile ("levelscope-m3", ".wav", false);
    Fixture (double first = -26, double second = -26, int channels = 2, int rate = 48000)
    {
        juce::WavAudioFormat f; auto stream = file.createOutputStream();
        expect (stream != nullptr, "fixture stream created");
        std::unique_ptr<juce::AudioFormatWriter> writer (f.createWriterFor (stream.get(), rate, static_cast<unsigned int> (channels), 24, {}, 0));
        expect (writer != nullptr, "fixture writer created"); stream.release();
        juce::AudioBuffer<float> audio (channels, rate * 8);
        for (int n = 0; n < audio.getNumSamples(); ++n)
            for (int c = 0; c < channels; ++c)
                audio.setSample (c, n, static_cast<float> (std::pow (10., (n < rate * 4 ? first : second) / 20)
                    * std::sin (2 * 3.141592653589793 * 1000 * n / rate)));
        expect (writer->writeFromAudioSampleBuffer (audio, 0, audio.getNumSamples()), "fixture written");
    }
    ~Fixture() { file.deleteFile(); }
};
void import (StandaloneSessionModel& s, const juce::File& file)
{
    JobControl control; auto v = validateMedia (file, control);
    expect (v.succeeded() && s.setValidatedSource (file, v.media), "source imported");
    s.setChannelRole (0, juce::AudioChannelSet::left); s.setChannelRole (1, juce::AudioChannelSet::right);
    expect (s.confirmLayout(), "roles confirmed");
    s.setSelectedTargetProfileIndex (1); s.setProgrammeApplicability (Applicability::preProducedLinearProgramme);
}
AnalysisSnapshot measure (StandaloneSessionModel& s, MeasurementDomain domain)
{
    expect (s.job.begin (JobStage::analysis, s.revision), "measurement begins");
    const auto request = s.measurementContext (domain); JobControl control;
    const auto result = OfflineAnalysisService().analyze (request, domain, control);
    expect (result.value.has_value(), "qualified worker returns actual measurement");
    expect (s.job.finish (request.ticket, s.revision, result.error, result.message), "measurement completes");
    expect (s.acceptMeasurement (*result.value), "measurement accepted");
    return *result.value;
}
TrimProposal draft (StandaloneSessionModel& s)
{
    const auto p = draftGlobalTrim (s.currentTrimRequest());
    expect (s.acceptTrimProposal (p), "fresh draft accepted");
    return p;
}
void lifecycle()
{
    Fixture f; StandaloneSessionModel s; import (s, f.file);
    measure (s, MeasurementDomain::source);
    expect (! s.getWorkflowReadiness().canDraftProposal(), "unmeasured Current State cannot draft");
    s.initializeCurrentStateFromSource(); measure (s, MeasurementDomain::currentState);
    expect (s.getWorkflowReadiness().canDraftProposal(), "measured matching inputs can draft despite failing rules");
    const auto sourceIdentity = s.evaluation (MeasurementDomain::source).identity;
    const auto currentIdentity = s.evaluation (MeasurementDomain::currentState).identity;
    const auto revision = s.revision; const auto job = s.job.ticket.id;
    const auto pointer = s.currentState.analysis.get();
    const auto p = draft (s);
    expect (p.status == TrimStatus::draft && s.canApplyTrim(), "apply enabled for current draft");
    expect (s.revision == revision && s.job.ticket.id == job && s.currentState.analysis.get() == pointer
            && s.evaluation (MeasurementDomain::currentState).identity == currentIdentity, "draft never mutates gain measurement or audio revision");
    s.referenceVideo.sha256 = "video"; s.referenceVideo.startOnAudio = video::Time (1, 1);
    s.sourceDocument.routing.fileToDevice[0] = 7; StandaloneAudition audition; audition.setMonitorGainDb (-40);
    expect (! audition.connected() && s.canApplyTrim() && draftGlobalTrim (s.currentTrimRequest()).identity == p.identity,
            "video offset routing monitor preserve draft and audio identity");
    s.discardTrim();
    expect (! s.canApplyTrim() && s.currentState.analysis.get() == pointer, "discard preserves measured audio");
    draft (s);
    s.trimDraft->totalGainDb = 24; // Application reconstructs from bound snapshots.
    expect (s.applyTrim() && s.currentStateDocument.staticGainDb == *p.totalGainDb, "explicit apply uses full-precision bound plan");
    const double applied = s.currentStateDocument.staticGainDb;
    expect (! s.currentState.hasMeasurement && s.currentStateNeedsRemeasurement()
            && s.evaluation (MeasurementDomain::currentState).verdict == Verdict::notEvaluated, "apply clears previous verdict and marks unverified");
    expect (s.evaluation (MeasurementDomain::source).identity == sourceIdentity, "apply preserves Source truth");
    expect (! s.applyTrim() && s.currentStateDocument.staticGainDb == applied, "double apply rejected without compounding");
    expect (! s.getWorkflowReadiness().canDraftProposal() && ! s.getWorkflowReadiness().canExport(), "prediction cannot draft again or enable export");
    measure (s, MeasurementDomain::currentState);
    expect (s.evaluation (MeasurementDomain::currentState).verdict == Verdict::pass, "actual processed measurement passes M2");
    const auto accepted = s.currentState.analysis.get(); const auto evaluated = s.evaluation (MeasurementDomain::currentState).identity;
    const auto beforeRevision = s.revision; const auto beforeJob = s.job.ticket.id;
    auto noop = draft (s);
    expect (noop.status == TrimStatus::noChange && ! s.canApplyTrim() && s.applyTrim(), "already passes gives no applicable audio change");
    expect (s.currentStateDocument.staticGainDb == applied && s.currentState.analysis.get() == accepted
            && s.evaluation (MeasurementDomain::currentState).identity == evaluated && s.revision == beforeRevision
            && s.job.ticket.id == beforeJob && ! s.currentStateNeedsRemeasurement(), "no-op preserves gain snapshot evaluation revisions and measurement count");
    s.setCurrentGainDb (-6); measure (s, MeasurementDomain::currentState);
    auto previous = draft (s); const auto measured = s.currentState.analysis.get(); const auto audioRevision = s.revision;
    s.setSelectedTargetProfileIndex (2);
    expect (! s.canApplyTrim() && ! s.acceptTrimProposal (previous)
            && s.currentState.analysis.get() == measured && s.revision == audioRevision, "profile change invalidates draft but reuses measurements");
    s.setSelectedTargetProfileIndex (1); previous = draft (s);
    s.setProgrammeApplicability (Applicability::unconfirmed);
    expect (! s.canApplyTrim() && ! s.acceptTrimProposal (previous), "applicability change rejects stale draft");
    s.setProgrammeApplicability (Applicability::preProducedLinearProgramme); previous = draft (s);
    s.setCurrentGainDb (-5);
    expect (! s.acceptTrimProposal (previous), "manual processing edit rejects stale draft");
    measure (s, MeasurementDomain::currentState); previous = draft (s);
    s.setChannelRole (0, juce::AudioChannelSet::right);
    expect (! s.acceptTrimProposal (previous), "role edit rejects draft");
    import (s, f.file);
    expect (! s.acceptTrimProposal (previous), "source replacement rejects draft");
    s.clearSourceFile();
    expect (! s.acceptTrimProposal (previous) && ! s.getWorkflowReadiness().canDraftProposal(), "source clear rejects draft");
}
void pairs()
{
    Fixture f; StandaloneSessionModel s; import (s, f.file);
    const auto source = measure (s, MeasurementDomain::source);
    s.initializeCurrentStateFromSource(); const auto current = measure (s, MeasurementDomain::currentState);
    for (int mismatch = 0; mismatch < 4; ++mismatch)
    {
        auto altered = current;
        if (mismatch == 0) altered.context.media.sha256 = altered.measuredContentHash = "other-content";
        if (mismatch == 1) std::swap (altered.context.fileLayout.roles[0], altered.context.fileLayout.roles[1]);
        if (mismatch == 2) { ++altered.context.media.frames; ++altered.context.range.end; ++altered.detail.frames; }
        if (mismatch == 3) { altered.context.media.sampleRate = 44100; altered.detail.sampleRate = 44100; }
        altered.identity = OfflineAnalysisService::identityFor (altered.context, altered.domain, altered.measuredContentHash);
        auto r = trimRequest (&source, &altered, source.context, altered.context, ebuR128ProgrammeQc(), Applicability::preProducedLinearProgramme);
        expect (r.source.currentIdentityMatches && r.current.currentIdentityMatches, "unit fixture snapshots individually valid against own expected context");
        expect (draftGlobalTrim (r).status == TrimStatus::unavailable, "individually valid but mismatched input pair rejected");
    }
    auto expected = current.context; expected.processing.staticGainDb = 1;
    auto r = trimRequest (&source, &current, source.context, expected, ebuR128ProgrammeQc(), Applicability::preProducedLinearProgramme);
    expect (! r.current.currentIdentityMatches && draftGlobalTrim (r).status == TrimStatus::unavailable, "Current State must match current processing plan");
    expected = current.context; expected.processing.graphJson = "{freeze:true}";
    r = trimRequest (&source, &current, source.context, expected, ebuR128ProgrammeQc(), Applicability::preProducedLinearProgramme);
    expect (draftGlobalTrim (r).status == TrimStatus::unavailable, "non-static graph or freeze rejected");
}
void unsuccessful()
{
    for (int kind = 0; kind < 3; ++kind)
    {
        Fixture f; StandaloneSessionModel s; import (s, f.file); measure (s, MeasurementDomain::source);
        s.initializeCurrentStateFromSource(); measure (s, MeasurementDomain::currentState); draft (s);
        expect (s.applyTrim(), "gain applied before unsuccessful verification");
        const auto gain = s.currentStateDocument.staticGainDb; const auto revision = s.revision;
        expect (s.job.begin (JobStage::analysis, s.revision), "verification job begins");
        auto context = s.measurementContext (MeasurementDomain::currentState); JobControl control;
        if (kind == 0) { s.job.requestCancel(); control.cancel.store (true); }
        if (kind == 1) f.file.deleteFile();
        if (kind == 2) f.file.appendText ("mutation");
        const auto result = OfflineAnalysisService().analyze (context, MeasurementDomain::currentState, control);
        expect (! result.value && result.error != ErrorCode::none, "actual worker cancellation or failure returned");
        expect (! s.job.finish (context.ticket, s.revision, result.error, result.message), "unsuccessful verification rejected");
        s.measurementUnsuccessful (MeasurementDomain::currentState, s.job.error);
        expect (s.currentStateDocument.isInitialized && s.currentStateDocument.staticGainDb == gain
                && s.currentStateNeedsRemeasurement() && ! s.currentState.analysis
                && s.evaluation (MeasurementDomain::currentState).verdict == Verdict::notEvaluated
                && s.currentStateDocument.statusText.contains ("UNVERIFIED"), "unsuccessful verification retains gain and no previous pass");
        expect (s.job.ticket.id == context.ticket.id && (kind == 2 || s.revision == revision)
                && ! s.canApplyTrim(), "no automatic retry correction or rollback");
    }
}
void gateDisagreement()
{
    // High/low blocks straddle the absolute gate at -50 dB Current State trim.
    // Normalizing brings both into the population; I+delta is only a prediction.
    Fixture f (-19, -21); StandaloneSessionModel s; import (s, f.file);
    measure (s, MeasurementDomain::source); s.initializeCurrentStateFromSource(); s.setCurrentGainDb (-50);
    measure (s, MeasurementDomain::currentState); const auto p = draft (s);
    expect (p.status == TrimStatus::draft && *p.predictedIntegrated >= -23.2 && *p.predictedIntegrated <= -22.8, "gate fixture scalar prediction inside allowance");
    expect (s.applyTrim(), "gate fixture explicitly applied");
    const auto gain = s.currentStateDocument.staticGainDb;
    measure (s, MeasurementDomain::currentState);
    expect (s.evaluation (MeasurementDomain::currentState).verdict == Verdict::fail
            && s.currentState.analysis && s.currentState.analysis->domain == MeasurementDomain::currentState,
            "actual gate-changing measurement overrides prediction with failure");
    expect (s.currentStateDocument.staticGainDb == gain && ! s.currentStateNeedsRemeasurement()
            && ! s.canApplyTrim(), "completed measured failure retains gain and stops without retry");
}

void ratesAndChannels()
{
    for (int rate : { 44100, 48000, 96000, 192000 })
    {
        const int channels = rate == 48000 ? 6 : rate == 44100 ? 1 : 2;
        Fixture f (-26, -26, channels, rate); StandaloneSessionModel s; JobControl control;
        const auto media = validateMedia (f.file, control);
        expect (media.succeeded() && s.setValidatedSource (f.file, media.media), "multirate source imports");
        const ChannelRole roles[] = { juce::AudioChannelSet::left, juce::AudioChannelSet::right,
            juce::AudioChannelSet::centre, juce::AudioChannelSet::LFE,
            juce::AudioChannelSet::leftSurround, juce::AudioChannelSet::rightSurround };
        for (int c = 0; c < channels; ++c) s.setChannelRole (c, channels == 1 ? juce::AudioChannelSet::centre : roles[c]);
        expect (s.confirmLayout(), "mono stereo or 5.1 roles confirmed");
        s.setSelectedTargetProfileIndex (1); s.setProgrammeApplicability (Applicability::preProducedLinearProgramme);
        measure (s, MeasurementDomain::source); s.initializeCurrentStateFromSource();
        const auto initial = measure (s, MeasurementDomain::currentState);
        const auto p = draft (s);
        expect (p.status == TrimStatus::draft && s.applyTrim(), "multirate trim applies");
        const auto actual = measure (s, MeasurementDomain::currentState);
        expect (s.evaluation (MeasurementDomain::currentState).verdict == Verdict::pass, "multirate measured result passes");
        for (int c = 0; c < channels; ++c)
        {
            const double observed = 20 * std::log10 (actual.detail.channelSamplePeaks[c] / initial.detail.channelSamplePeaks[c]);
            expect (std::abs (observed - *p.deltaDb) < 0.00001, "gain reaches every file channel including LFE");
        }
    }
}

juce::Array<juce::var> scalarOutcomes()
{
    juce::Array<juce::var> result;
    const double cases[][3] = { {-26, -5, 0}, {-20, 1, 0}, {-29, -8, -6}, {-26, -3.9, 0},
        {-26, -3, 0}, {-23.2, -1, -6.125}, {-22.8, -1, 2.125}, {-26.123456789, -10, -.123456789} };
    for (const auto& c : cases)
    {
        const auto p = draftGlobalTrim (input (c[0], c[1], c[2]));
        auto* row = new juce::DynamicObject(); row->setProperty ("identity", p.identity);
        row->setProperty ("status", static_cast<int> (p.status));
        row->setProperty ("total_gain", p.totalGainDb ? juce::var (*p.totalGainDb) : juce::var());
        row->setProperty ("delta", p.deltaDb ? juce::var (*p.deltaDb) : juce::var());
        row->setProperty ("predicted_I", p.predictedIntegrated ? juce::var (*p.predictedIntegrated) : juce::var());
        row->setProperty ("predicted_TP", p.predictedTruePeak ? juce::var (*p.predictedTruePeak) : juce::var());
        result.add (juce::var (row));
    }
    return result;
}
}
int main (int argc, char** argv)
{
    bool passed = false;
    try { scalar(); lifecycle(); pairs(); unsuccessful(); gateDisagreement(); ratesAndChannels(); passed = true; }
    catch (const std::exception& e) { std::cerr << "FAIL: " << e.what() << '\n'; }
    auto* r = new juce::DynamicObject();
    r->setProperty ("suite", "M3 global trim"); r->setProperty ("passed", passed); r->setProperty ("checks", checks);
    r->setProperty ("method", trimMethodVersion); r->setProperty ("profile_digest", ebuR128ProgrammeQc().digest());
    r->setProperty ("scalar_outcomes", scalarOutcomes());
    r->setProperty ("build", LEVELSCOPE_MEASUREMENT_BUILD); r->setProperty ("engine", measurement::version);
    r->setProperty ("os", juce::SystemStats::getOperatingSystemName());
#if defined(__aarch64__) || defined(__arm64__)
    r->setProperty ("architecture", "arm64");
#else
    r->setProperty ("architecture", "x86_64");
#endif
    if (argc == 3 && juce::String (argv[1]) == "--report")
        if (! juce::File (juce::String (argv[2])).replaceWithText (juce::JSON::toString (juce::var (r)))) return 2;
    std::cout << (passed ? "PASS: " : "FAIL: ") << checks.size() << " M3 checks; no audio device opened\n";
    return passed ? 0 : 1;
}
