#include "Standalone/StandaloneSessionModel.h"
#include "Standalone/StandaloneAudition.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <type_traits>

using namespace levelscope;
using namespace levelscope::standalone;
using namespace levelscope::proposal;
namespace
{
juce::Array<juce::var> records;
void expect (bool condition, const char* name)
{
    auto* item = new juce::DynamicObject(); item->setProperty ("name", name); item->setProperty ("passed", condition);
    records.add (juce::var (item));
    if (! condition) throw std::runtime_error (name);
}
MeasuredProgramme measured()
{
    MeasuredProgramme m;
    m.domain = SnapshotDomain::source; m.currentIdentityMatches = true;
    m.snapshotIdentity = "synthetic-boundary-fixture"; m.measurementVersion = measurement::version;
    m.applicability = Applicability::preProducedLinearProgramme;
    m.integrated = { -23.0, measurement::Status::available };
    m.truePeak = { -1.0, measurement::Status::available };
    return m;
}
void rules()
{
    const auto p = ebuR128ProgrammeQc();
    const auto base = measured();
    const auto eval = [&] (MeasuredProgramme m) { return evaluateMeasuredProgramme (p, m); };
    expect (p.validationError().isEmpty(), "built-in profile valid");
    expect (p.limits.integratedLufs == -23 && p.limits.integratedToleranceLu == .2
            && p.limits.truePeakCeilingDbtp == -1 && ! p.limits.hasMaxLoudnessRange, "pinned R128 scalar rules");
    expect (p.truePeakMeterToleranceDb == .3 && p.truePeakMeterBandwidthHz == 20000
            && p.meterToleranceNote.isNotEmpty(), "documentary meter tolerance has explicit bandwidth");
    expect (eval (base).verdict == Verdict::pass && eval (base).rules.size() == 2, "missing informative LRA allows pass");
    expect (eval (base).message == "Passes selected programme loudness and true-peak rules.", "narrow success wording");
    expect (eval (base).identity == eval (base).identity, "deterministic evaluation");
    const double infinity = std::numeric_limits<double>::infinity();
    for (double edge : { -23.2, -22.8 })
    {
        auto m = base; m.integrated.value = edge;
        expect (eval (m).verdict == Verdict::pass, "integrated exact boundary inclusive");
        m.integrated.value = std::nextafter (edge, -23.0);
        expect (eval (m).verdict == Verdict::pass, "integrated one ULP inside");
        m.integrated.value = std::nextafter (edge, edge < -23 ? -infinity : infinity);
        expect (eval (m).verdict == Verdict::fail, "integrated one ULP outside with no epsilon");
        expect (juce::String (*m.integrated.value, 6) == juce::String (edge, 6), "display-equivalent integrated values retain distinct verdicts");
    }
    for (double value : { -1.0, std::nextafter (-1.0, -infinity), std::nextafter (-1.0, infinity), -0.7 })
    {
        auto m = base; m.truePeak.value = value;
        expect (eval (m).verdict == (value <= -1 ? Verdict::pass : Verdict::fail), "strict peak boundary including meter-tolerance trap");
    }
    auto m = base; m.integrated.value = -22.75;
    expect (eval (m).verdict == Verdict::fail, "engine test tolerance cannot extend QC allowance");
    for (bool peak : { false, true })
    {
        for (double invalid : { infinity, -infinity, std::numeric_limits<double>::quiet_NaN() })
        {
            m = base; (peak ? m.truePeak : m.integrated).value = invalid;
            expect (eval (m).verdict == Verdict::indeterminate, "non-finite required metric cannot pass");
        }
        for (auto status : { measurement::Status::silence, measurement::Status::belowGate, measurement::Status::insufficientDuration })
        {
            m = base; (peak ? m.truePeak : m.integrated) = { {}, status };
            expect (eval (m).verdict == Verdict::indeterminate, "missing required metric cannot pass");
        }
        m = base; (peak ? m.truePeak : m.integrated).status = measurement::Status::silence;
        expect (eval (m).verdict == Verdict::indeterminate, "inconsistent metric status cannot pass");
    }
    m = base; m.integrated = {}; m.truePeak.value = 0;
    expect (eval (m).verdict == Verdict::fail && eval (m).rules[0].verdict == Verdict::indeterminate,
            "known failure takes precedence without hiding indeterminate rule");
    m = base; m.lra = { std::numeric_limits<double>::quiet_NaN(), measurement::Status::insufficientDuration };
    expect (eval (m).verdict == Verdict::pass, "unused non-finite LRA is not a required rule");
    m = base; m.domain = SnapshotDomain::predicted;
    expect (eval (m).verdict == Verdict::notEvaluated, "predictions cannot satisfy measured rules");
    m.domain = SnapshotDomain::renderedOutput;
    expect (eval (m).verdict == Verdict::notEvaluated, "M2 cannot certify an encoded export");
    m = base; m.domain = SnapshotDomain::currentState;
    expect (eval (m).domain == SnapshotDomain::currentState && eval (m).identity != eval (base).identity, "domain is explicit and identity bound");
    m = base; m.currentIdentityMatches = false;
    expect (eval (m).verdict == Verdict::notEvaluated, "stale snapshot cannot pass");
    m = base; m.snapshotIdentity.clear();
    expect (eval (m).verdict == Verdict::notEvaluated, "missing snapshot binding cannot pass");
    m = base; m.measurementVersion = "plugin-history";
    expect (eval (m).verdict == Verdict::notEvaluated, "unsupported measurement semantics refused");
    for (auto applicability : { Applicability::unconfirmed, Applicability::liveProgramme,
            Applicability::deliberatelyLowerTarget, Applicability::shortForm, Applicability::distributionSpecific })
    {
        m = base; m.applicability = applicability;
        expect (eval (m).verdict == Verdict::notEvaluated, "unsupported or unconfirmed programme context refused");
    }
    m = base; m.snapshotIdentity = "different-snapshot";
    expect (eval (m).identity != eval (base).identity, "snapshot change changes evaluation binding");
    m = base; m.truePeak.value = -2;
    expect (eval (m).identity != eval (base).identity, "actual values participate in evaluation binding");

    auto optional = p; optional.id = "test-only-lra"; optional.limits.hasMaxLoudnessRange = true;
    optional.limits.maxLoudnessRangeLu = 10; optional.lraProvenance = "Synthetic test constraint, not an EBU default";
    expect (evaluateMeasuredProgramme (optional, base).verdict == Verdict::indeterminate, "required optional LRA cannot be missing");
    m = base; m.lra = { 10, measurement::Status::available }; m.lraStable = true;
    expect (evaluateMeasuredProgramme (optional, m).verdict == Verdict::pass, "explicit LRA ceiling inclusive");
    m.lra.value = std::nextafter (10.0, infinity);
    expect (evaluateMeasuredProgramme (optional, m).verdict == Verdict::fail, "explicit LRA ceiling one ULP outside");
    m.lra.value = 9; m.lraStable = false;
    expect (evaluateMeasuredProgramme (optional, m).verdict == Verdict::indeterminate, "required unstable LRA refused");
}
void profiles()
{
    const auto base = ebuR128ProgrammeQc();
    for (int index : { 0, 2, 3 })
        expect (evaluateMeasuredProgramme (StandaloneTargetProfileCatalog::getProfile (index), measured()).verdict == Verdict::notEvaluated,
                "placeholder profile never passes");
    auto p = base; p.definitionVersion = "test-v2";
    expect (p.digest() != base.digest(), "definition version changes profile identity");
    p = base; p.standardEdition = "test-edition";
    expect (p.digest() != base.digest(), "standard edition changes identity");
    p = base; p.integratedProvenance += " changed";
    expect (p.digest() != base.digest(), "provenance changes identity");
    p = base; p.limits.integratedToleranceLu = .1;
    expect (p.digest() != base.digest(), "limits change identity");
    p = base; p.applicability = Applicability::liveProgramme;
    expect (p.digest() != base.digest() && ! p.validationError().isEmpty(), "unsupported profile applicability refused");
    for (int invalid = 0; invalid < 12; ++invalid)
    {
        p = base;
        switch (invalid)
        {
            case 0: p.id.clear(); break;
            case 1: p.definitionVersion.clear(); break;
            case 2: p.standardEdition.clear(); break;
            case 3: p.referenceUrl.clear(); break;
            case 4: p.integratedProvenance.clear(); break;
            case 5: p.truePeakProvenance.clear(); break;
            case 6: p.method = MeasurementMethod::unsupported; break;
            case 7: p.limits.integratedToleranceLu = -1; break;
            case 8: p.limits.truePeakCeilingDbtp = std::numeric_limits<double>::quiet_NaN(); break;
            case 9: p.limits.integratedLufs = std::numeric_limits<double>::infinity(); break;
            case 10: p.limits.hasTruePeakCeiling = false; break;
            case 11: p.schemaVersion = 999; break;
        }
        expect (! p.validationError().isEmpty() && evaluateMeasuredProgramme (p, measured()).verdict == Verdict::notEvaluated,
                "invalid profile rejected before evaluation");
    }
    static_assert (! std::is_convertible<ProposalDraft, MeasuredProgramme>::value, "Predictions require an explicit measured adapter");
    static_assert (! std::is_convertible<ProposalDraft, AnalysisSnapshot>::value, "Predictions cannot enter the snapshot API");
}

AnalysisSnapshot measureSession (StandaloneSessionModel& s, MeasurementDomain domain)
{
    expect (s.job.begin (JobStage::analysis, s.revision), "measurement job begins");
    const auto request = s.measurementContext (domain); JobControl control;
    const auto result = OfflineAnalysisService().analyze (request, domain, control);
    expect (result.value.has_value(), "actual snapshot measurement succeeds");
    expect (s.job.finish (request.ticket, s.revision, result.error, result.message), "measurement completes");
    expect (s.acceptMeasurement (*result.value), "session accepts current worker result");
    return *result.value;
}
void lifecycle()
{
    auto file = juce::File::getSpecialLocation (juce::File::tempDirectory).getNonexistentChildFile ("levelscope-m2", ".wav", false);
    struct Cleanup { juce::File file; ~Cleanup() { file.deleteFile(); } } cleanup { file };
    {
        juce::WavAudioFormat format; auto stream = file.createOutputStream();
        expect (stream != nullptr, "M2 temporary fixture created");
        std::unique_ptr<juce::AudioFormatWriter> writer (format.createWriterFor (stream.get(), 48000, 2, 24, {}, 0));
        expect (writer != nullptr, "M2 PCM writer created"); stream.release();
        juce::AudioBuffer<float> audio (2, 48000 * 4);
        for (int i = 0; i < audio.getNumSamples(); ++i)
            for (int c = 0; c < 2; ++c) audio.setSample (c, i, static_cast<float> (std::pow (10., -23. / 20) * std::sin (2 * 3.141592653589793 * 1000 * i / 48000)));
        expect (writer->writeFromAudioSampleBuffer (audio, 0, audio.getNumSamples()), "M2 fixture encoded");
    }
    JobControl control; const auto validation = validateMedia (file, control);
    expect (validation.succeeded(), "M2 source validates");
    StandaloneSessionModel session;
    expect (session.setValidatedSource (file, validation.media), "M2 source imports");
    session.setChannelRole (0, juce::AudioChannelSet::left); session.setChannelRole (1, juce::AudioChannelSet::right);
    expect (session.confirmLayout(), "roles explicitly confirmed");
    session.setSelectedTargetProfileIndex (1);
    expect (session.measurementContext (MeasurementDomain::source).profileVersion == "1.0.0", "profile version populated in request");
    const auto source = measureSession (session, MeasurementDomain::source);
    expect (session.evaluation (MeasurementDomain::source).verdict == Verdict::notEvaluated, "unconfirmed programme context cannot pass");
    session.setProgrammeApplicability (Applicability::preProducedLinearProgramme);
    const auto initial = session.evaluation (MeasurementDomain::source);
    expect (initial.verdict == Verdict::pass, "real Source measurements pass selected rules");
    expect (initial.snapshotIdentity.isNotEmpty() && initial.profileDigest == ebuR128ProgrammeQc().digest(), "exact snapshot and profile bound");
    expect (evaluationText (initial).startsWith ("Source evaluation: PASS"), "Source display domain identified");
    const auto pointer = session.source.analysis.get(); const auto jobId = session.job.ticket.id;
    const auto audioRevision = session.revision;
    for (int profile : { 0, 2, 3, 1 })
    {
        session.setSelectedTargetProfileIndex (profile);
        expect (session.source.analysis.get() == pointer && session.job.ticket.id == jobId && session.revision == audioRevision,
                "profile selection reuses measurement with no new audio job");
        expect (session.evaluation (MeasurementDomain::source).verdict == (profile == 1 ? Verdict::pass : Verdict::notEvaluated),
                "profile selection reevaluates immediately");
    }
    expect (session.evaluation (MeasurementDomain::source).identity == initial.identity, "returning to exact profile restores deterministic evaluation");
    auto compatible = ebuR128ProgrammeQc(); compatible.id = "test-compatible"; compatible.definitionVersion = "2";
    compatible.limits.integratedLufs = -24;
    const auto changedProfile = evaluateSnapshot (pointer, session.measurementContext (MeasurementDomain::source), MeasurementDomain::source,
                                                  compatible, Applicability::preProducedLinearProgramme);
    expect (changedProfile.verdict == Verdict::fail && changedProfile.snapshotIdentity == initial.snapshotIdentity
            && changedProfile.identity != initial.identity, "compatible different profile reuses exact snapshot and changes verdict binding");
    session.initializeCurrentStateFromSource();
    const auto current = measureSession (session, MeasurementDomain::currentState);
    const auto currentEvaluation = session.evaluation (MeasurementDomain::currentState);
    expect (currentEvaluation.verdict == Verdict::pass && current.measured.integratedLufs == source.measured.integratedLufs, "0 dB measured Current State parity");
    expect (evaluationText (currentEvaluation).startsWith ("Current State evaluation: PASS"), "Current State display domain identified");
    session.referenceVideo.file = file; session.referenceVideo.sha256 = "video-only";
    session.referenceVideo.startOnAudio = video::Time (-1001, 30000); session.referenceVideo.streamId = 2;
    session.sourceDocument.routing.fileToDevice[0] = 7;
    StandaloneAudition audition; audition.setMonitorGainDb (-12);
    expect (! audition.connected() && session.evaluation (MeasurementDomain::source).identity == initial.identity
            && session.evaluation (MeasurementDomain::currentState).identity == currentEvaluation.identity,
            "video offset/stream, routing and monitor gain preserve measured evaluation identity with device closed");
    expect (session.setCurrentGainDb (-6), "gain edit accepted");
    expect (session.evaluation (MeasurementDomain::currentState).verdict == Verdict::notEvaluated && ! session.currentState.analysis,
            "gain edit removes Current State measurement and evaluation");
    expect (session.evaluation (MeasurementDomain::source).identity == initial.identity, "processing edit preserves Source evaluation");
    expect (! session.acceptMeasurement (current), "stale worker result rejected after gain edit");
    measureSession (session, MeasurementDomain::currentState);
    expect (session.evaluation (MeasurementDomain::currentState).verdict == Verdict::fail, "remeasured gain change evaluates actual samples");
    expect (! session.getWorkflowReadiness().canDraftProposal() && ! session.getWorkflowReadiness().canExport(), "evaluation does not enable proposal/export");
    for (int invalid = 0; invalid < 6; ++invalid)
    {
        auto context = session.measurementContext (MeasurementDomain::source);
        auto snapshot = source;
        switch (invalid)
        {
            case 0: context.sessionId = "other-session"; break;
            case 1: context.media.sha256 = "different-content"; break;
            case 2: context.fileLayout.roles[0] = juce::AudioChannelSet::centre; break;
            case 3: context.range.end--; break;
            case 4: snapshot.measured.truePeakDbtp = 0; break;
            case 5: snapshot.measurementVersion = "unsupported"; break;
        }
        expect (evaluateSnapshot (&snapshot, context, MeasurementDomain::source, ebuR128ProgrammeQc(), Applicability::preProducedLinearProgramme).verdict != Verdict::pass,
                "incompatible snapshot/request rejected");
    }
    expect (session.job.begin (JobStage::analysis, session.revision), "cancelled measurement job begins");
    auto pending = source; pending.context = session.measurementContext (MeasurementDomain::source);
    session.job.requestCancel();
    expect (! session.job.finish (pending.context.ticket, session.revision, ErrorCode::none, "done") && ! session.acceptMeasurement (pending),
            "cancelled result cannot replace accepted snapshot");
    expect (session.evaluation (MeasurementDomain::source).identity == initial.identity, "cancel retains previous valid measurement");
    session.setChannelRole (0, juce::AudioChannelSet::right);
    expect (session.evaluation (MeasurementDomain::source).verdict == Verdict::notEvaluated
            && session.evaluation (MeasurementDomain::currentState).verdict == Verdict::notEvaluated, "role edit invalidates both domains");
    expect (session.setValidatedSource (file, validation.media), "replacement accepted");
    expect (! session.source.analysis && ! session.currentState.analysis
            && session.programmeApplicability == Applicability::unconfirmed, "replacement clears snapshots and applicability confirmation");
    expect (! session.acceptMeasurement (source), "old source result rejected after replacement");
    session.clearSourceFile();
    expect (session.evaluation (MeasurementDomain::source).verdict == Verdict::notEvaluated, "clear source removes evaluation");
}
}
int main (int argc, char** argv)
{
    bool passed = false;
    try { profiles(); rules(); lifecycle(); passed = true; }
    catch (const std::exception& e) { std::cerr << "FAIL: " << e.what() << '\n'; }
    auto* report = new juce::DynamicObject();
    report->setProperty ("suite", "M2 measured compliance"); report->setProperty ("passed", passed);
    report->setProperty ("checks", records); report->setProperty ("profile_digest", ebuR128ProgrammeQc().digest());
    report->setProperty ("engine", measurement::version); report->setProperty ("evaluator", evaluationVersion);
    report->setProperty ("os", juce::SystemStats::getOperatingSystemName());
    report->setProperty ("build", LEVELSCOPE_MEASUREMENT_BUILD);
#if defined(__aarch64__) || defined(__arm64__)
    report->setProperty ("architecture", "arm64");
#else
    report->setProperty ("architecture", "x86_64");
#endif
    juce::var value (report);
    if (argc == 3 && juce::String (argv[1]) == "--report")
    {
        if (! juce::File (juce::String (argv[2])).replaceWithText (juce::JSON::toString (value)))
        { std::cerr << "Cannot write report\n"; return 2; }
    }
    std::cout << (passed ? "PASS: " : "FAIL: ") << records.size() << " M2 checks; no audio device opened\n";
    return passed ? 0 : 1;
}
