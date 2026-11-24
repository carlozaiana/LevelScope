// Same worker driver links against the qualified baseline and new production objects.
#include "CompositionTestFixtures.h"
#include "Standalone/StandaloneRegionProposal.h"
#include "Standalone/ProposalOperation.h"
#include <juce_cryptography/juce_cryptography.h>

namespace
{
struct PreparationCounts : JobFeedback
{
    JobPhase last = JobPhase::starting;
    int copies = 0, validations = 0, hashPasses = 0;
    void publish(JobPhase phase, double fraction) override
    {
        if (phase != last)
        {
            copies += phase == JobPhase::copy;
            validations += phase == JobPhase::samples;
            hashPasses += phase == JobPhase::sourceHash;
        }
        last = phase;
        JobFeedback::publish(phase, fraction);
    }
};
void optional(juce::MemoryOutputStream& bytes, std::optional<double> x)
{
    bytes.writeBool(bool(x));
    if (x) bytes.writeDouble(*x);
}
juce::String numericPayload(const measurement::Result& r, juce::String* encoded = nullptr)
{
    juce::MemoryOutputStream b;
    b.writeInt(r.sampleRate); b.writeInt64(r.frames);
    for (auto metric : {r.integrated, r.lra, r.momentaryMax, r.shortTermMax, r.samplePeak, r.truePeak})
    { b.writeInt(int(metric.status)); optional(b, metric.value); }
    for (auto value : {r.integratedRelativeGate, r.lraRelativeGate, r.p10, r.p95}) optional(b, value);
    b.writeInt64(r.integratedIncluded); b.writeInt64(r.lraIncluded); b.writeBool(r.lraStable);
    for (const auto* windows : {&r.momentary, &r.shortTerm, &r.lraWindows})
    {
        b.writeInt64(windows->size());
        for (const auto& w : *windows) { b.writeInt64(w.endFrame); b.writeDouble(w.energy); }
    }
    for (const auto* values : {&r.channelSamplePeaks, &r.channelTruePeaks})
    { b.writeInt64(values->size()); for (auto x : *values) b.writeDouble(x); }
    for (const auto* values : {&r.samplePeakFrames, &r.truePeakEighthFrames})
    { b.writeInt64(values->size()); for (auto x : *values) b.writeInt64(x); }
    if (encoded) *encoded = juce::Base64::toBase64(b.getData(), b.getDataSize());
    return juce::SHA256(b.getData(), b.getDataSize()).toHexString();
}
void compareProgramme(const juce::File& file, juce::DynamicObject& record)
{
    JobControl job;
    const auto media = validateMedia(file, job);
    expect(media.succeeded() && media.media.sha256 == "7ba6404041c94c909df290a6829a0fe92ce594329a1d575490c2bb4a5b6cc8f0", "agreed stereo development input");
    StandaloneSessionModel s;
    expect(s.setValidatedSource(file, media.media), "adopt supplied input");
    s.sourceDocument.layout = layoutFor(2);
    s.setSelectedTargetProfileIndex(1);
    s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
    auto source = measuredComposition(s, MeasurementDomain::source);
    s.synchronizeSections();
    record.setProperty("media_sha256", media.media.sha256);
    record.setProperty("source_numeric_sha256", numericPayload(source.detail));
    auto req = s.jobRequirements;
    req.integratedTenths = -230; req.integratedToleranceTenths = 2;
    req.truePeakTenths = -10; req.maximumLraTenths = int(std::floor((*source.detail.lra.value - 1.) * 10.));
    req.maximumMomentaryTenths.reset(); req.maximumShortTermTenths.reset();
    juce::String error;
    expect(s.setJobRequirements(req, error), "set comparison requirements");
    auto settings = proposal::RegionPlanningSettings::innerDefaults(proposal::RegionMode::automatic, sectionInputIdentity(s.measurementContext(MeasurementDomain::source)));
    settings.allowAutomaticMargins = true;
    settings.leftMarginMs = settings.rightMarginMs = 1000;
    expect(s.setRegionPlanning(settings, error), "select qualified competing-boundary method");
    juce::Array<juce::var> runs;
    SectionProposalOutcome last;
    for (int n = 0; n < 2; ++n)
    {
        JobControl worker;
        auto counter = std::make_shared<PreparationCounts>();
        worker.feedback = counter;
        const auto start = juce::Time::getMillisecondCounterHiRes();
        auto result = generateSectionProposal(s.sectionProposalRequest(), worker);
        const double total = juce::Time::getMillisecondCounterHiRes() - start;
        expect(result.value && (*result.value)->regionReview(), "proposal completes on unchanged baseline");
        const auto& review = *(*result.value)->regionReview();
        expect(!review.attempts.empty(), "comparison includes a measured candidate");
        juce::var row(new juce::DynamicObject());
        auto* run = row.getDynamicObject();
        run->setProperty("elapsed_ms", total); run->setProperty("preparation_ms", review.preparationSeconds * 1000);
        run->setProperty("copy_passes", counter->copies); run->setProperty("validation_passes", counter->validations);
        run->setProperty("validation_hash_passes", counter->hashPasses);
        run->setProperty("retained_attempt", int(review.retainedAttempt)); run->setProperty("outcome", (*result.value)->outcomeTitle());
        juce::Array<juce::var> alternatives;
        for (const auto& a : review.attempts)
        {
            juce::var v(new juce::DynamicObject()); auto* o = v.getDynamicObject();
            o->setProperty("plan", a.plan.envelope.audioKey()); o->setProperty("global_db", juce::String(a.plan.globalDb, 17));
            o->setProperty("verdict", int(a.evaluation.verdict));
            if (a.receipt) o->setProperty("candidate_numeric_sha256", numericPayload(a.receipt->measured()));
            alternatives.add(v);
        }
        run->setProperty("attempts", alternatives); runs.add(row);
        last = std::move(result);
    }
    record.setProperty("runs", runs);
    ProposalOperation op;
    expect(op.begin(s) && op.candidateComplete(s, last) && op.apply(s), "apply selected measured alternative through the accepted workflow");
    const auto measured = measuredComposition(s, MeasurementDomain::currentState);
    op.measurementComplete(s, true, {});
    expect(s.evaluation(MeasurementDomain::currentState).verdict == proposal::Verdict::pass, "independently measured applied Current passes");
    record.setProperty("applied_plan", proposalPlanIdentity(s.currentStateDocument.processingPlan()));
    record.setProperty("applied_choice", int(op.selectedChoice));
    record.setProperty("processed_current_numeric_sha256", numericPayload(measured.detail));
    auto output = ex::render(ex::capture(s), juce::File::getSpecialLocation(juce::File::tempDirectory), job);
    expect(bool(output.value), "render PCM and independently measure encoded bytes");
    record.setProperty("rendered_pcm_sha256", (*output.value)->hash());
    juce::String encoded;
    record.setProperty("rendered_numeric_sha256", numericPayload((*output.value)->analysis().detail, &encoded));
    record.setProperty("rendered_numeric_base64", encoded);
}
}
int sessionSourceBenchmark(int argc, char** argv)
{
    if (argc != 3 || juce::File(argv[2]).exists()) return 2;
    juce::var result(new juce::DynamicObject()); auto* record = result.getDynamicObject();
    juce::String failure;
    try { compareProgramme(juce::File(argv[1]), *record); }
    catch (const std::exception& e) { failure = e.what(); }
    record->setProperty("build", LEVELSCOPE_MEASUREMENT_BUILD); record->setProperty("failure", failure);
    record->setProperty("passed", failure.isEmpty()); record->setProperty("checks", checks);
    std::cout << failure << std::endl;
    return juce::File(argv[2]).replaceWithText(juce::JSON::toString(result)) && failure.isEmpty() ? 0 : 1;
}

#ifndef LEVELSCOPE_EMBEDDED_SOURCE_BENCHMARK
int main(int argc, char** argv) { return sessionSourceBenchmark(argc, argv); }
#endif
