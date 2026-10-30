#pragma once
#include "Pcm24Wave.h"
#include "StandaloneSessionModel.h"
#include "VerifiedSource.h"

namespace levelscope::standalone::exporting
{
inline constexpr const char* version = "verified-pcm-export-m10.1";
struct Request
{
    ServiceContext source, current;
    std::shared_ptr<const AnalysisSnapshot> sourceAnalysis, currentAnalysis;
    proposal::Applicability applicability = proposal::Applicability::unconfirmed;
    juce::String documentId, programmeId, mediaIssue;
    uint64_t documentRevision = 0;
    bool unsaved = true;
    juce::File sessionFile;
    proposal::SectionMap sections;
    proposal::JobRequirements jobRequirements;
    proposal::ProposalPolicy proposalPolicy;
    std::optional<proposal::ProposalPolicy> lastAppliedProposalPolicy;
    juce::String lastAppliedProposalBinding,proposalInspection,lastAppliedProposalInspection;
    std::optional<proposal::markers::Receipt> markerImport;
    video::Reference video;
};
Request capture (const StandaloneSessionModel&);
struct Specification
{
    pcm24::Format format;
    std::vector<int> outputToInput;
    ChannelLayout layout;
    uint64_t seed = 0;
    juce::String identity, seedHex;
};
juce::String specification (const Request&, Specification&);
juce::String readiness (const Request&);
juce::String preview (const Request&, bool includeReadiness = true);
juce::String audioIdentity (const Request&);
juce::String requestIdentity (const Request&);

class Candidate
{
public:
    ~Candidate();
    const AnalysisSnapshot& analysis() const { return *measured; }
    const Specification& spec() const { return output; }
    const juce::String& hash() const { return sha256; }
    juce::File file() const { return directory.getChildFile ("programme.wav"); }
    bool matchesAudio (const Request&) const;
    bool isPublished() const { return published; }
    juce::File folder() const { return directory; }
    std::shared_ptr<const VerifiedSource> verifiedSource() const { return sourceCopy; }
    proposal::MeasuredProgramme measuredInput() const { return artifactInput; }
private:
    Candidate() = default;
    juce::File directory;
    uint64_t directoryDevice = 0, directoryInode = 0;
    bool published = false;
    Specification output;
    juce::String originAudioIdentity, sha256, runtimeId;
    std::shared_ptr<const AnalysisSnapshot> measured;
    std::shared_ptr<const VerifiedSource> sourceCopy;
    proposal::MeasuredProgramme artifactInput;
    friend ServiceResult<std::shared_ptr<Candidate>> render (const Request&, const juce::File&, JobControl&,
                                                            const pcm24::Hooks&, ProcessedAudioObserver*, std::shared_ptr<const VerifiedSource>);
    friend ServiceResult<juce::File> publish (Candidate&, const Request&, const juce::File&, const juce::String&,
                                             JobControl&, const pcm24::Hooks&);
};
proposal::Evaluation artifactEvaluation (const Candidate&, const Request&);
bool rulesPermitPublication(const proposal::Evaluation& current,const proposal::Evaluation& artifact);
juce::String approvalIdentity (const Candidate&, const Request&);
// Both services are worker-only. The caller locks edits and checks runtime/job identity on adoption.
ServiceResult<std::shared_ptr<Candidate>> render (const Request&, const juce::File& destinationParent,
                                                JobControl&, const pcm24::Hooks& = {}, ProcessedAudioObserver* = nullptr,
                                                std::shared_ptr<const VerifiedSource> reusableSource = {});
ServiceResult<juce::File> publish (Candidate&, const Request& fresh, const juce::File& destination,
                                 const juce::String& explicitApproval, JobControl&, const pcm24::Hooks& = {});
ServiceResult<juce::File> measurementReport (const Request&, const juce::File& destination,
                                           JobControl&, const pcm24::Hooks& = {});
ServiceResult<juce::File> artifactReport (const Candidate&, const Request&, const juce::File& destination,
                                        JobControl&, const pcm24::Hooks& = {});
// Derived model, never read back as measurement authority. Null candidate means no artifact fields.
juce::var report (const Request&, const Candidate* = nullptr, bool publishing = false);
juce::String reportHtml (const juce::var&);
juce::String reportJson (const juce::var&);
}
