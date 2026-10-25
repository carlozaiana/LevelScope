#pragma once
#include "StandaloneSessionModel.h"

namespace levelscope::standalone::sessionfile
{
inline constexpr uint64_t maximumBytes = 128 * 1024 * 1024;
inline constexpr int schemaVersion = 2, programmeAssemblyVersion = 2;
struct ProfileReference { juce::String id, version, digest; };
// Value snapshot: each save owns one coherent revision; shared analyses are immutable.
struct Document
{
    juce::String id, programmeId, runtimeId;
    uint64_t revision = 0;
    int loadedSchema = schemaVersion;
    SourceDocument source;
    CurrentStateDocument current;
    ProfileReference profile;
    proposal::Applicability applicability = proposal::Applicability::unconfirmed;
    proposal::SectionMap sections;
    video::Reference reference;
    MeasurementDomain audition = MeasurementDomain::source;
    std::shared_ptr<const AnalysisSnapshot> sourceCache, currentCache;
    bool audioResolved = false, videoResolved = false;
    int64_t verifiedVideoBytes = 0, verifiedVideoModified = 0;
    juce::StringArray notices, verifiedAudioMembers;
};
struct Prepared
{
    StandaloneSessionModel model;
    juce::String notice;
    bool consumed = false;
};
struct Outcome
{
    std::optional<Document> document;
    ErrorCode error = ErrorCode::none;
    juce::String message;
    std::shared_ptr<Prepared> prepared = {}; // built/evaluated on the opening worker
    bool succeeded() const { return error == ErrorCode::none; }
};
// Failure injection belongs to tests; production uses an empty callback.
struct FileHooks
{
    std::function<bool(const char*)> fail;
    std::function<void(const char*)> reached;
    std::function<void(const juce::File&)> beforeCreate;
};
Document capture (StandaloneSessionModel&);
Outcome save (const Document&, const juce::File&, JobControl&, const FileHooks& = {});
Outcome read (const juce::File&, JobControl&); // core and independent caches; no media I/O
Outcome verify (Document, JobControl&);       // per-file verification, no devices
bool install (const Document&, StandaloneSessionModel&, juce::String&); // caller stops/closes audition first
bool adopt (const Document&, Prepared&, StandaloneSessionModel&, juce::String&);
bool saved (StandaloneSessionModel&, const Document&, const juce::File&);
// Relink changes only a locator; verification must still match saved bytes/format.
bool relinkAudio (Document&, const juce::File&);
bool relinkAudioMember (Document&, const juce::String& slotId, const juce::File&);
bool relinkVideo (Document&, const juce::File&);
}
