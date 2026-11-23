#pragma once
#include "StandaloneServices.h"
#include <functional>
#include <condition_variable>
#include <map>
#include <mutex>
#include <memory>

namespace levelscope::standalone
{
// Worker-only storage boundary, injectable for deterministic disk failure tests.
struct SnapshotWriter
{
    virtual ~SnapshotWriter() = default;
    virtual bool write (const void*, size_t) = 0;
    virtual bool flush() = 0;
};
struct SnapshotStorage
{
    juce::File directory = juce::File::getSpecialLocation (juce::File::tempDirectory);
    std::function<juce::int64 (const juce::File&)> availableBytes;
    std::function<std::unique_ptr<SnapshotWriter> (const juce::File&)> createWriter;

    bool usesSessionStorage() const
    {
        return !availableBytes && !createWriter && directory == juce::File::getSpecialLocation(juce::File::tempDirectory);
    }
};
class VerifiedSource : public std::enable_shared_from_this<VerifiedSource>
{
public:
    ~VerifiedSource();
    const juce::File& file() const { return ownedFile; }
    bool matches (const MediaInfo&) const;
    bool intact() const;
    std::unique_ptr<ProgrammeReader> createReader() const;
    std::vector<juce::File> files() const;
private:
    VerifiedSource (juce::File f, MediaInfo m) : ownedFile (std::move (f)), media (std::move (m)) {}
    juce::File ownedFile;
    MediaInfo media;
    SourceFileVersion sealedVersion;
    std::vector<std::shared_ptr<const VerifiedSource>> members;
    std::optional<ProgrammeComposition> composition;
    friend class SessionVerifiedAssets;
    static ServiceResult<std::shared_ptr<const VerifiedSource>> copyAsset(const ServiceContext&, JobControl&, const SnapshotStorage&);
    friend ServiceResult<std::shared_ptr<const VerifiedSource>> prepareVerifiedSource (const ServiceContext&, JobControl&, const SnapshotStorage&);
};
using VerifiedSourceResult = ServiceResult<std::shared_ptr<const VerifiedSource>>;
// The owner changes on New/Open or source replacement. Completed, still-valid
// assets may be carried across component/role edits within a live session.
class SessionVerifiedAssets
{
public:
    struct Counts { uint64_t copies = 0, copiedBytes = 0, cacheHits = 0, waits = 0; };
    Counts counts() const;
    uint64_t missingBytes(const std::vector<ProgrammeMember>&) const;
    std::shared_ptr<SessionVerifiedAssets> retaining(const std::vector<ProgrammeMember>&) const;
    VerifiedSourceResult acquire(const ServiceContext&, JobControl&);
private:
    struct Asset
    {
        bool preparing = true;
        std::condition_variable changed;
        VerifiedSourceResult result;
    };
    mutable std::mutex mutex;
    std::map<juce::String, std::shared_ptr<Asset>> assets;
    Counts totals;
};
VerifiedSourceResult prepareVerifiedSource (const ServiceContext&, JobControl&, const SnapshotStorage& = {});
}
