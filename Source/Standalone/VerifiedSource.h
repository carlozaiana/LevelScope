#pragma once
#include "StandaloneServices.h"
#include <functional>
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
};
class VerifiedSource
{
public:
    ~VerifiedSource();
    const juce::File& file() const { return ownedFile; }
    bool matches (const MediaInfo&) const;
    std::unique_ptr<ProgrammeReader> createReader() const;
    std::vector<juce::File> files() const;
private:
    VerifiedSource (juce::File f, MediaInfo m) : ownedFile (std::move (f)), media (std::move (m)) {}
    juce::File ownedFile;
    MediaInfo media;
    std::vector<std::shared_ptr<const VerifiedSource>> members;
    friend ServiceResult<std::shared_ptr<const VerifiedSource>> prepareVerifiedSource (const ServiceContext&, JobControl&, const SnapshotStorage&);
};
using VerifiedSourceResult = ServiceResult<std::shared_ptr<const VerifiedSource>>;
VerifiedSourceResult prepareVerifiedSource (const ServiceContext&, JobControl&, const SnapshotStorage& = {});
}
