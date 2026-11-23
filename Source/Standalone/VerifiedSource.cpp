#include "VerifiedSource.h"
#include <fcntl.h>
#include <unistd.h>

namespace levelscope::standalone
{
namespace
{
struct Failure { ErrorCode code; juce::String message; };
void checkCancel (JobControl& c)
{
    if (c.isCancelled()) throw Failure { ErrorCode::cancelled, "Source preparation cancelled" };
}
struct FileWriter final : SnapshotWriter
{
    explicit FileWriter (const juce::File& file) : stream (file.createOutputStream()) {}
    bool write (const void* data, size_t n) override { return stream && stream->write (data, n); }
    bool flush() override { if (! stream) return false; stream->flush(); return stream->getStatus().wasOk(); }
    std::unique_ptr<juce::FileOutputStream> stream;
};

}
VerifiedSource::~VerifiedSource()
{
    // Never chmod or unlink a replacement at a formerly owned pathname.
    if (sealedVersion.sameFile(SourceFileVersion::read(ownedFile, false))) ownedFile.deleteFile();
}
bool VerifiedSource::intact() const
{
    if (!members.empty())
        return std::all_of(members.begin(), members.end(), [](const auto& m) { return m->intact(); });
    return sealedVersion == SourceFileVersion::read(ownedFile, false);
}
bool ProgrammeReader::sourceIntact() const { return !sourceLease || sourceLease->intact(); }
std::vector<juce::File> VerifiedSource::files() const
{
    if (members.empty()) return { ownedFile };
    std::vector<juce::File> result; for (const auto& m : members) {const auto f=m->files();result.insert(result.end(),f.begin(),f.end());} return result;
}
std::unique_ptr<ProgrammeReader> VerifiedSource::createReader() const
{
    if (!intact()) return {};
    auto reader = composition
        ? ProgrammeReader::compose(members[0]->createReader(), members[1]->createReader(),
                                   composition->map(false), composition->map(true))
        : ProgrammeReader::open(files(), media);
    if (!reader || !intact()) return {};
    reader->sourceLease = shared_from_this();
    return reader;
}
bool VerifiedSource::matches (const MediaInfo& other) const
{
    return media.sha256 == other.sha256 && media.bytes == other.bytes && media.frames == other.frames
        && media.channels == other.channels && media.sampleRate == other.sampleRate
        && media.bitsPerSample == other.bitsPerSample && media.format == other.format
        && intact();
}
namespace
{
juce::String assetKey(const juce::File& file, const MediaInfo& media)
{
    return media.sha256 + ":" + file.getFullPathName();
}
}
SessionVerifiedAssets::Counts SessionVerifiedAssets::counts() const
{
    std::lock_guard<std::mutex> lock(mutex);
    return totals;
}
uint64_t SessionVerifiedAssets::missingBytes(const std::vector<ProgrammeMember>& members) const
{
    std::lock_guard<std::mutex> lock(mutex);
    uint64_t missing = 0;
    for (const auto& member : members)
    {
        const auto found = assets.find(assetKey(member.file, member.media));
        if (found == assets.end() || found->second->preparing || !found->second->result.value
            || !(*found->second->result.value)->matches(member.media))
            missing += uint64_t(member.media.bytes);
    }
    return missing;
}
std::shared_ptr<SessionVerifiedAssets> SessionVerifiedAssets::retaining(const std::vector<ProgrammeMember>& members) const
{
    auto next = std::make_shared<SessionVerifiedAssets>();
    std::lock_guard<std::mutex> lock(mutex);
    for (const auto& member : members)
    {
        const auto key = assetKey(member.file, member.media);
        const auto found = assets.find(key);
        if (found == assets.end() || found->second->preparing || !found->second->result.value) continue;
        if (!member.media.matchesValidatedFile(member.file) || !(*found->second->result.value)->matches(member.media)) continue;
        auto retained = std::make_shared<Asset>();
        retained->preparing = false;
        retained->result = found->second->result;
        next->assets.emplace(key, std::move(retained));
    }
    return next;
}
VerifiedSourceResult SessionVerifiedAssets::acquire(const ServiceContext& context, JobControl& control)
{
    const auto key = assetKey(context.sourceFile, context.media);
    std::shared_ptr<Asset> asset;
    {
        std::unique_lock<std::mutex> lock(mutex);
        for (;;)
        {
            checkCancel(control);
            if (!context.media.matchesValidatedFile(context.sourceFile))
                return {{}, ErrorCode::sourceChanged, "Original source missing or changed; revalidate or relink"};
            auto found = assets.find(key);
            if (found == assets.end())
            {
                asset = std::make_shared<Asset>();
                assets.emplace(key, asset);
                break;
            }
            asset = found->second;
            if (asset->preparing)
            {
                ++totals.waits;
                while (asset->preparing && !control.isCancelled())
                    asset->changed.wait_for(lock, std::chrono::milliseconds(20));
                checkCancel(control);
                // A cancelled producer discarded its partial file. Another live
                // request can now take over, with at most one copier per asset.
                if (!asset->result.value && asset->result.error != ErrorCode::cancelled) return asset->result;
                continue;
            }
            if (asset->result.value && (*asset->result.value)->matches(context.media))
            {
                ++totals.cacheHits;
                control.setProgress(0.5f);
                return asset->result;
            }
            assets.erase(found); // private copy disappeared or changed; verify a fresh copy
        }
    }
    auto result = VerifiedSource::copyAsset(context, control, {});
    {
        std::lock_guard<std::mutex> lock(mutex);
        asset->result = result;
        asset->preparing = false;
        if (result.value)
        {
            ++totals.copies;
            totals.copiedBytes += uint64_t(context.media.bytes);
        }
        else assets.erase(key);
    }
    asset->changed.notify_all();
    return result;
}
VerifiedSourceResult prepareVerifiedSource (const ServiceContext& context, JobControl& control, const SnapshotStorage& storage)
{
    try
    {
        checkCancel (control);
        if (context.media.bytes <= 0 || context.media.sha256.isEmpty())
            throw Failure { ErrorCode::invalidRequest, "Imported source identity required" };
        if ((context.composition || context.assembly) && (!context.verifiedAssets || storage.availableBytes
            || storage.createWriter || storage.directory != juce::File::getSpecialLocation(juce::File::tempDirectory)))
        {
            const auto free = storage.availableBytes ? storage.availableBytes(storage.directory) : storage.directory.getBytesFreeOnVolume();
            if (free >= 0 && (free < context.media.bytes || free - context.media.bytes < 65536))
                throw Failure { ErrorCode::insufficientDiskSpace, "Insufficient temporary space for the programme" };
        }
        if(context.composition)
        {
            if(!context.isValid())throw Failure{ErrorCode::invalidRequest,"Exact confirmed composition required"};
            if(!context.sourceFilesMatch())throw Failure{ErrorCode::sourceChanged,"A required composition member is missing, changed or duplicated"};
            auto group=std::shared_ptr<VerifiedSource>(new VerifiedSource({},context.media));group->composition=context.composition;
            juce::int64 completeBytes=0;
            for(const auto* component:{&context.composition->me,&context.composition->dialogue})
            {
                ServiceContext part=context;part.composition.reset();part.sourceFile=component->file;part.media=component->media;
                part.assembly=component->assembly;part.fileLayout=component->layout;part.range={0,part.media.frames};
                JobControl child;child.parent=&control;child.parentProgress=&control.progress;
                child.progressBase=0.5f*static_cast<float>(static_cast<double>(completeBytes)/context.media.bytes);
                child.progressScale=static_cast<float>(static_cast<double>(component->media.bytes)/context.media.bytes);
                auto copy=prepareVerifiedSource(part,child,storage);
                if(!copy.value)throw Failure{copy.error,component->purpose+": "+copy.message};
                group->members.push_back(*copy.value);completeBytes+=component->media.bytes;
            }
            checkCancel(control);
            if(!context.sourceFilesMatch())throw Failure{ErrorCode::sourceChanged,"A composition member changed during preparation"};
            control.setProgress(0.5f);
            return {std::shared_ptr<const VerifiedSource>(std::move(group)),ErrorCode::none,"Both complete composition components verified"};
        }
        if (context.assembly)
        {
            if (!context.isValid()) throw Failure { ErrorCode::invalidRequest, "Confirmed exact programme assembly required" };
            if (!context.sourceFilesMatch()) throw Failure { ErrorCode::sourceChanged, "A programme member changed, disappeared or aliases another member" };
            auto group = std::shared_ptr<VerifiedSource>(new VerifiedSource({}, context.media));
            juce::int64 completeBytes = 0;
            for (const auto& member : context.assembly->members)
            {
                checkCancel(control);
                ServiceContext part = context; part.assembly.reset(); part.sourceFile = member.file; part.media = member.media;
                JobControl child; child.parent = &control; child.parentProgress = &control.progress;
                child.progressBase = 0.5f * static_cast<float>(static_cast<double>(completeBytes) / context.media.bytes);
                child.progressScale = static_cast<float>(static_cast<double>(member.media.bytes) / context.media.bytes);
                const auto copy = prepareVerifiedSource(part, child, storage);
                if (!copy.value) throw Failure { copy.error, member.file.getFileName() + ": " + copy.message };
                group->members.push_back(*copy.value); completeBytes += member.media.bytes;
                control.setProgress(0.5f * static_cast<float>(static_cast<double>(completeBytes) / context.media.bytes));
            }
            checkCancel(control);
            if (!context.sourceFilesMatch()) throw Failure { ErrorCode::sourceChanged, "A programme member changed during preparation" };
            return { std::shared_ptr<const VerifiedSource>(std::move(group)), ErrorCode::none, "Complete programme group verified" };
        }
        if (context.verifiedAssets && !storage.availableBytes && !storage.createWriter
            && storage.directory == juce::File::getSpecialLocation(juce::File::tempDirectory))
        {
            auto owner = context.verifiedAssets->lock();
            if (!owner) return {{}, ErrorCode::staleRevision, "Source session closed"};
            return owner->acquire(context, control);
        }
        return VerifiedSource::copyAsset(context, control, storage);
    }
    catch (const Failure& e) { return { {}, e.code, e.message }; }
    catch (const std::exception& e) { return { {}, ErrorCode::unavailable, e.what() }; }
}
VerifiedSourceResult VerifiedSource::copyAsset(const ServiceContext& context, JobControl& control, const SnapshotStorage& storage)
{
    try
    {
        checkCancel(control);
        if (!context.media.matchesValidatedFile(context.sourceFile))
            throw Failure { ErrorCode::sourceChanged, "Original source missing or changed; revalidate or relink" };
        if (! storage.directory.isDirectory())
            throw Failure { ErrorCode::writeFailed, "Temporary snapshot directory is unavailable" };
        const auto free = storage.availableBytes ? storage.availableBytes (storage.directory) : storage.directory.getBytesFreeOnVolume();
        if (free >= 0 && (free < context.media.bytes || free - context.media.bytes < 65536))
            throw Failure { ErrorCode::insufficientDiskSpace, "Insufficient temporary disk space for a verified source copy" };
        auto input = context.sourceFile.createInputStream();
        if (! input) throw Failure { ErrorCode::missingFile, "Source is unavailable" };
        if (input->getTotalLength() != context.media.bytes)
            throw Failure { ErrorCode::sourceChanged, "Source length changed; import again" };
        auto owned = std::shared_ptr<VerifiedSource> (new VerifiedSource (
            storage.directory.getChildFile ("levelscope-" + juce::Uuid().toString() + ".wav"), context.media));
        const auto path = owned->file().getFullPathName();
        const int fd = ::open(path.toRawUTF8(), O_CREAT | O_EXCL | O_WRONLY, 0600);
        if (fd < 0) throw Failure { ErrorCode::writeFailed, "Cannot create private source copy" };
        ::close(fd);
        owned->sealedVersion = SourceFileVersion::read(owned->file(), false);
        auto output = storage.createWriter ? storage.createWriter (owned->file()) : std::make_unique<FileWriter> (owned->file());
        if (! output) throw Failure { ErrorCode::writeFailed, "Cannot create private source copy" };
        control.reportPhase(JobPhase::copy, 0);
        std::array<char, 65536> buffer;
        juce::int64 copied = 0;
        while (copied < context.media.bytes)
        {
            checkCancel (control);
            const auto wanted = static_cast<int> (std::min<juce::int64> (buffer.size(), context.media.bytes - copied));
            if (input->read (buffer.data(), wanted) != wanted)
                throw Failure { ErrorCode::decodeFailed, "Snapshot read failed" };
            if (! output->write (buffer.data(), static_cast<size_t> (wanted)))
                throw Failure { ErrorCode::writeFailed, "Temporary source copy write failed; check available disk space" };
            copied += wanted;
            control.reportPhase(JobPhase::copy, double(copied) / context.media.bytes);
            control.setProgress (0.15f * static_cast<float> (static_cast<double> (copied) / context.media.bytes));
        }
        if (! output->flush()) throw Failure { ErrorCode::writeFailed, "Temporary source copy flush failed" };
        output.reset();
        if (! input->isExhausted() || input->getTotalLength() != context.media.bytes)
            throw Failure { ErrorCode::sourceChanged, "Source length changed during preparation" };
        if (!owned->sealedVersion.sameFile(SourceFileVersion::read(owned->file(), false)))
            throw Failure { ErrorCode::sourceChanged, "Private source copy replaced during preparation" };
        if (! owned->file().setReadOnly (true)) throw Failure { ErrorCode::writeFailed, "Cannot protect private source copy" };
        owned->sealedVersion = SourceFileVersion::read(owned->file(), false);
        JobControl validation; validation.parent = &control;
        validation.parentProgress = &control.progress; validation.progressBase = 0.15f; validation.progressScale = 0.35f;
        const auto verified = validateMedia (owned->file(), validation);
        checkCancel (control);
        if (! verified.succeeded()) throw Failure { verified.error, verified.message };
        if (!context.media.matchesValidatedFile(context.sourceFile) || ! owned->matches (verified.media))
            throw Failure { ErrorCode::sourceChanged, "Verified copy differs from the imported source; revalidate or relink" };
        control.setProgress (0.5f);
        return { std::shared_ptr<const VerifiedSource> (std::move (owned)), ErrorCode::none, "Source copy verified" };
    }
    catch (const Failure& e) { return { {}, e.code, e.message }; }
    catch (const std::exception& e) { return { {}, ErrorCode::unavailable, e.what() }; }
}
}
