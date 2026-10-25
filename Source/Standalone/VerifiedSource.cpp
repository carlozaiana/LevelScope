#include "VerifiedSource.h"
#include <juce_cryptography/juce_cryptography.h>

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
class HashInput final : public juce::InputStream
{
public:
    HashInput (juce::FileInputStream& s, JobControl& c) : stream (s, 65536), control (c), length (s.getTotalLength()) {}
    juce::int64 getTotalLength() override { return length; }
    juce::int64 getPosition() override { return stream.getPosition(); }
    bool setPosition (juce::int64 p) override { return stream.setPosition (p); }
    bool isExhausted() override { return control.isCancelled() || stream.isExhausted(); }
    int read (void* dst, int n) override
    {
        checkCancel (control);
        const auto got = stream.read (dst, n);
        control.setProgress (0.15f + 0.15f * static_cast<float> (static_cast<double> (stream.getPosition()) / std::max<juce::int64> (1, length)));
        return got;
    }
private:
    juce::BufferedInputStream stream;
    JobControl& control;
    const juce::int64 length;
};
}
VerifiedSource::~VerifiedSource() { if (ownedFile != juce::File()) { ownedFile.setReadOnly (false); ownedFile.deleteFile(); } }
std::vector<juce::File> VerifiedSource::files() const
{
    if (members.empty()) return { ownedFile };
    std::vector<juce::File> result; for (const auto& m : members) result.push_back(m->file()); return result;
}
std::unique_ptr<ProgrammeReader> VerifiedSource::createReader() const { return ProgrammeReader::open(files(), media); }
bool VerifiedSource::matches (const MediaInfo& other) const
{
    return media.sha256 == other.sha256 && media.bytes == other.bytes && media.frames == other.frames
        && media.channels == other.channels && media.sampleRate == other.sampleRate
        && media.bitsPerSample == other.bitsPerSample && media.format == other.format
        && (members.empty() ? ownedFile.existsAsFile() : std::all_of(members.begin(), members.end(), [](const auto& m) { return m->file().existsAsFile(); }));
}
VerifiedSourceResult prepareVerifiedSource (const ServiceContext& context, JobControl& control, const SnapshotStorage& storage)
{
    try
    {
        checkCancel (control);
        if (context.media.bytes <= 0 || context.media.sha256.isEmpty())
            throw Failure { ErrorCode::invalidRequest, "Imported source identity required" };
        if (! storage.directory.isDirectory())
            throw Failure { ErrorCode::writeFailed, "Temporary snapshot directory is unavailable" };
        const auto free = storage.availableBytes ? storage.availableBytes (storage.directory) : storage.directory.getBytesFreeOnVolume();
        if (free >= 0 && (free < context.media.bytes || free - context.media.bytes < 65536))
            throw Failure { ErrorCode::insufficientDiskSpace, "Insufficient temporary disk space for a verified source copy" };
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
        auto input = context.sourceFile.createInputStream();
        if (! input) throw Failure { ErrorCode::missingFile, "Source is unavailable" };
        if (input->getTotalLength() != context.media.bytes)
            throw Failure { ErrorCode::sourceChanged, "Source length changed; import again" };
        auto owned = std::shared_ptr<VerifiedSource> (new VerifiedSource (
            storage.directory.getChildFile ("levelscope-" + juce::Uuid().toString() + ".wav"), context.media));
        auto output = storage.createWriter ? storage.createWriter (owned->file()) : std::make_unique<FileWriter> (owned->file());
        if (! output) throw Failure { ErrorCode::writeFailed, "Cannot create private source copy" };
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
            control.setProgress (0.15f * static_cast<float> (static_cast<double> (copied) / context.media.bytes));
        }
        if (! output->flush()) throw Failure { ErrorCode::writeFailed, "Temporary source copy flush failed" };
        output.reset();
        if (! input->isExhausted() || input->getTotalLength() != context.media.bytes)
            throw Failure { ErrorCode::sourceChanged, "Source length changed during preparation" };
        auto hashStream = owned->file().createInputStream();
        if (! hashStream) throw Failure { ErrorCode::missingFile, "Cannot reopen private source copy" };
        HashInput hashInput (*hashStream, control);
        const auto hash = juce::SHA256 (hashInput).toHexString();
        checkCancel (control);
        if (hashStream->getStatus().failed() || hashStream->getPosition() != context.media.bytes)
            throw Failure { ErrorCode::decodeFailed, "Incomplete source copy hash" };
        if (hash != context.media.sha256)
            throw Failure { ErrorCode::sourceChanged, "Source content changed; import again" };
        hashStream.reset();
        if (! owned->file().setReadOnly (true)) throw Failure { ErrorCode::writeFailed, "Cannot protect private source copy" };
        JobControl validation; validation.parent = &control;
        validation.parentProgress = &control.progress; validation.progressBase = 0.3f; validation.progressScale = 0.2f;
        const auto verified = validateMedia (owned->file(), validation);
        checkCancel (control);
        if (! verified.succeeded()) throw Failure { verified.error, verified.message };
        if (! owned->matches (verified.media))
            throw Failure { ErrorCode::sourceChanged, "Source copy metadata differs from imported source" };
        control.setProgress (0.5f);
        return { std::shared_ptr<const VerifiedSource> (std::move (owned)), ErrorCode::none, "Source copy verified" };
    }
    catch (const Failure& e) { return { {}, e.code, e.message }; }
    catch (const std::exception& e) { return { {}, ErrorCode::unavailable, e.what() }; }
}
}
