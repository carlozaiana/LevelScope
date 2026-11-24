#include "SourceFileVersion.h"

#if JUCE_MAC
#include <cerrno>
#include <fcntl.h>
#include <mutex>
#include <sys/event.h>
#include <unistd.h>
#endif

namespace levelscope::standalone
{
SourceFileVersion SourceFileVersion::fromStat(const struct stat& info)
{
    SourceFileVersion version;
    if (!S_ISREG(info.st_mode)) return version;
#if JUCE_MAC
    const auto change = info.st_ctimespec;
#else
    const auto change = info.st_ctim;
#endif
    version.device = info.st_dev;
    version.inode = info.st_ino;
    version.size = info.st_size;
    version.seconds = change.tv_sec;
    version.nanos = change.tv_nsec;
    version.valid = true;
    return version;
}

SourceFileVersion SourceFileVersion::read(const juce::File& file, bool followLinks)
{
    const auto* path = file.getFullPathName().toRawUTF8();
    struct stat info {};
    return (followLinks ? ::stat(path, &info) : ::lstat(path, &info)) == 0
        ? fromStat(info) : SourceFileVersion{};
}

#if JUCE_MAC
// Copies of a validated MediaInfo share this watch. The kernel records changes
// between checks, including writes that leave HFS+ size and timestamps equal.
struct SourceFileVersion::Watch
{
    int file = -1, queue = -1;
    std::mutex mutex;
    bool invalidated = false;

    ~Watch()
    {
        if (queue >= 0) ::close(queue);
        if (file >= 0) ::close(file);
    }

    bool arm(const juce::File& path, bool followLinks)
    {
        file = ::open(path.getFullPathName().toRawUTF8(), O_EVTONLY | O_CLOEXEC | (followLinks ? 0 : O_NOFOLLOW));
        if (file < 0) return false;
        queue = ::kqueue();
        if (queue < 0 || ::fcntl(queue, F_SETFD, FD_CLOEXEC) < 0) return false;
        struct kevent change {};
        EV_SET(&change, static_cast<uintptr_t>(file), EVFILT_VNODE, EV_ADD | EV_ONESHOT,
               NOTE_WRITE | NOTE_EXTEND | NOTE_ATTRIB | NOTE_LINK | NOTE_RENAME | NOTE_DELETE | NOTE_REVOKE, 0, nullptr);
        return ::kevent(queue, &change, 1, nullptr, 0, nullptr) == 0;
    }

    bool unchanged()
    {
        std::lock_guard<std::mutex> lock(mutex);
        if (!invalidated)
        {
            struct kevent event {};
            const struct timespec immediate {};
            int result;
            do { result = ::kevent(queue, nullptr, 0, &event, 1, &immediate); }
            while (result < 0 && errno == EINTR);
            invalidated = result != 0;
        }
        return !invalidated;
    }
};
#endif

SourceFileVersion SourceFileVersion::capture(const juce::File& file, bool followLinks)
{
#if JUCE_MAC
    auto monitor = std::make_shared<Watch>();
    if (!monitor->arm(file, followLinks)) return {};
    struct stat info {};
    if (::fstat(monitor->file, &info) != 0) return {};
    auto version = fromStat(info);
    version.watch = std::move(monitor);
    return version;
#else
    return read(file, followLinks);
#endif
}

bool SourceFileVersion::matches(const juce::File& file, bool followLinks) const
{
    const auto observed = read(file, followLinks);
    if (!sameFile(observed) || size != observed.size || seconds != observed.seconds || nanos != observed.nanos) return false;
#if JUCE_MAC
    return watch && watch->unchanged();
#else
    return true;
#endif
}
}
