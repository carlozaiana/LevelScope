#pragma once
#include <juce_core/juce_core.h>
#include <sys/stat.h>

namespace levelscope::standalone
{
// Runtime change detector. Only full verification establishes content identity.
// The change time also detects writes followed by restoration of the old mtime.
struct SourceFileVersion
{
    uint64_t device = 0, inode = 0;
    int64_t size = 0, seconds = 0, nanos = 0;
    bool valid = false;

    static SourceFileVersion read(const juce::File& file, bool followLinks = true)
    {
        struct stat s {};
        const char* path = file.getFullPathName().toRawUTF8();
        if ((followLinks ? ::stat(path, &s) : ::lstat(path, &s)) != 0 || !S_ISREG(s.st_mode)) return {};
#if JUCE_MAC
        const auto changed = s.st_ctimespec;
#else
        const auto changed = s.st_ctim;
#endif
        return {uint64_t(s.st_dev), uint64_t(s.st_ino), s.st_size, changed.tv_sec, changed.tv_nsec, true};
    }
    bool sameFile(const SourceFileVersion& other) const
    {
        return valid && other.valid && device == other.device && inode == other.inode;
    }
    bool operator==(const SourceFileVersion& other) const
    {
        return sameFile(other) && size == other.size && seconds == other.seconds && nanos == other.nanos;
    }
};
}
