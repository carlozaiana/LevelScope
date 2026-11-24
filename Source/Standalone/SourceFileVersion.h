#pragma once
#include <juce_core/juce_core.h>
#include <memory>
#include <sys/stat.h>

namespace levelscope::standalone
{
// Runtime change detector. Only full verification establishes content identity.
// A captured version monitors subsequent changes even on coarse-timestamp volumes.
struct SourceFileVersion
{
    uint64_t device = 0, inode = 0;
    int64_t size = 0, seconds = 0, nanos = 0;
    bool valid = false;

    // read() is a one-time identity observation for pathname/cleanup checks.
    // capture() begins monitoring before validating an original or sealed copy.
    static SourceFileVersion read(const juce::File&, bool followLinks = true);
    static SourceFileVersion capture(const juce::File&, bool followLinks = true);
    bool matches(const juce::File&, bool followLinks = true) const;
    bool sameFile(const SourceFileVersion& other) const
    {
        return valid && other.valid && device == other.device && inode == other.inode;
    }
private:
    static SourceFileVersion fromStat(const struct stat&);
#if JUCE_MAC
    struct Watch;
    std::shared_ptr<Watch> watch;
#endif
};
}
