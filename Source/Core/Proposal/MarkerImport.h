#pragma once
#include "Model/ProgrammeSections.h"
#include <juce_cryptography/juce_cryptography.h>

namespace levelscope::proposal::markers
{
inline constexpr const char* version = "point-marker-import-m15.1";
inline constexpr size_t maximumBytes = 8 * 1024 * 1024, maximumRows = 10000;
enum class Format { nuendoCsv, proToolsText, pastedPositions };
enum class Unit { audioFrames, seconds, videoFrames, nonDropTimecode };
struct Options
{
    Format format = Format::nuendoCsv;
    Unit unit = Unit::videoFrames;
    char separator = ',';
    bool header = true, confirmedPointMarkers = false, confirmedUnits = false;
    int positionColumn = 0, nameColumn = 1, typeColumn = -1;
    int rateNumerator = 25, rateDenominator = 1;
    juce::String origin = "0", initialName = "Programme start";
};
struct Row
{
    int line = 0;
    juce::String position, name, seconds, quantization;
    int64_t frame = 0;
    bool coalesced = false;
};
struct Receipt
{
    Options options;
    juce::String basename, inputHash;
    std::vector<Row> rows;
    juce::StringArray notices;
    SectionMap imported;
    juce::String identity() const;
};
struct Preview
{
    Receipt receipt;
    juce::String previousMap, runtime;
    juce::StringArray diff;
    bool mapChanged = false;
};
using Cancel = std::function<bool()>;
// Throws a row-specific invalid_argument on malformed, unsupported or ambiguous input.
// Input is a bounded immutable byte snapshot; this module never opens media/devices.
Preview prepare (const std::string& utf8, const juce::String& basename, const Options&,
                 const SectionMap&, const juce::String& runtime, const Cancel& = {});
// Used by persistence and reports. All integers are exact decimal strings.
juce::var encode (const Receipt&);
Receipt decode (const juce::var&);
juce::String contentHash (const std::string&);
juce::String mapDescription (const SectionMap&);
}
