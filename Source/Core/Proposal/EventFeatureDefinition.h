#pragma once
#include <juce_core/juce_core.h>

namespace levelscope::proposal
{
inline constexpr const char* eventConstraintsVersion="region-constraints-v2";
inline constexpr const char* eventFormationVersion="region-formation-v2";
inline constexpr const char* eventSearchVersion="region-gain-search-v2";
inline constexpr const char* eventRankVersion="candidate-review-rank-v3";
inline constexpr const char* eventFeaturesVersion="event-boundary-features-v1";
inline constexpr const char* eventCoverageVersion="event-treatment-coverage-v1";
inline constexpr const char* eventScheduleVersion="event-assessment-schedule-v1";
// Initial evaluation hypotheses, not a musical-onset standard or user DSP controls.
struct EventFeatureDefinition
{
    int hopMs=10;
    static constexpr int hpHz=20, lowHz=160, fastMs=20, lowMs=80, backgroundMs=250;
    static constexpr int floorDb=-80, fastRiseDb=6, lowRiseDb=3, fastPersistMs=20, lowPersistMs=40;
    static constexpr int agreementMs=80, groupMs=1000, lastAttackMs=100;
    static constexpr int decayDb=6, decayPersistMs=200, returnDb=3, returnPersistMs=300, contextMs=3000;
    static constexpr size_t maximumBytes=10*1024*1024, maximumGroupsPerRoot=64;
    juce::String validationError()const;
    juce::var encode()const;
    static EventFeatureDefinition decode(const juce::var&);
    juce::String digest()const;
};
}
