#pragma once

// [BEGIN MTDM-PARAM-GROUPS-HEADER]
// Non-RT metadata only.
// Purpose:
// - Provide a stable grouping of MTDM parameters for UI organization and future policy/standalone tooling.
// - Does not affect DSP, persistence, or parameter IDs.
// - Safe to use from message thread / UI thread / offline analysis.
//
// IMPORTANT:
// - Parameter IDs are the compatibility contract. UI labels/grouping may change freely,
//   but IDs in MultiThresholdDynamicsParamIDs.h must remain stable.
// [END MTDM-PARAM-GROUPS-HEADER]

#include <juce_core/juce_core.h>
#include <vector>

#include "MultiThresholdDynamicsParamIDs.h"

namespace levelscope::mtdm
{
    // [BEGIN MTDM-PARAM-GROUPS-DECL]
    struct ParamGroup
    {
        juce::String groupKey;                 // stable-ish internal key (not a persistence contract)
        juce::String displayName;              // user-facing group label (may change)
        std::vector<const char*> paramIDs;     // ParamIDs::* entries (stable)
    };

    // Returns default MTDM parameter groups.
    // Non-RT: may allocate.
    inline std::vector<ParamGroup> createDefaultParamGroups()
    {
        using namespace ParamIDs;

        std::vector<ParamGroup> groups;
        groups.reserve (8);

        // Master / Mode selection
        // [BEGIN MTDM-PARAMGROUPS-ADVANCED-ROUTING-ADD-MC]
        groups.push_back (ParamGroup {
            "mtdm.routing.advanced",
            "Advanced Routing",
            {
                mcPolicyChoice,
                dialogDetectorChoice,
                dialogApplyChoice,

                lfeInDetector,
                lfeInApply
            }
        });
        // [END MTDM-PARAMGROUPS-ADVANCED-ROUTING-ADD-MC]

        // Upward common controls (shared by Spectral + Broadband modes)
        groups.push_back (ParamGroup {
            "mtdm.upward.common",
            "Upward: Common",
            {
                t0Lufs,
                t1Lufs,

                sucAmount01,
                sucMaxBoostDb,
                sucCurve,
                sucCurveTypeChoice,

                sucLowKneeDb,
                sucHighKneeDb,

                sucAttackMs,
                sucReleaseMs,

                sucCalTrimDb
            }
        });

        // Upward spectral advanced controls (Spectral mode only)
        groups.push_back (ParamGroup {
            "mtdm.upward.spectral.advanced",
            "Upward: Spectral Advanced",
            {
                sucFftSizeChoice,
                sucBandsPerOctChoice,
                sucMinFreqHz,
                sucMaxFreqHz
            }
        });

        // Downward compressor (T2–T3)
        groups.push_back (ParamGroup {
            "mtdm.downward",
            "Downward Compressor",
            {
                t2Lufs,
                t3Lufs,

                downEnabled01,
                downRatio,
                downKneeDb,
                downAttackMs,
                downReleaseMs,
                downMakeupDb
            }
        });

        // [BEGIN MTDM-PARAMGROUPS-LIMITER]
        groups.push_back (ParamGroup {
            "mtdm.limiter",
            "Limiter (Safety)",
            {
                limEnabled01,
                limCeilingDb,
                limDriveDb,
                limLookaheadMs,
                limAttackMs,
                limReleaseMs,
                limOversamplingChoice
            }
        });
        // [END MTDM-PARAMGROUPS-LIMITER]

        // Advanced routing / policies
        groups.push_back (ParamGroup {
            "mtdm.routing.advanced",
            "Advanced Routing",
            {
                lfeInDetector,
                lfeInApply
            }
        });

        // [BEGIN MTDM-PARAMGROUPS-ZONE-MUTES]
        groups.push_back (ParamGroup {
            "mtdm.zones",
            "Zones (A/B)",
            {
                zoneSoloChoice,
                zoneUpwardMute01,
                zoneDownwardMute01,
                zoneLimiterMute01,
                zoneUntouchedMute01
            }
        });
        // [END MTDM-PARAMGROUPS-ZONE-MUTES]

        // [BEGIN MTDM-PARAMGROUPS-ZONE-AUDITION]
        groups.push_back (ParamGroup {
            "mtdm.zoneAudition",
            "Zone Audition (T0–T3)",
            {
                zoneAudBelowT0_01,
                zoneAudT0T1_01,
                zoneAudT1T2_01,
                zoneAudT2T3_01,
                zoneAudAboveT3_01
            }
        });
        // [END MTDM-PARAMGROUPS-ZONE-AUDITION]

        return groups;
    }
    // [END MTDM-PARAM-GROUPS-DECL]
} // namespace levelscope::mtdm