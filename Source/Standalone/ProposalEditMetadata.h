#pragma once
#include "ManualEditHistory.h"
#include "AppliedProposalExplanation.h"

namespace levelscope::standalone
{
// Optional disk-journal supplement for proposal Apply/recovery. Documentary
// history and processing provenance only: never a measurement or live authority.
struct ProposalEditMetadata
{
    proposal::SectionMap provenance;
    std::shared_ptr<const AppliedProposalExplanation> explanation;
    std::optional<proposal::ProposalPolicy> policy;
    juce::String binding,legacyText;
    juce::MemoryBlock bytes() const
    {
        juce::MemoryOutputStream out;out.writeInt(1);
        const auto map=provenance.valid()?ManualEditState{provenance,{},0}.bytes():juce::MemoryBlock();
        out.writeInt(int(map.getSize()));out.write(map.getData(),map.getSize());
        const auto history=explanation?encodeAppliedExplanation(*explanation):std::vector<uint8_t>();
        out.writeInt(int(history.size()));out.write(history.data(),history.size());
        out.writeString(policy?juce::JSON::toString(policy->encode(),true):juce::String());
        out.writeString(binding);out.writeString(explanation?juce::String():legacyText);
        if(out.getDataSize()>20*1024*1024)throw std::runtime_error("Proposal Undo metadata exceeds its storage bound");
        return out.getMemoryBlock();
    }
    static std::optional<ProposalEditMetadata> read(const juce::MemoryBlock& data)
    {
        try
        {
            juce::MemoryInputStream in(data,false);if(in.readInt()!=1)return {};
            ProposalEditMetadata out;auto n=in.readInt();if(n<0||n>8*1024*1024)return {};
            if(n){juce::MemoryBlock map;if(in.readIntoMemoryBlock(map,size_t(n))!=size_t(n))return {};auto parsed=ManualEditState::read(map);if(!parsed)return {};out.provenance=parsed->map;}
            n=in.readInt();if(n<0||n>int(AppliedProposalExplanation::maximumBytes))return {};
            if(n){std::vector<uint8_t> bytes(size_t(n),0);if(in.read(bytes.data(),n)!=n)return {};JobControl job;out.explanation=decodeAppliedExplanation(bytes,job);if(!out.explanation)return {};}
            const auto policy=in.readString();if(policy.isNotEmpty()){out.policy=proposal::ProposalPolicy::decode(juce::JSON::parse(policy));if(out.policy->validationError().isNotEmpty())return {};}
            out.binding=in.readString();out.legacyText=in.readString();
            if(in.getPosition()!=int64_t(data.getSize()))return {};return out;
        }
        catch(...){return {};}
    }
};
}
