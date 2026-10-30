#pragma once
#include <juce_core/juce_core.h>
#include <juce_cryptography/juce_cryptography.h>
#include <optional>
#include <stdexcept>

namespace levelscope::proposal
{
// Session metadata, never a delivery profile or processing-plan component.
struct ProposalPolicy
{
    static constexpr const char* version="section-preferences-v2";
    static constexpr const char* accounting="connected-envelope-support-v2";
    std::optional<int> regionLimit; // absent explicitly means no count restriction
    int durationThousandths=200, additionalTenths=30, cumulativeTenths=120, reserveTenths=5;
    juce::String unavailable; // recover an unknown/invalid policy without inventing defaults
    juce::String validationError() const
    {
        if(unavailable.isNotEmpty())return unavailable;
        if(regionLimit&&(*regionLimit<1||*regionLimit>10000))return "Region limit must be 1-10,000 or No region-count restriction";
        if(durationThousandths<0||durationThousandths>1000)return "Affected duration must be 0-100%";
        if(additionalTenths<0||additionalTenths>60)return "Additional attenuation must be 0-6 dB";
        if(cumulativeTenths<0||cumulativeTenths>120)return "Cumulative local attenuation must be 0-12 dB";
        if(reserveTenths<0||reserveTenths>10||reserveTenths>additionalTenths)return "Reserve must be 0-1 dB and no greater than additional attenuation";
        return {};
    }
    juce::var encode() const
    {
        auto* o=new juce::DynamicObject;
        o->setProperty("version",version);o->setProperty("accounting",accounting);
        o->setProperty("region_limit",regionLimit?juce::var(*regionLimit):juce::var());
        o->setProperty("duration_thousandths",durationThousandths);
        o->setProperty("additional_tenths_db",additionalTenths);o->setProperty("cumulative_tenths_db",cumulativeTenths);
        o->setProperty("reserve_tenths_db",reserveTenths);o->setProperty("unavailable",unavailable);
        return juce::var(o);
    }
    juce::String digest() const
    {
        const auto s=juce::JSON::toString(encode(),true);
        return juce::SHA256(s.toRawUTF8(),s.getNumBytesAsUTF8()).toHexString();
    }
    static ProposalPolicy decode(const juce::var& v)
    {
        ProposalPolicy p;
        try
        {
            const auto* o=v.getDynamicObject();if(!o||o->getProperties().size()!=8)throw std::runtime_error("Invalid proposal policy fields");
            const auto field=[&](const char* key) -> juce::var {if(!o->hasProperty(key))throw std::runtime_error("Missing proposal policy field");return o->getProperty(key);};
            if(field("version").toString()!=version||field("accounting").toString()!=accounting)throw std::runtime_error("Unsupported proposal policy version/accounting");
            const auto integer=[&](const char* key){const auto a=field(key);if(!a.isInt())throw std::runtime_error("Policy values require exact integer units");return int(a);};
            const auto limit=field("region_limit");if(!limit.isVoid()){if(!limit.isInt())throw std::runtime_error("Invalid region limit");p.regionLimit=int(limit);}
            p.durationThousandths=integer("duration_thousandths");p.additionalTenths=integer("additional_tenths_db");
            p.cumulativeTenths=integer("cumulative_tenths_db");p.reserveTenths=integer("reserve_tenths_db");
            if(!field("unavailable").isString()||field("unavailable").toString().length()>512)throw std::runtime_error("Invalid policy recovery status");
            p.unavailable=field("unavailable").toString();if(p.unavailable.isNotEmpty())return p;
            const auto error=p.validationError();if(error.isNotEmpty())throw std::runtime_error(error.toStdString());
        }
        catch(const std::exception& e){p.unavailable="Proposal preferences unavailable: "+juce::String(e.what()).substring(0,400)+". Explicitly replace these limits; audio settings preserved.";}
        return p;
    }
    juce::String text() const
    {
        if(validationError().isNotEmpty())return validationError();
        return juce::String(version)+" / "+accounting+"\nPolicy digest: "+digest()
            +"\n"+(regionLimit?"Maximum changed processing regions: "+juce::String(*regionLimit):"No region-count restriction")
            +"; affected duration "+juce::String(durationThousandths/10.,1)+"% (preservation preference, not broadcast rule).\n"
            +"Additional attenuation "+juce::String(additionalTenths/10.,1)+" dB; cumulative local attenuation "+juce::String(cumulativeTenths/10.,1)
            +" dB; requested reserve "+juce::String(reserveTenths/10.,1)+" dB inside both caps.\n";
    }
};
}
