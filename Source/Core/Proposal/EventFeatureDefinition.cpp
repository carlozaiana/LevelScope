#include "EventFeatureDefinition.h"
#include <juce_cryptography/juce_cryptography.h>
#include <stdexcept>
namespace levelscope::proposal
{
juce::String EventFeatureDefinition::validationError()const
{return hopMs==5||hopMs==10||hopMs==20?juce::String():"Unsupported explicit event cadence";}
juce::var EventFeatureDefinition::encode()const
{
    if(validationError().isNotEmpty())throw std::invalid_argument(validationError().toStdString());
    auto v=juce::var(new juce::DynamicObject());auto* o=v.getDynamicObject();
    o->setProperty("version",innerEntrances?innerFeaturesVersion:eventFeaturesVersion);o->setProperty("filters","continuous-first-order-hp-lp-v1");
    o->setProperty("channels","non-lfe-energy-mean-and-role-witness-v1");
    o->setProperty("silence","separate-floor-onset-and-return-v1");
    for(auto p:{std::pair<const char*,int>{"hopMs",hopMs},{"hpHz",hpHz},{"lowHz",lowHz},
        {"fastMs",fastMs},{"lowMs",lowMs},{"backgroundMs",backgroundMs},{"floorDb",floorDb},
        {"fastRiseDb",fastRiseDb},{"lowRiseDb",lowRiseDb},{"fastPersistMs",fastPersistMs},
        {"lowPersistMs",lowPersistMs},{"agreementMs",agreementMs},{"groupMs",groupMs},
        {"lastAttackMs",lastAttackMs},{"decayDb",decayDb},{"decayPersistMs",decayPersistMs},
        {"returnDb",returnDb},{"returnPersistMs",returnPersistMs},{"contextMs",contextMs}})
        o->setProperty(p.first,p.second);
    if(innerEntrances)
    {
        o->setProperty("predecessors",predecessors);
        o->setProperty("referenceStepDb",referenceStepDb);
        o->setProperty("momentaryStepLu",momentaryStepLu);
        o->setProperty("neighbourhoodMs",neighbourhoodMs);
        o->setProperty("comparableDb",comparableDb);
        o->setProperty("recentReferences",recentReferences);
        o->setProperty("attribution","included-window-excess-weighted-bin-share-v1");
        o->setProperty("pairing","one-live-inner-directed-edges-v1");
        o->setProperty("missing","retain-outer-report-incomplete-v1");
        o->setProperty("rolePath","same-role-tentative-v1");
        o->setProperty("eventLocations",6);
        o->setProperty("innerProxySlots",8);
    }
    return v;
}
EventFeatureDefinition EventFeatureDefinition::decode(const juce::var& v)
{
    if(!v.isObject()||!v["hopMs"].isInt())throw std::invalid_argument("Missing explicit event feature definition");
    EventFeatureDefinition out;out.hopMs=int(v["hopMs"]);out.innerEntrances=v["version"].toString()==innerFeaturesVersion;const auto expected=out.encode();
    if(v.getDynamicObject()->getProperties().size()!=expected.getDynamicObject()->getProperties().size())throw std::invalid_argument("Unknown event feature fields");
    for(const auto& p:expected.getDynamicObject()->getProperties())
        if(!v.hasProperty(p.name)||!v[p.name].equalsWithSameType(p.value))throw std::invalid_argument("Unsupported event feature hypothesis");
    return out;
}
juce::String EventFeatureDefinition::digest()const
{auto text=juce::JSON::toString(encode(),true).toStdString();return juce::SHA256(text.data(),text.size()).toHexString();}
}
