#pragma once
#include "OutputDeviceSettings.h"
#include <map>

namespace levelscope::standalone
{
struct MonitoringRoutes
{
    juce::StringArray outputs; // complete ordered interface output list, not only routed outputs
    std::map<ChannelRole,int> roles; // -1 is an intentional mute, including all-muted configurations
};
struct MonitoringRestore
{
    OutputRouting routing;
    juce::String explanation;
    bool matched = false;
};
// Local to a Mac. Never part of a programme, measurement context or export.
class MonitoringSettings
{
public:
    explicit MonitoringSettings(const juce::File& connectionSettingsFile);
    static juce::String identityKey(const OutputDeviceIdentity&);
    static juce::String layoutKey(const ChannelLayout&);
    static ChannelRole speakerRole(ChannelRole);
    std::optional<double> volume(const OutputDeviceIdentity&) const;
    MonitoringRestore restore(const OutputDeviceIdentity&,const juce::StringArray&,const ChannelLayout&) const;
    bool rememberVolume(const OutputDeviceIdentity&,double db);
    bool rememberRouting(const OutputDeviceIdentity&,const juce::StringArray&,const ChannelLayout&,const OutputRouting&);
    const juce::String& warning() const { return notice; }
    const juce::File& file() const { return path; }
private:
    struct Interface { OutputDeviceIdentity identity; std::optional<double> db; std::map<juce::String,MonitoringRoutes> layouts; };
    std::map<juce::String,Interface> interfaces;
    juce::File path;
    juce::String machine,notice;
    Interface* entry(const OutputDeviceIdentity&);
    bool save();
};
}
