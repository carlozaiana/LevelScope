#pragma once
#include "StandaloneAudition.h"
#include <optional>

namespace levelscope::standalone
{
struct OutputDeviceIdentity
{
    juce::String type, name, uid;
    bool operator==(const OutputDeviceIdentity& other) const
    { return type == other.type && name == other.name && uid == other.uid; }
};
struct OutputDeviceChoice
{
    OutputDeviceIdentity identity;
    juce::String problem;
};

// Discovery does not open devices. Tests inject the complete catalog and factory.
class OutputDeviceCatalog
{
public:
    virtual ~OutputDeviceCatalog() = default;
    virtual std::vector<OutputDeviceChoice> scan() = 0;
    virtual std::unique_ptr<juce::AudioIODevice> create(const OutputDeviceIdentity&) = 0;
    static std::shared_ptr<OutputDeviceCatalog> system();
};

class OutputDeviceSettings
{
public:
    explicit OutputDeviceSettings(juce::File directory = {}, juce::String machineId = juce::SystemStats::getUniqueDeviceID());
    std::optional<OutputDeviceIdentity> lastOutput;
    bool automaticReconnect = false;
    juce::String warning;
    void save();
    const juce::File& file() const { return path; }
private:
    juce::File path;
    juce::String machineKey;
};

// Message-thread connection intent; independent of programme roles/processing.
// There is no timer-driven retry and this object never starts an audio callback.
class OutputDeviceConnection
{
public:
    OutputDeviceConnection(StandaloneAudition&, std::shared_ptr<OutputDeviceCatalog>, OutputDeviceSettings);
    void rescan();
    void select(const OutputDeviceIdentity&); // explicit intent, may await a known rate
    void connect();                          // explicit retry/intent
    void disconnect();                      // suppress auto until application restart
    void setAutomaticReconnect(bool);
    void programmeChanged(double validatedRate);
    const std::vector<OutputDeviceChoice>& choices() const { return available; }
    const std::optional<OutputDeviceIdentity>& selection() const { return selected; }
    bool automaticReconnect() const { return settings.automaticReconnect; }
    bool automaticSuppressed() const { return suppressed; }
    bool awaitingRate() const { return explicitPending && programmeRate <= 0; }
    juce::String status() const;
    const juce::String& settingsWarning() const { return settings.warning; }
    uint64_t deviceRevision() const { return revision; }
private:
    StandaloneAudition& audio;
    std::shared_ptr<OutputDeviceCatalog> catalog;
    OutputDeviceSettings settings;
    std::vector<OutputDeviceChoice> available;
    std::optional<OutputDeviceIdentity> selected;
    double programmeRate = 0;
    bool suppressed = false, explicitPending = false;
    uint64_t revision = 0;
    juce::String message;
    std::optional<OutputDeviceIdentity> resolve();
    void attempt();
    void describeDisconnected();
};
}
