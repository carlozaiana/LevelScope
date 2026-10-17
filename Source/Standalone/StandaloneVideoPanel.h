#pragma once
#include "StandaloneVideo.h"
#include "StandaloneSessionModel.h"
#include "StandaloneAudition.h"
#include <juce_gui_basics/juce_gui_basics.h>

namespace levelscope::standalone
{
class StandaloneVideoPanel final : public juce::Component, private juce::Timer
{
public:
    StandaloneVideoPanel (StandaloneSessionModel&, StandaloneAudition&);
    ~StandaloneVideoPanel() override;
    bool prerollReady() const { return player.prerollReady(); }
    void synchronize() { timerCallback(); }
    void paint (juce::Graphics&) override;
    void resized() override;
private:
    StandaloneSessionModel& session;
    StandaloneAudition& audio;
    video::Player player;
    video::Player::View view;
    juce::TextButton attach { "Attach video..." }, clear { "Clear video" }, apply { "Apply offset" };
    juce::ComboBox streams;
    juce::TextEditor offset;
    juce::Label offsetLabel, status;
    juce::TextEditor metadata;
    juce::Rectangle<int> pictureBounds;
    std::unique_ptr<juce::FileChooser> chooser;
    const video::Asset* displayedMetadata = nullptr;
    juce::String actionError;
    void timerCallback() override;
};
}
