#pragma once
#include "StandaloneSessionModel.h"
#include <juce_gui_basics/juce_gui_basics.h>
namespace levelscope::standalone
{
// Bounded page view prevents 10,000-row proposal text from flooding a UI editor.
class ProposalInspectionPanel final : public juce::Component
{
public:
    explicit ProposalInspectionPanel(juce::String);
    void resized() override;
    int pageCount() const { return int(pages.size()); }
    juce::String allText() const;
    void nextPage();
private:
    std::vector<juce::String> pages;int page=0;
    juce::TextEditor view;juce::Label location;
    juce::TextButton previous{"Previous page"},next{"Next page"};
    void refresh();
};
class ProposalPreferencesPanel final : public juce::Component
{
public:
    ProposalPreferencesPanel(StandaloneSessionModel&,std::function<void()>);
    void resized() override;
    static std::optional<int> exactUnits(juce::String,bool tenths);
private:
    StandaloneSessionModel& session;std::function<void()> changed;juce::String runtime,baseDigest;
    juce::ToggleButton unrestricted{"No region-count restriction"};
    std::array<juce::Label,5> labels;
    std::array<juce::TextEditor,5> fields;
    juce::TextButton use{"Use these limits"},cancel{"Cancel edits"},defaults{"Draft defaults"};
    juce::TextEditor status;
    void load(const proposal::ProposalPolicy&);
};
}
