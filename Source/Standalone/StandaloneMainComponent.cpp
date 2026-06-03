#include "StandaloneMainComponent.h"

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-WORKFLOW-SHELL]

StandaloneMainComponent::StandaloneMainComponent()
{
    configureLabels();
    configureTargetProfileBox();
    configureNavigationButtons();
    configureSourceControls();

    configureReadOnlyBox (sourceStateBox);
    configureReadOnlyBox (targetStateBox);
    configureReadOnlyBox (currentStateBox);
    configureReadOnlyBox (pageDetailBox);

    addAndMakeVisible (sourceStateBox);
    addAndMakeVisible (targetStateBox);
    addAndMakeVisible (currentStateBox);
    addAndMakeVisible (pageDetailBox);

    refreshFromSession();

    setSize (1120, 760);
}

void StandaloneMainComponent::configureLabels()
{
    titleLabel.setText ("LevelScope Standalone", juce::dontSendNotification);
    titleLabel.setJustificationType (juce::Justification::centredLeft);
    titleLabel.setFont (juce::Font (28.0f, juce::Font::bold));
    titleLabel.setColour (juce::Label::textColourId, juce::Colours::white);
    addAndMakeVisible (titleLabel);

    subtitleLabel.setText ("Workflow shell only: import, measurement, proposal, render, and export are not wired yet.",
                           juce::dontSendNotification);
    subtitleLabel.setJustificationType (juce::Justification::centredLeft);
    subtitleLabel.setColour (juce::Label::textColourId, juce::Colour (0xffc6ccd7));
    addAndMakeVisible (subtitleLabel);

    workflowStateLabel.setJustificationType (juce::Justification::centredLeft);
    workflowStateLabel.setColour (juce::Label::textColourId, juce::Colour (0xffc6ccd7));
    addAndMakeVisible (workflowStateLabel);

    targetProfileLabel.setText ("Target profile:", juce::dontSendNotification);
    targetProfileLabel.setJustificationType (juce::Justification::centredRight);
    targetProfileLabel.setColour (juce::Label::textColourId, juce::Colour (0xffc6ccd7));
    addAndMakeVisible (targetProfileLabel);
}

void StandaloneMainComponent::configureTargetProfileBox()
{
    for (int i = 0; i < StandaloneSessionModel::numTargetProfiles; ++i)
        targetProfileBox.addItem (StandaloneSessionModel::getTargetProfile (i).displayName, i + 1);

    targetProfileBox.setSelectedId (session.selectedTargetProfileIndex + 1, juce::dontSendNotification);

    targetProfileBox.onChange = [this]
    {
        const auto selectedId = targetProfileBox.getSelectedId();

        if (selectedId > 0)
            session.setSelectedTargetProfileIndex (selectedId - 1);

        refreshFromSession();
    };

    addAndMakeVisible (targetProfileBox);
}

void StandaloneMainComponent::configureNavigationButtons()
{
    importButton.onClick = [this] { setPage (WorkflowPage::importSource); };
    analyzeButton.onClick = [this] { setPage (WorkflowPage::analyze); };
    currentStateButton.onClick = [this] { setPage (WorkflowPage::currentState); };
    exportButton.onClick = [this] { setPage (WorkflowPage::exportResults); };

    addAndMakeVisible (importButton);
    addAndMakeVisible (analyzeButton);
    addAndMakeVisible (currentStateButton);
    addAndMakeVisible (exportButton);
}

void StandaloneMainComponent::configureSourceControls()
{
    chooseSourceButton.onClick = [this] { chooseSourceFile(); };
    clearSourceButton.onClick = [this] { clearSourceFile(); };

    measureSourceButton.setEnabled (false);

    addAndMakeVisible (chooseSourceButton);
    addAndMakeVisible (clearSourceButton);
    addAndMakeVisible (measureSourceButton);
}

void StandaloneMainComponent::configureReadOnlyBox (juce::TextEditor& box)
{
    box.setMultiLine (true);
    box.setReadOnly (true);
    box.setCaretVisible (false);
    box.setScrollbarsShown (true);

    box.setColour (juce::TextEditor::backgroundColourId, juce::Colour (0xff20242b));
    box.setColour (juce::TextEditor::textColourId, juce::Colour (0xffeef1f7));
    box.setColour (juce::TextEditor::outlineColourId, juce::Colour (0xff3a3f4a));
    box.setColour (juce::TextEditor::focusedOutlineColourId, juce::Colour (0xff6f7788));
}

void StandaloneMainComponent::chooseSourceFile()
{
    sourceFileChooser = std::make_unique<juce::FileChooser> (
        "Choose an audio source file",
        juce::File(),
        getSourceFileWildcard());

    juce::Component::SafePointer<StandaloneMainComponent> safeThis (this);

    sourceFileChooser->launchAsync (
        juce::FileBrowserComponent::openMode | juce::FileBrowserComponent::canSelectFiles,
        [safeThis] (const juce::FileChooser& chooser)
        {
            if (safeThis == nullptr)
                return;

            const auto chosenFile = chooser.getResult();

            if (! chosenFile.existsAsFile())
                return;

            safeThis->session.setSourceFile (chosenFile);
            safeThis->session.selectedPage = WorkflowPage::analyze;
            safeThis->refreshFromSession();
        });
}

void StandaloneMainComponent::clearSourceFile()
{
    session.clearSourceFile();
    refreshFromSession();
}

void StandaloneMainComponent::setPage (WorkflowPage page)
{
    session.selectedPage = page;
    refreshFromSession();
}

void StandaloneMainComponent::refreshFromSession()
{
    workflowStateLabel.setText ("Workflow page: " + getPageName (session.selectedPage)
                                    + "   |   Proposal engine: not implemented",
                                juce::dontSendNotification);

    targetProfileBox.setSelectedId (session.selectedTargetProfileIndex + 1, juce::dontSendNotification);

    clearSourceButton.setEnabled (session.sourceDocument.hasSource);

    sourceStateBox.setText (buildSourceStateText(), false);
    targetStateBox.setText (buildTargetStateText(), false);
    currentStateBox.setText (buildCurrentStateText(), false);
    pageDetailBox.setText (buildPageDetailText(), false);
}

juce::String StandaloneMainComponent::getSourceFileWildcard()
{
    return "*.wav;*.wave;*.aif;*.aiff;*.flac;*.mp3;*.m4a";
}

juce::String StandaloneMainComponent::getPageName (WorkflowPage page) const
{
    switch (page)
    {
        case WorkflowPage::importSource:  return "Import";
        case WorkflowPage::analyze:       return "Analyze";
        case WorkflowPage::currentState:  return "Current State";
        case WorkflowPage::exportResults: return "Export";
        default:                          return "Unknown";
    }
}

juce::String StandaloneMainComponent::buildSourceStateText() const
{
    juce::String text;

    text << "SOURCE\n\n";
    text << "File: " << session.sourceDocument.displayName << "\n";
    text << "Import status: " << session.sourceDocument.importStatus << "\n";
    text << "Measurement status: " << session.source.statusText << "\n\n";

    if (session.sourceDocument.hasSource)
    {
        text << "Extension: " << session.sourceDocument.extension << "\n";
        text << "Size: " << session.sourceDocument.fileSizeBytes << " bytes\n";
        text << "Path: " << session.sourceDocument.fullPath << "\n\n";
    }
    else
    {
        text << "Extension: --\n";
        text << "Size: --\n";
        text << "Path: --\n\n";
    }

    if (session.source.hasMeasurement)
    {
        text << "Integrated: " << session.source.integratedLufs << " LUFS\n";
        text << "LRA: " << session.source.loudnessRangeLu << " LU\n";
        text << "True Peak: " << session.source.truePeakDbtp << " dBTP\n";
    }
    else
    {
        text << "Integrated: --\n";
        text << "LRA: --\n";
        text << "True Peak: --\n";
    }

    return text;
}

juce::String StandaloneMainComponent::buildTargetStateText() const
{
    const auto profile = session.getSelectedTargetProfile();

    juce::String text;

    text << "TARGET\n\n";
    text << "Profile: " << profile.displayName << "\n";
    text << "ID: " << profile.id << "\n\n";
    text << profile.description << "\n\n";
    text << "Compliance values: not loaded\n";
    text << "Profile validation: not implemented\n";

    return text;
}

juce::String StandaloneMainComponent::buildCurrentStateText() const
{
    juce::String text;

    text << "CURRENT STATE\n\n";
    text << "Status: " << session.currentState.statusText << "\n\n";

    if (session.currentState.hasMeasurement)
    {
        text << "Integrated: " << session.currentState.integratedLufs << " LUFS\n";
        text << "LRA: " << session.currentState.loudnessRangeLu << " LU\n";
        text << "True Peak: " << session.currentState.truePeakDbtp << " dBTP\n";
    }
    else
    {
        text << "Integrated: --\n";
        text << "LRA: --\n";
        text << "True Peak: --\n";
    }

    text << "\nRe-measure flow: not implemented\n";

    return text;
}

juce::String StandaloneMainComponent::buildPageDetailText() const
{
    juce::String text;

    text << getPageName (session.selectedPage).toUpperCase() << " VIEW\n\n";

    switch (session.selectedPage)
    {
        case WorkflowPage::importSource:
            text << "Source import scaffold.\n\n";
            text << "Implemented now:\n";
            text << "- Choose an audio source file.\n";
            text << "- Store file name, path, extension, and size in standalone session state.\n";
            text << "- Clear the selected source.\n\n";
            text << "Accepted picker patterns:\n";
            text << getSourceFileWildcard() << "\n\n";

            if (session.sourceDocument.hasSource)
            {
                text << "Selected source:\n";
                text << session.sourceDocument.displayName << "\n\n";
                text << "Next step remains measurement implementation.\n";
            }
            else
            {
                text << "No source selected yet.\n";
            }

            break;

        case WorkflowPage::analyze:
            text << "Measurement view scaffold.\n\n";

            if (session.sourceDocument.hasSource)
            {
                text << "Source is ready for future measurement:\n";
                text << session.sourceDocument.displayName << "\n\n";
                text << "The Measure Source button is intentionally disabled in this patch.\n\n";
                text << "Future implementation:\n";
                text << "- Decode the selected source file.\n";
                text << "- Run loudness/peak measurement.\n";
                text << "- Write measured values into Source state.\n";
            }
            else
            {
                text << "Import a source file before measurement can be implemented.\n";
            }

            break;

        case WorkflowPage::currentState:
            text << "Current State scaffold.\n\n";
            text << "Current State remains separate from Source.\n\n";
            text << "Future implementation:\n";
            text << "- After processing changes are possible, re-measure Current State.\n";
            text << "- Keep real measured Current State separate from any future proposal prediction.\n";
            text << "- Do not add proposal calculations yet.\n";
            break;

        case WorkflowPage::exportResults:
            text << "Export readiness scaffold.\n\n";
            text << "Export remains blocked until future steps provide:\n";
            text << "- imported source\n";
            text << "- measured source state\n";
            text << "- selected real target profile\n";
            text << "- measured current state\n\n";
            text << "No render or export logic is implemented in this patch.\n";
            break;

        default:
            text << "Unknown workflow page.\n";
            break;
    }

    text << "\n\nGuardrail: this standalone shell must not depend on PluginProcessor, PluginEditor, or APVTS.";

    return text;
}

void StandaloneMainComponent::paint (juce::Graphics& g)
{
    g.fillAll (juce::Colour (0xff15171c));

    g.setColour (juce::Colour (0xff2a2d35));
    g.drawRoundedRectangle (getLocalBounds().toFloat().reduced (12.0f), 10.0f, 1.0f);
}

void StandaloneMainComponent::resized()
{
    auto bounds = getLocalBounds().reduced (24);

    titleLabel.setBounds (bounds.removeFromTop (40));
    subtitleLabel.setBounds (bounds.removeFromTop (28));

    bounds.removeFromTop (10);

    auto targetRow = bounds.removeFromTop (32);
    targetProfileLabel.setBounds (targetRow.removeFromLeft (120));
    targetProfileBox.setBounds (targetRow.removeFromLeft (360));
    workflowStateLabel.setBounds (targetRow.reduced (16, 0));

    bounds.removeFromTop (14);

    auto navRow = bounds.removeFromTop (34);
    importButton.setBounds (navRow.removeFromLeft (120));
    navRow.removeFromLeft (8);
    analyzeButton.setBounds (navRow.removeFromLeft (120));
    navRow.removeFromLeft (8);
    currentStateButton.setBounds (navRow.removeFromLeft (140));
    navRow.removeFromLeft (8);
    exportButton.setBounds (navRow.removeFromLeft (120));

    bounds.removeFromTop (10);

    auto sourceActionRow = bounds.removeFromTop (34);
    chooseSourceButton.setBounds (sourceActionRow.removeFromLeft (190));
    sourceActionRow.removeFromLeft (8);
    clearSourceButton.setBounds (sourceActionRow.removeFromLeft (120));
    sourceActionRow.removeFromLeft (8);
    measureSourceButton.setBounds (sourceActionRow.removeFromLeft (240));

    bounds.removeFromTop (14);

    auto summaryArea = bounds.removeFromTop (190);
    const auto gap = 12;
    const auto columnWidth = (summaryArea.getWidth() - (2 * gap)) / 3;

    sourceStateBox.setBounds (summaryArea.removeFromLeft (columnWidth));
    summaryArea.removeFromLeft (gap);
    targetStateBox.setBounds (summaryArea.removeFromLeft (columnWidth));
    summaryArea.removeFromLeft (gap);
    currentStateBox.setBounds (summaryArea);

    bounds.removeFromTop (14);

    pageDetailBox.setBounds (bounds);
}

// [END LS-STANDALONE-WORKFLOW-SHELL]
}