#include "LevelerCurveComponent.h"
#include "PluginProcessor.h"
#include "Core/Processing/Modules/LevelerParamIDs.h"

#include <cmath>

LevelerCurveComponent::LevelerCurveComponent (LevelScopeAudioProcessor& p)
    : processor (p)
{
    setOpaque (false);
    startTimerHz (15);
}

void LevelerCurveComponent::timerCallback()
{
    repaint();
}

void LevelerCurveComponent::paint (juce::Graphics& g)
{
    auto bounds = getLocalBounds().toFloat();
    g.setColour (juce::Colours::black.withMultipliedAlpha (0.14f));
    g.fillRoundedRectangle (bounds, 6.0f);

    g.setColour (juce::Colours::white.withMultipliedAlpha (0.10f));
    g.drawRoundedRectangle (bounds, 6.0f, 1.0f);

    auto r = bounds.reduced (10.0f);
    if (r.getWidth() < 90.0f || r.getHeight() < 60.0f)
        return;

    auto& apvts = processor.getAPVTS();
    using namespace levelscope::lvlr;

    const auto loadParam = [&] (const char* id, float fallback) -> float
    {
        if (auto* a = apvts.getRawParameterValue (id))
            return a->load (std::memory_order_relaxed);
        return fallback;
    };

    const float target   = loadParam (ParamIDs::targetLufs,       Defaults::targetLufs);
    const float maxBoost = loadParam (ParamIDs::maxBoostDb,       Defaults::maxBoostDb);
    const float maxCut   = loadParam (ParamIDs::maxCutDb,         Defaults::maxCutDb);
    const int measChoice = (int) std::lround (loadParam (ParamIDs::measChoice, (float) Defaults::measChoice));
    const int modeChoice = (int) std::lround (loadParam (ParamIDs::modeChoice, (float) Defaults::modeChoice));

    const int controlMode = (int) std::lround (loadParam (ParamIDs::controlModeChoice,
                                                          (float) Defaults::controlModeChoice));
    const float hostGainDb = loadParam (ParamIDs::hostGainDb, Defaults::hostGainDb);
    const bool captureArmed = (loadParam (ParamIDs::captureToHost01, Defaults::captureToHost01) >= 0.5f);

    const bool hostGainMode = (controlMode == 1);

    constexpr float xMin = levelscope::lvlr::Ranges::targetMinLufs;
    constexpr float xMax = levelscope::lvlr::Ranges::targetMaxLufs;
    constexpr float yMin = levelscope::lvlr::Ranges::hostGainMinDb;
    constexpr float yMax = levelscope::lvlr::Ranges::hostGainMaxDb;

    auto mapX = [&] (float x) -> float
    {
        const float n = (x - xMin) / (xMax - xMin);
        return r.getX() + juce::jlimit (0.0f, 1.0f, n) * r.getWidth();
    };

    auto mapY = [&] (float y) -> float
    {
        const float n = (y - yMin) / (yMax - yMin);
        return r.getBottom() - juce::jlimit (0.0f, 1.0f, n) * r.getHeight();
    };

    // Grid
    g.setColour (juce::Colours::white.withMultipliedAlpha (0.08f));
    for (float xTick : { -48.0f, -42.0f, -36.0f, -30.0f, -24.0f, -18.0f, -12.0f })
        g.drawVerticalLine ((int) std::round (mapX (xTick)), r.getY(), r.getBottom());

    for (float yTick : { -24.0f, -12.0f, 0.0f, 12.0f, 24.0f })
        g.drawHorizontalLine ((int) std::round (mapY (yTick)), r.getX(), r.getRight());

    // Dim/normal overlay alpha depending on control mode
    const float activeAlpha = hostGainMode ? 0.22f : 0.90f;
    const float dimAlpha    = hostGainMode ? 0.10f : 0.75f;

    // Clamp lines
    g.setColour (juce::Colours::limegreen.withMultipliedAlpha (activeAlpha));
    g.drawHorizontalLine ((int) std::round (mapY ( maxBoost)), r.getX(), r.getRight());

    g.setColour (juce::Colours::deepskyblue.withMultipliedAlpha (activeAlpha));
    g.drawHorizontalLine ((int) std::round (mapY (-maxCut)), r.getX(), r.getRight());

    // Target marker
    g.setColour (juce::Colours::white.withMultipliedAlpha (hostGainMode ? 0.18f : 0.45f));
    g.drawVerticalLine ((int) std::round (mapX (target)), r.getY(), r.getBottom());

    // Internal conceptual mapping
    {
        juce::Path p;
        bool started = false;

        const int numSteps = juce::jlimit (80, 240, (int) std::round (r.getWidth()));
        for (int i = 0; i <= numSteps; ++i)
        {
            const float a = (float) i / (float) numSteps;
            const float x = xMin + a * (xMax - xMin);

            float y = target - x;
            y = juce::jlimit (-maxCut, maxBoost, y);

            const float px = mapX (x);
            const float py = mapY (y);

            if (! started) { p.startNewSubPath (px, py); started = true; }
            else           { p.lineTo (px, py); }
        }

        g.setColour (juce::Colours::white.withMultipliedAlpha (hostGainMode ? 0.22f : 0.92f));
        g.strokePath (p, juce::PathStrokeType (2.0f));
    }

    // Host gain overlay if active
    if (hostGainMode)
    {
        const float yHost = mapY (hostGainDb);
        g.setColour (juce::Colours::orange.withMultipliedAlpha (0.88f));
        g.drawHorizontalLine ((int) std::round (yHost), r.getX(), r.getRight());

        g.setFont (12.0f);
        g.drawFittedText ("Host Gain active",
                          ((int) r.getX()) + 2, (int) r.getY() + 2,
                          (int) r.getWidth() - 4, 14,
                          juce::Justification::topLeft, 1);

        g.setColour (juce::Colours::orange.withMultipliedAlpha (0.70f));
        g.drawFittedText ("Host " + juce::String (hostGainDb, 1) + " dB",
                          ((int) r.getX()) + 2, (int) r.getY() + 18,
                          (int) r.getWidth() - 4, 14,
                          juce::Justification::topLeft, 1);
    }

    // Actual applied gain marker(s)
    const auto met = processor.getLevelerMeteringSnapshot();
    const float yCur  = mapY (juce::jlimit (yMin, yMax, met.gainDbCurrent));
    const float yHold = mapY (juce::jlimit (yMin, yMax, met.gainDbHold));

    g.setColour (met.gainDbCurrent >= 0.0f
                    ? juce::Colours::limegreen.withMultipliedAlpha (0.92f)
                    : juce::Colours::deepskyblue.withMultipliedAlpha (0.88f));
    g.drawLine (r.getRight() - 18.0f, yCur, r.getRight() - 2.0f, yCur, 1.6f);

    g.setColour (juce::Colours::white.withMultipliedAlpha (0.85f));
    g.drawLine (r.getRight() - 18.0f, yHold, r.getRight() - 2.0f, yHold, 1.2f);

    // Internal-mode capture badge
    if (! hostGainMode && captureArmed)
    {
        auto badge = juce::Rectangle<float> (r.getRight() - 58.0f, r.getY() + 4.0f, 54.0f, 16.0f);

        g.setColour (juce::Colours::red.withMultipliedAlpha (0.80f));
        g.fillRoundedRectangle (badge, 4.0f);

        g.setColour (juce::Colours::white.withMultipliedAlpha (0.95f));
        g.setFont (11.0f);
        g.drawFittedText ("CAPTURE", badge.toNearestInt(), juce::Justification::centred, 1);
    }

    // Labels
    const char* measText = "Auto";
    if (measChoice == 1) measText = "Momentary";
    else if (measChoice == 2) measText = "Short-term";

    const char* modeText = (modeChoice == 1 ? "Learn-Hold" : "Adaptive");

    g.setFont (11.0f);
    g.setColour (juce::Colours::white.withMultipliedAlpha (dimAlpha));
    g.drawText ("Measured LUFS", (int) r.getX(), (int) r.getBottom() - 14, 90, 14, juce::Justification::left);
    g.drawText ("Gain", (int) r.getRight() - 30, (int) r.getY(), 28, 14, juce::Justification::right);

    g.setColour (juce::Colours::white.withMultipliedAlpha (0.62f));
    const juce::String info = juce::String (measText) + "   " + juce::String (modeText);
    g.drawFittedText (info,
                      ((int) r.getX()) + 2, (int) r.getBottom() - 30,
                      (int) r.getWidth() - 4, 14,
                      juce::Justification::bottomLeft, 1);

    g.setColour (juce::Colours::white.withMultipliedAlpha (hostGainMode ? 0.45f : 0.68f));
    g.drawText ("T", (int) mapX (target) - 8, (int) r.getY() + 2, 16, 12, juce::Justification::centred);
}