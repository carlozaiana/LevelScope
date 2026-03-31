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
    if (r.getWidth() < 110.0f || r.getHeight() < 80.0f)
        return;

    auto topArea    = r.removeFromTop (28.0f);
    auto bottomArea = r.removeFromBottom (18.0f);
    auto rightArea  = r.removeFromRight (32.0f);
    auto plot       = r;

    if (plot.getWidth() < 60.0f || plot.getHeight() < 40.0f)
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
        return plot.getX() + juce::jlimit (0.0f, 1.0f, n) * plot.getWidth();
    };

    auto mapY = [&] (float y) -> float
    {
        const float n = (y - yMin) / (yMax - yMin);
        return plot.getBottom() - juce::jlimit (0.0f, 1.0f, n) * plot.getHeight();
    };

    // Grid
    g.setColour (juce::Colours::white.withMultipliedAlpha (0.08f));
    for (float xTick : { -48.0f, -42.0f, -36.0f, -30.0f, -24.0f, -18.0f, -12.0f })
        g.drawVerticalLine ((int) std::round (mapX (xTick)), plot.getY(), plot.getBottom());

    for (float yTick : { -24.0f, -12.0f, 0.0f, 12.0f, 24.0f })
        g.drawHorizontalLine ((int) std::round (mapY (yTick)), plot.getX(), plot.getRight());

    // Right-side dB ruler
    g.setFont (10.0f);
    g.setColour (juce::Colours::white.withMultipliedAlpha (0.55f));
    for (float yTick : { -24.0f, -12.0f, 0.0f, 12.0f, 24.0f })
    {
        const float y = mapY (yTick);
        g.drawLine (plot.getRight(), y, plot.getRight() + 4.0f, y, 1.0f);
        g.drawText (juce::String ((int) yTick),
                    rightArea.toNearestInt().withY ((int) std::round (y - 7.0f)).withHeight (14),
                    juce::Justification::centredRight, false);
    }

    // Bottom LUFS ruler
    for (float xTick : { -48.0f, -36.0f, -24.0f, -12.0f })
    {
        const float x = mapX (xTick);
        g.setColour (juce::Colours::white.withMultipliedAlpha (0.40f));
        g.drawLine (x, plot.getBottom(), x, plot.getBottom() + 4.0f, 1.0f);

        g.setColour (juce::Colours::white.withMultipliedAlpha (0.55f));
        g.drawText (juce::String ((int) xTick),
                    juce::Rectangle<int> ((int) std::round (x - 18.0f),
                                          (int) bottomArea.getY(),
                                          36,
                                          (int) bottomArea.getHeight()),
                    juce::Justification::centred, false);
    }

    // Dim/normal overlay alpha depending on control mode
    const float activeAlpha = hostGainMode ? 0.22f : 0.90f;
    const float dimAlpha    = hostGainMode ? 0.10f : 0.75f;

    // Clamp lines
    g.setColour (juce::Colours::limegreen.withMultipliedAlpha (activeAlpha));
    g.drawHorizontalLine ((int) std::round (mapY ( maxBoost)), plot.getX(), plot.getRight());

    g.setColour (juce::Colours::deepskyblue.withMultipliedAlpha (activeAlpha));
    g.drawHorizontalLine ((int) std::round (mapY (-maxCut)), plot.getX(), plot.getRight());

    // Target marker
    g.setColour (juce::Colours::white.withMultipliedAlpha (hostGainMode ? 0.18f : 0.45f));
    g.drawVerticalLine ((int) std::round (mapX (target)), plot.getY(), plot.getBottom());

    // Internal conceptual mapping
    {
        juce::Path p;
        bool started = false;

        const int numSteps = juce::jlimit (80, 240, (int) std::round (plot.getWidth()));
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
        g.drawHorizontalLine ((int) std::round (yHost), plot.getX(), plot.getRight());
    }

    // Actual applied gain marker(s)
    const auto met = processor.getLevelerMeteringSnapshot();
    const float yCur  = mapY (juce::jlimit (yMin, yMax, met.gainDbCurrent));
    const float yHold = mapY (juce::jlimit (yMin, yMax, met.gainDbHold));

    g.setColour (met.gainDbCurrent >= 0.0f
                    ? juce::Colours::limegreen.withMultipliedAlpha (0.92f)
                    : juce::Colours::deepskyblue.withMultipliedAlpha (0.88f));
    g.drawLine (plot.getRight() - 18.0f, yCur, plot.getRight() - 2.0f, yCur, 1.6f);

    g.setColour (juce::Colours::white.withMultipliedAlpha (0.85f));
    g.drawLine (plot.getRight() - 18.0f, yHold, plot.getRight() - 2.0f, yHold, 1.2f);

    // Top status area
    const char* measText = "Auto";
    if (measChoice == 1) measText = "Momentary";
    else if (measChoice == 2) measText = "Short-term";

    const char* modeText = (modeChoice == 1 ? "Learn-Hold" : "Adaptive");

    g.setFont (11.0f);

    if (hostGainMode)
    {
        g.setColour (juce::Colours::orange.withMultipliedAlpha (0.90f));
        g.drawFittedText ("Host Gain active", topArea.removeFromTop (14).toNearestInt(),
                          juce::Justification::centredLeft, 1);

        g.setColour (juce::Colours::orange.withMultipliedAlpha (0.72f));
        g.drawFittedText ("Host " + juce::String (hostGainDb, 1) + " dB",
                          topArea.toNearestInt(),
                          juce::Justification::centredLeft, 1);
    }
    else
    {
        auto line = topArea.removeFromTop (14);
        g.setColour (juce::Colours::white.withMultipliedAlpha (0.65f));
        g.drawFittedText (juce::String (measText) + "   " + juce::String (modeText),
                          line.toNearestInt(),
                          juce::Justification::centredLeft, 1);

        if (captureArmed)
        {
            auto badge = juce::Rectangle<float> (line.getRight() - 54.0f, line.getY(), 54.0f, 14.0f);
            g.setColour (juce::Colours::red.withMultipliedAlpha (0.80f));
            g.fillRoundedRectangle (badge, 4.0f);

            g.setColour (juce::Colours::white.withMultipliedAlpha (0.95f));
            g.setFont (10.0f);
            g.drawFittedText ("CAPTURE", badge.toNearestInt(), juce::Justification::centred, 1);
        }
    }

    // Axis labels / target label
    g.setFont (10.0f);
    g.setColour (juce::Colours::white.withMultipliedAlpha (dimAlpha));
    g.drawText ("Measured LUFS",
                juce::Rectangle<int> ((int) plot.getX(), (int) bottomArea.getY(), 90, (int) bottomArea.getHeight()),
                juce::Justification::centredLeft, false);

    g.drawText ("Gain",
                juce::Rectangle<int> ((int) rightArea.getX(), (int) plot.getY(), (int) rightArea.getWidth(), 12),
                juce::Justification::centredRight, false);

    g.setColour (juce::Colours::white.withMultipliedAlpha (hostGainMode ? 0.45f : 0.68f));
    g.drawText ("T", (int) mapX (target) - 8, (int) plot.getY() + 2, 16, 12, juce::Justification::centred);
}