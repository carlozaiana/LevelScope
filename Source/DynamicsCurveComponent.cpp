#include "DynamicsCurveComponent.h"
#include "PluginProcessor.h"
#include "Core/Processing/Modules/MultiThresholdDynamicsParamIDs.h"

#include <cmath>

DynamicsCurveComponent::DynamicsCurveComponent (LevelScopeAudioProcessor& p, CurveKind k)
    : processor (p), kind (k)
{
    setOpaque (false);
    startTimerHz (15);
}

void DynamicsCurveComponent::timerCallback()
{
    repaint();
}

void DynamicsCurveComponent::paint (juce::Graphics& g)
{
    auto bounds = getLocalBounds().toFloat();
    g.setColour (juce::Colours::black.withMultipliedAlpha (0.14f));
    g.fillRoundedRectangle (bounds, 6.0f);

    g.setColour (juce::Colours::white.withMultipliedAlpha (0.10f));
    g.drawRoundedRectangle (bounds, 6.0f, 1.0f);

    auto r = bounds.reduced (10.0f);
    if (r.getWidth() < 80.0f || r.getHeight() < 60.0f)
        return;

    if (kind != CurveKind::downward)
    {
        g.setColour (juce::Colours::white.withMultipliedAlpha (0.45f));
        g.setFont (12.0f);
        g.drawFittedText ("Curve pending", r.toNearestInt(), juce::Justification::centred, 1);
        return;
    }

    auto& apvts = processor.getAPVTS();
    using namespace levelscope::mtdm;

    const auto loadParam = [&] (const char* id, float fallback) -> float
    {
        if (auto* a = apvts.getRawParameterValue (id))
            return a->load (std::memory_order_relaxed);
        return fallback;
    };

    const float t2    = loadParam (ParamIDs::t2Lufs,      Defaults::t2Lufs);
    const float t3    = loadParam (ParamIDs::t3Lufs,      Defaults::t3Lufs);
    const float ratio = loadParam (ParamIDs::downRatio,   Defaults::downRatio);
    const float knee  = loadParam (ParamIDs::downKneeDb,  Defaults::downKneeDb);

    constexpr float xMin = -60.0f;
    constexpr float xMax =   0.0f;
    constexpr float yMin =   0.0f;
    constexpr float yMax =  24.0f;

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

    // Background grid
    g.setColour (juce::Colours::white.withMultipliedAlpha (0.08f));
    for (float xTick : { -60.0f, -48.0f, -36.0f, -24.0f, -12.0f, 0.0f })
        g.drawVerticalLine ((int) std::round (mapX (xTick)), r.getY(), r.getBottom());

    for (float yTick : { 0.0f, 6.0f, 12.0f, 18.0f, 24.0f })
        g.drawHorizontalLine ((int) std::round (mapY (yTick)), r.getX(), r.getRight());

    // Shade the T2..T3 zone
    {
        const float xA = mapX (juce::jmin (t2, t3));
        const float xB = mapX (juce::jmax (t2, t3));

        g.setColour (juce::Colours::deepskyblue.withMultipliedAlpha (0.08f));
        g.fillRect (juce::Rectangle<float> (xA, r.getY(), xB - xA, r.getHeight()));
    }

    // Knee band around T2
    {
        const float xA = mapX (t2 - 0.5f * knee);
        const float xB = mapX (t2 + 0.5f * knee);

        g.setColour (juce::Colours::white.withMultipliedAlpha (0.06f));
        g.fillRect (juce::Rectangle<float> (xA, r.getY(), xB - xA, r.getHeight()));
    }

    // Threshold markers
    g.setColour (juce::Colours::white.withMultipliedAlpha (0.35f));
    g.drawVerticalLine ((int) std::round (mapX (t2)), r.getY(), r.getBottom());
    g.drawVerticalLine ((int) std::round (mapX (t3)), r.getY(), r.getBottom());

    // Conceptual downward GR curve
    juce::Path p;
    bool started = false;

    const float safeRatio = juce::jmax (1.0f, ratio);
    const float safeKnee  = juce::jmax (0.0f, knee);
    const float T         = t2;
    const float K         = safeKnee;

    auto gainReductionForInput = [&] (float xDb) -> float
    {
        float yOut = xDb;

        if (K > 1.0e-6f)
        {
            const float d = xDb - T;

            if (2.0f * d < -K)
            {
                yOut = xDb;
            }
            else if (2.0f * std::abs (d) <= K)
            {
                yOut = xDb + (1.0f / safeRatio - 1.0f)
                             * (d + K * 0.5f) * (d + K * 0.5f)
                             / (2.0f * K);
            }
            else
            {
                yOut = T + (xDb - T) / safeRatio;
            }
        }
        else
        {
            if (xDb > T)
                yOut = T + (xDb - T) / safeRatio;
        }

        const float gr = juce::jmax (0.0f, xDb - yOut);
        return juce::jlimit (yMin, yMax, gr);
    };

    const int numSteps = juce::jlimit (80, 240, (int) std::round (r.getWidth()));
    for (int i = 0; i <= numSteps; ++i)
    {
        const float a = (float) i / (float) numSteps;
        const float x = xMin + a * (xMax - xMin);
        const float y = gainReductionForInput (x);

        const float px = mapX (x);
        const float py = mapY (y);

        if (! started) { p.startNewSubPath (px, py); started = true; }
        else           { p.lineTo (px, py); }
    }

    g.setColour (juce::Colours::deepskyblue.withMultipliedAlpha (0.95f));
    g.strokePath (p, juce::PathStrokeType (2.0f));

    // Existing GR snapshot -> magnitude-only indicator (not true x/y operating point)
    const auto met = processor.getDownwardMeteringSnapshot();
    const float grCur  = juce::jmax (0.0f, met.grDbCurrent);
    const float grHold = juce::jmax (0.0f, met.grDbHold);

    const float yCur  = mapY (grCur);
    const float yHold = mapY (grHold);

    g.setColour (juce::Colours::deepskyblue.withMultipliedAlpha (0.85f));
    g.drawLine (r.getRight() - 18.0f, yCur, r.getRight() - 2.0f, yCur, 1.6f);

    g.setColour (juce::Colours::white.withMultipliedAlpha (0.85f));
    g.drawLine (r.getRight() - 18.0f, yHold, r.getRight() - 2.0f, yHold, 1.2f);

    // Labels
    g.setFont (11.0f);
    g.setColour (juce::Colours::white.withMultipliedAlpha (0.65f));
    g.drawText ("In (LUFS)", (int) r.getX(), (int) r.getBottom() - 14, 70, 14, juce::Justification::left);
    g.drawText ("GR", (int) r.getRight() - 26, (int) r.getY(), 24, 14, juce::Justification::right);

    g.setColour (juce::Colours::white.withMultipliedAlpha (0.70f));
    g.drawText ("T2", (int) mapX (t2) - 12, (int) r.getY() + 2, 24, 12, juce::Justification::centred);
    g.drawText ("T3", (int) mapX (t3) - 12, (int) r.getY() + 16, 24, 12, juce::Justification::centred);

    g.setColour (juce::Colours::white.withMultipliedAlpha (0.55f));
    const juce::String info = "Ratio " + juce::String (safeRatio, 2)
                            + "   Knee " + juce::String (safeKnee, 1) + " dB";
    g.drawFittedText (info, ((int) r.getX()) + 2, (int) r.getY() + 2, (int) r.getWidth() - 40, 14,
                      juce::Justification::topLeft, 1);
}