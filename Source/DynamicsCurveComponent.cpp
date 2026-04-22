#include "DynamicsCurveComponent.h"
#include "PluginProcessor.h"
#include "Core/Processing/Modules/MultiThresholdDynamicsParamIDs.h"
// [BEGIN UI-CURVE-UPWARD-INCLUDE-GAINLAW]
#include "Core/Processing/DSP/UpwardGainLaw.h"
// [END UI-CURVE-UPWARD-INCLUDE-GAINLAW]

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

// [BEGIN UI-CURVE-DOWN-DRAG-IMPL]
bool DynamicsCurveComponent::getDownwardInteractionGeometry (juce::Rectangle<float>& plotOut,
                                                             float& t2XOut,
                                                             float& t3XOut) const
{
    if (kind != CurveKind::downward)
        return false;

    auto bounds = getLocalBounds().toFloat();
    auto r = bounds.reduced (10.0f);
    if (r.getWidth() < 90.0f || r.getHeight() < 70.0f)
        return false;

    // [BEGIN UI-CURVE-TEXTSCALE-GEOMETRY]
    auto bottomArea = r.removeFromBottom (28.0f);
    auto rightArea  = r.removeFromRight (32.0f);
    juce::ignoreUnused (bottomArea, rightArea);
    // [END UI-CURVE-TEXTSCALE-GEOMETRY]

    auto plot = r;
    if (plot.getWidth() < 60.0f || plot.getHeight() < 40.0f)
        return false;

    auto& apvts = processor.getAPVTS();
    using namespace levelscope::mtdm;

    auto loadParam = [&] (const char* id, float fallback) -> float
    {
        if (auto* a = apvts.getRawParameterValue (id))
            return a->load (std::memory_order_relaxed);
        return fallback;
    };

    const float t2 = loadParam (ParamIDs::t2Lufs, Defaults::t2Lufs);
    const float t3 = loadParam (ParamIDs::t3Lufs, Defaults::t3Lufs);

    // [BEGIN UI-CURVE-XRANGE-NEG90]
    constexpr float xMin = -90.0f;
    constexpr float xMax =   0.0f;
    // [END UI-CURVE-XRANGE-NEG90]

    auto mapX = [&] (float x) -> float
    {
        const float n = (x - xMin) / (xMax - xMin);
        return plot.getX() + juce::jlimit (0.0f, 1.0f, n) * plot.getWidth();
    };

    plotOut = plot;
    t2XOut = mapX (t2);
    t3XOut = mapX (t3);
    return true;
}

void DynamicsCurveComponent::updateMouseCursorForThresholds (juce::Point<float> pos)
{
    if (kind != CurveKind::downward)
        return;

    juce::Rectangle<float> plot;
    float t2X = 0.0f, t3X = 0.0f;

    if (! getDownwardInteractionGeometry (plot, t2X, t3X))
    {
        const int oldHover = hoverThresholdIndex;
        hoverThresholdIndex = -1;
        if (oldHover != hoverThresholdIndex)
            repaint();

        if (! thresholdDragging)
            setMouseCursor (juce::MouseCursor::NormalCursor);
        return;
    }

    const auto hit = plot.expanded (6.0f, 0.0f);

    int newHover = -1;
    if (hit.contains (pos))
    {
        const float d2 = std::abs (pos.x - t2X);
        const float d3 = std::abs (pos.x - t3X);

        const bool hit2 = (d2 <= 6.0f);
        const bool hit3 = (d3 <= 6.0f);

        if (hit2 && hit3)
            newHover = (d2 <= d3 ? 2 : 3);
        else if (hit2)
            newHover = 2;
        else if (hit3)
            newHover = 3;
    }

    if (hoverThresholdIndex != newHover)
    {
        hoverThresholdIndex = newHover;
        repaint();
    }

    if (! thresholdDragging)
        setMouseCursor (hoverThresholdIndex >= 0 ? juce::MouseCursor::LeftRightResizeCursor
                                                 : juce::MouseCursor::NormalCursor);
}

void DynamicsCurveComponent::computeOrderedThresholdsWithPush (int changedIndex,
                                                               float newValueLufs,
                                                               float outVals[4]) const noexcept
{
    if (changedIndex < 0 || changedIndex > 3)
        return;

    auto& apvts = processor.getAPVTS();
    using namespace levelscope::mtdm;

    const char* ids[4] =
    {
        ParamIDs::t0Lufs,
        ParamIDs::t1Lufs,
        ParamIDs::t2Lufs,
        ParamIDs::t3Lufs
    };

    juce::RangedAudioParameter* params[4] =
    {
        apvts.getParameter (ids[0]),
        apvts.getParameter (ids[1]),
        apvts.getParameter (ids[2]),
        apvts.getParameter (ids[3])
    };

    auto clampSnap = [&] (int idx, float v) -> float
    {
        if (auto* p = params[idx])
        {
            const auto r = p->getNormalisableRange();
            v = juce::jlimit ((float) r.start, (float) r.end, v);
            v = r.snapToLegalValue (v);
        }
        return v;
    };

    for (int i = 0; i < 4; ++i)
    {
        float v = 0.0f;
        if (auto* a = apvts.getRawParameterValue (ids[i]))
            v = a->load (std::memory_order_relaxed);
        else
            v = (i == 0 ? Defaults::t0Lufs
                 : i == 1 ? Defaults::t1Lufs
                 : i == 2 ? Defaults::t2Lufs
                          : Defaults::t3Lufs);

        outVals[i] = clampSnap (i, v);
    }

    outVals[changedIndex] = clampSnap (changedIndex, newValueLufs);

    for (int i = changedIndex + 1; i < 4; ++i)
    {
        const float minAllowed = outVals[i - 1] + kThreshMinGapLu;
        if (outVals[i] < minAllowed)
            outVals[i] = minAllowed;

        outVals[i] = clampSnap (i, outVals[i]);
    }

    for (int i = changedIndex - 1; i >= 0; --i)
    {
        const float maxAllowed = outVals[i + 1] - kThreshMinGapLu;
        if (outVals[i] > maxAllowed)
            outVals[i] = maxAllowed;

        outVals[i] = clampSnap (i, outVals[i]);
    }
}

void DynamicsCurveComponent::applyThresholdValuesDuringDrag (const float newVals[4])
{
    auto& apvts = processor.getAPVTS();
    using namespace levelscope::mtdm;

    const char* ids[4] =
    {
        ParamIDs::t0Lufs,
        ParamIDs::t1Lufs,
        ParamIDs::t2Lufs,
        ParamIDs::t3Lufs
    };

    for (int i = 0; i < 4; ++i)
    {
        auto* p = apvts.getParameter (ids[i]);
        auto* a = apvts.getRawParameterValue (ids[i]);
        if (p == nullptr || a == nullptr)
            continue;

        const float oldV = a->load (std::memory_order_relaxed);
        const float newV = newVals[i];
        if (std::abs (newV - oldV) < 1.0e-6f)
            continue;

        if (! thresholdGestureActive[(size_t) i])
        {
            p->beginChangeGesture();
            thresholdGestureActive[(size_t) i] = true;
        }

        p->setValueNotifyingHost (juce::jlimit (0.0f, 1.0f, p->convertTo0to1 (newV)));
    }
}

void DynamicsCurveComponent::endAllThresholdGestures()
{
    auto& apvts = processor.getAPVTS();
    using namespace levelscope::mtdm;

    const char* ids[4] =
    {
        ParamIDs::t0Lufs,
        ParamIDs::t1Lufs,
        ParamIDs::t2Lufs,
        ParamIDs::t3Lufs
    };

    for (int i = 0; i < 4; ++i)
    {
        if (! thresholdGestureActive[(size_t) i])
            continue;

        if (auto* p = apvts.getParameter (ids[i]))
            p->endChangeGesture();

        thresholdGestureActive[(size_t) i] = false;
    }
}

void DynamicsCurveComponent::mouseMove (const juce::MouseEvent& e)
{
    updateMouseCursorForThresholds (e.position);
}

void DynamicsCurveComponent::mouseExit (const juce::MouseEvent& e)
{
    juce::ignoreUnused (e);

    if (! thresholdDragging)
        setMouseCursor (juce::MouseCursor::NormalCursor);

    if (hoverThresholdIndex != -1)
    {
        hoverThresholdIndex = -1;
        repaint();
    }
}

void DynamicsCurveComponent::mouseDown (const juce::MouseEvent& e)
{
    if (kind != CurveKind::downward || ! e.mods.isLeftButtonDown())
        return;

    juce::Rectangle<float> plot;
    float t2X = 0.0f, t3X = 0.0f;

    if (! getDownwardInteractionGeometry (plot, t2X, t3X))
        return;

    const auto hit = plot.expanded (6.0f, 0.0f);
    if (! hit.contains (e.position))
        return;

    const float d2 = std::abs (e.position.x - t2X);
    const float d3 = std::abs (e.position.x - t3X);

    const bool hit2 = (d2 <= 6.0f);
    const bool hit3 = (d3 <= 6.0f);

    int hitIndex = -1;
    if (hit2 && hit3)      hitIndex = (d2 <= d3 ? 2 : 3);
    else if (hit2)         hitIndex = 2;
    else if (hit3)         hitIndex = 3;

    if (hitIndex < 0)
        return;

    auto& apvts = processor.getAPVTS();
    using namespace levelscope::mtdm;
    const char* ids[4] =
    {
        ParamIDs::t0Lufs,
        ParamIDs::t1Lufs,
        ParamIDs::t2Lufs,
        ParamIDs::t3Lufs
    };

    if (auto* p = apvts.getParameter (ids[hitIndex]))
    {
        p->beginChangeGesture();
        thresholdGestureActive[(size_t) hitIndex] = true;
        thresholdDragging = true;
        activeThresholdIndex = hitIndex;
        setMouseCursor (juce::MouseCursor::LeftRightResizeCursor);
        repaint();
    }
}

void DynamicsCurveComponent::mouseDrag (const juce::MouseEvent& e)
{
    if (kind != CurveKind::downward || ! thresholdDragging || activeThresholdIndex < 0)
        return;

    juce::Rectangle<float> plot;
    float t2X = 0.0f, t3X = 0.0f;

    if (! getDownwardInteractionGeometry (plot, t2X, t3X))
        return;

    // [BEGIN UI-CURVE-XRANGE-NEG90-DRAG]
    constexpr float xMin = -90.0f;
    constexpr float xMax =   0.0f;
    // [END UI-CURVE-XRANGE-NEG90-DRAG]

    const float x = juce::jlimit (plot.getX(), plot.getRight(), e.position.x);
    const float n = (plot.getWidth() > 1.0f ? (x - plot.getX()) / plot.getWidth() : 0.0f);
    const float lufs = xMin + juce::jlimit (0.0f, 1.0f, n) * (xMax - xMin);

    float newVals[4] {};
    computeOrderedThresholdsWithPush (activeThresholdIndex, lufs, newVals);
    applyThresholdValuesDuringDrag (newVals);

    hoverThresholdIndex = activeThresholdIndex;
    repaint();
}

void DynamicsCurveComponent::mouseUp (const juce::MouseEvent& e)
{
    juce::ignoreUnused (e);

    if (! thresholdDragging)
        return;

    endAllThresholdGestures();

    thresholdDragging = false;
    activeThresholdIndex = -1;
    setMouseCursor (juce::MouseCursor::NormalCursor);
    repaint();
}
// [END UI-CURVE-DOWN-DRAG-IMPL]

void DynamicsCurveComponent::paint (juce::Graphics& g)
{
    auto bounds = getLocalBounds().toFloat();
    g.setColour (juce::Colours::black.withMultipliedAlpha (0.14f));
    g.fillRoundedRectangle (bounds, 6.0f);

    g.setColour (juce::Colours::white.withMultipliedAlpha (0.10f));
    g.drawRoundedRectangle (bounds, 6.0f, 1.0f);

    auto r = bounds.reduced (10.0f);
    if (r.getWidth() < 90.0f || r.getHeight() < 70.0f)
        return;

    // [BEGIN UI-CURVE-UPWARD-PAINT]
    if (kind == CurveKind::upwardConceptual)
    {
        // [BEGIN UI-CURVE-TEXTSCALE-UPWARD-AREAS]
        auto topArea    = r.removeFromTop (20.0f);
        auto bottomArea = r.removeFromBottom (28.0f);
        auto rightArea  = r.removeFromRight (32.0f);
        // [END UI-CURVE-TEXTSCALE-UPWARD-AREAS]
        auto plot       = r;

        if (plot.getWidth() < 60.0f || plot.getHeight() < 40.0f)
            return;

        auto& apvts = processor.getAPVTS();
        using namespace levelscope::mtdm;

        const auto loadParam = [&] (const char* id, float fallback) -> float
        {
            if (auto* a = apvts.getRawParameterValue (id))
                return a->load (std::memory_order_relaxed);
            return fallback;
        };

        const float t0        = loadParam (ParamIDs::t0Lufs,          Defaults::t0Lufs);
        const float t1        = loadParam (ParamIDs::t1Lufs,          Defaults::t1Lufs);
        const float amount01  = loadParam (ParamIDs::sucAmount01,     Defaults::sucAmount01);
        const float maxBoost  = loadParam (ParamIDs::sucMaxBoostDb,   Defaults::sucMaxBoostDb);
        const float curve01   = loadParam (ParamIDs::sucCurve,        Defaults::sucCurve);
        const float lowKnee   = loadParam (ParamIDs::sucLowKneeDb,    Defaults::sucLowKneeDb);
        const float highKnee  = loadParam (ParamIDs::sucHighKneeDb,   Defaults::sucHighKneeDb);
        const float calTrimDb = loadParam (ParamIDs::sucCalTrimDb,    Defaults::sucCalTrimDb);

        const int curveTypeChoice =
            (int) std::lround (loadParam (ParamIDs::sucCurveTypeChoice, (float) Defaults::sucCurveTypeChoice));

        const int upwardModeChoice =
            (int) std::lround (loadParam (ParamIDs::upwardModeChoice, (float) Defaults::upwardModeChoice));

        // [BEGIN UI-CURVE-UPWARD-AUTO-XRANGE-NEG90]
        const float tLo = juce::jmin (t0, t1);
        const float tHi = juce::jmax (t0, t1);

        const float safeLowKnee  = juce::jmax (0.0f, lowKnee);
        const float safeHighKnee = juce::jmax (0.0f, highKnee);

        const float lowerKneeStart = tLo - safeLowKnee;

        // Desired window: centered-ish around the zone, but must include:
        // - lower knee start (below T0)
        // - some space above T1 (so you see the return to 0 boost)
        float xMin = lowerKneeStart - 12.0f;
        float xMax = tHi + 18.0f;

        // Ensure a minimum readable span
        const float minSpan = 60.0f;
        if ((xMax - xMin) < minSpan)
        {
            const float mid = 0.5f * (tLo + tHi);
            xMin = mid - 0.5f * minSpan;
            xMax = mid + 0.5f * minSpan;

            // re-ensure knee + headroom constraints
            xMin = juce::jmin (xMin, lowerKneeStart - 6.0f);
            xMax = juce::jmax (xMax, tHi + 12.0f);
        }

        // Clamp to sensible loudness domain (extended to -90)
        xMin = juce::jmax (-90.0f, xMin);
        xMax = juce::jmin (  0.0f, xMax);

        // Final safety: avoid degenerate spans after clamping
        if ((xMax - xMin) < 40.0f)
            xMin = juce::jmax (-90.0f, xMax - 40.0f);
        // [END UI-CURVE-UPWARD-AUTO-XRANGE-NEG90]

        const float safeMaxBoost = juce::jlimit (0.0f, 24.0f, maxBoost);
        const float yMin = 0.0f;
        const float yMax = juce::jmax (1.0f, safeMaxBoost);

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

        // Grid (x like other curves; y is Boost increasing upward)
        g.setColour (juce::Colours::white.withMultipliedAlpha (0.08f));
        for (float xTick : { -60.0f, -48.0f, -36.0f, -24.0f, -12.0f, 0.0f })
            g.drawVerticalLine ((int) std::round (mapX (xTick)), plot.getY(), plot.getBottom());

        // Y ticks: use a stable set + always include yMax
        auto shouldDrawTick = [&] (float v) { return v >= yMin - 1.0e-6f && v <= yMax + 1.0e-6f; };

        juce::Array<float> yTicks;
        for (float v : { 0.0f, 3.0f, 6.0f, 9.0f, 12.0f, 18.0f, 24.0f })
            if (shouldDrawTick (v))
                yTicks.add (v);

        if (yTicks.isEmpty() || std::abs (yTicks.getLast() - yMax) > 0.25f)
            yTicks.add (yMax);

        for (auto v : yTicks)
            g.drawHorizontalLine ((int) std::round (mapY (v)), plot.getX(), plot.getRight());

        // Right-side Boost ruler (0 at bottom, increasing upward)
        g.setFont (14.0f);
        g.setColour (juce::Colours::white.withMultipliedAlpha (0.55f));

        for (auto v : yTicks)
        {
            const float y = mapY (v);
            g.drawLine (plot.getRight(), y, plot.getRight() + 4.0f, y, 1.0f);

            const juce::String s = (v >= 9.95f ? juce::String ((int) std::lround (v))
                                               : juce::String (v, 1));

            g.drawText (s,
                        rightArea.toNearestInt().withY ((int) std::round (y - 10.0f)).withHeight (20),
                        juce::Justification::centredRight, false);
        }

        // Bottom LUFS ruler: ticks/labels + axis title
        const auto bottomTicksArea = bottomArea.removeFromTop (14.0f);
        const auto bottomTitleArea = bottomArea;

        for (float xTick : { -60.0f, -48.0f, -36.0f, -24.0f, -12.0f, 0.0f })
        {
            const float x = mapX (xTick);

            g.setColour (juce::Colours::white.withMultipliedAlpha (0.40f));
            g.drawLine (x, plot.getBottom(), x, plot.getBottom() + 4.0f, 1.0f);

            g.setColour (juce::Colours::white.withMultipliedAlpha (0.55f));
            g.drawText (juce::String ((int) xTick),
                        juce::Rectangle<int> ((int) std::round (x - 18.0f),
                                              (int) bottomTicksArea.getY(),
                                              36,
                                              (int) bottomTicksArea.getHeight()),
                        juce::Justification::centred, false);
        }

        // Shade active band T0..T1
        {
            const float a = juce::jmin (t0, t1);
            const float b = juce::jmax (t0, t1);

            const float xA = mapX (a);
            const float xB = mapX (b);

            g.setColour (juce::Colours::orange.withMultipliedAlpha (0.10f));
            g.fillRect (juce::Rectangle<float> (juce::jmin (xA, xB), plot.getY(),
                                                std::abs (xB - xA), plot.getHeight()));
        }

        // One-sided knee hints (below thresholds)
        {
            const float a = juce::jmin (t0, t1);
            const float b = juce::jmax (t0, t1);

            const float x0A = mapX (a - juce::jmax (0.0f, lowKnee));
            const float x0B = mapX (a);

            const float x1A = mapX (b - juce::jmax (0.0f, highKnee));
            const float x1B = mapX (b);

            // [BEGIN UI-CURVE-UPWARD-KNEE-CONTRAST]
            g.setColour (juce::Colours::orange.withMultipliedAlpha (0.12f));
            // [END UI-CURVE-UPWARD-KNEE-CONTRAST]
            g.fillRect (juce::Rectangle<float> (juce::jmin (x0A, x0B), plot.getY(),
                                                std::abs (x0B - x0A), plot.getHeight()));
            g.fillRect (juce::Rectangle<float> (juce::jmin (x1A, x1B), plot.getY(),
                                                std::abs (x1B - x1A), plot.getHeight()));
        }

        // Threshold markers
        const float t0X = mapX (t0);
        const float t1X = mapX (t1);

        g.setColour (juce::Colours::white.withMultipliedAlpha (0.35f));
        g.drawVerticalLine ((int) std::round (t0X), plot.getY(), plot.getBottom());
        g.drawVerticalLine ((int) std::round (t1X), plot.getY(), plot.getBottom());

        // Sample curve from shared DSP law (DSP-matched)
        levelscope::dsp::UpwardGainLaw::Params gp;
        gp.t0Db        = t0;
        gp.t1Db        = t1;
        gp.lowKneeDb   = juce::jmax (0.0f, lowKnee);
        gp.highKneeDb  = juce::jmax (0.0f, highKnee);
        gp.maxBoostDb  = juce::jmax (0.0f, safeMaxBoost);
        gp.curve01     = juce::jlimit (0.0f, 1.0f, curve01);
        gp.curveType   = (curveTypeChoice == 1
                            ? levelscope::dsp::UpwardGainLaw::CurveType::bell
                            : levelscope::dsp::UpwardGainLaw::CurveType::monotonic);

        juce::Path p;
        bool started = false;

        constexpr int N = 192;
        for (int i = 0; i <= N; ++i)
        {
            const float a = (float) i / (float) N;
            const float x = xMin + a * (xMax - xMin);

            // UI incorporates cal trim as an x-axis shift (conceptual mapping)
            const float xAdjusted = x + calTrimDb;

            const float boostDb =
                levelscope::dsp::UpwardGainLaw::computeUpwardGainDb (xAdjusted, gp, amount01);

            const float px = mapX (x);
            const float py = mapY (juce::jlimit (yMin, yMax, boostDb));

            if (! started) { p.startNewSubPath (px, py); started = true; }
            else           { p.lineTo (px, py); }
        }

        g.setColour (juce::Colours::orange.withMultipliedAlpha (0.95f));
        g.strokePath (p, juce::PathStrokeType (2.0f));

        // Upward metering snapshot -> side indicators (current + hold)
        {
            const auto met = processor.getUpwardMeteringSnapshot();
            const float cur  = juce::jlimit (yMin, yMax, juce::jmax (0.0f, met.boostDbCurrent));
            const float hold = juce::jlimit (yMin, yMax, juce::jmax (0.0f, met.boostDbHold));

            const float yCur  = mapY (cur);
            const float yHold = mapY (hold);

            g.setColour (juce::Colours::orange.withMultipliedAlpha (0.85f));
            g.drawLine (plot.getRight() - 18.0f, yCur, plot.getRight() - 2.0f, yCur, 1.6f);

            g.setColour (juce::Colours::white.withMultipliedAlpha (0.85f));
            g.drawLine (plot.getRight() - 18.0f, yHold, plot.getRight() - 2.0f, yHold, 1.2f);
        }

        // Top threshold labels
        g.setFont (14.0f);
        g.setColour (juce::Colours::white.withMultipliedAlpha (0.70f));
        g.drawText ("T0", (int) t0X - 12, (int) plot.getY() + 2, 24, 12, juce::Justification::centred);
        g.drawText ("T1", (int) t1X - 12, (int) plot.getY() + 16, 24, 12, juce::Justification::centred);

        // Top row: curve type + amount + mode tag
        const juce::String typeStr = (curveTypeChoice == 1 ? "Bell" : "Monotonic");
        const juce::String modeStr = (upwardModeChoice == 1 ? "Broadband" : "Spectral");

        g.setFont (14.0f);
        g.setColour (juce::Colours::white.withMultipliedAlpha (0.70f));
        // [BEGIN UI-CURVE-UPWARD-TOPLABELS-NO-OVERLAP]
        const juce::String topLeft =
            modeStr + "   " + typeStr + "   Amt " + juce::String (juce::jlimit (0.0f, 1.0f, amount01), 2);

        g.setFont (14.0f);
        g.setColour (juce::Colours::white.withMultipliedAlpha (0.72f));
        g.drawFittedText (topLeft,
                          topArea.toNearestInt(),
                          juce::Justification::centredLeft,
                          1);
        // [END UI-CURVE-UPWARD-TOPLABELS-NO-OVERLAP]

        // Footer axis labels + info
        g.setFont (14.0f);
        const int axisLabelW = 60;

        g.setColour (juce::Colours::white.withMultipliedAlpha (0.58f));
        g.drawText ("LUFS",
                    juce::Rectangle<int> ((int) bottomTitleArea.getX(),
                                          (int) bottomTitleArea.getY(),
                                          axisLabelW,
                                          (int) bottomTitleArea.getHeight()),
                    juce::Justification::centredLeft, false);

        const juce::String info =
            "Max " + juce::String (safeMaxBoost, 1) + " dB"
            + "   Knees " + juce::String (juce::jmax (0.0f, lowKnee), 1)
            + "/" + juce::String (juce::jmax (0.0f, highKnee), 1) + " dB"
            + "   Trim " + juce::String (calTrimDb, 1) + " dB";

        g.setColour (juce::Colours::white.withMultipliedAlpha (0.55f));
        g.drawFittedText (info,
                          juce::Rectangle<int> ((int) bottomTitleArea.getX() + axisLabelW,
                                                (int) bottomTitleArea.getY(),
                                                (int) bottomTitleArea.getWidth() - axisLabelW,
                                                (int) bottomTitleArea.getHeight()),
                          juce::Justification::centredLeft, 1);

        g.setColour (juce::Colours::white.withMultipliedAlpha (0.58f));
        g.drawText ("Boost",
                    juce::Rectangle<int> ((int) rightArea.getX(), (int) topArea.getY(),
                                          (int) rightArea.getWidth(), 12),
                    juce::Justification::centredRight, false);

        return;
    }
    // [END UI-CURVE-UPWARD-PAINT]

    if (kind != CurveKind::downward)
    {
        g.setColour (juce::Colours::white.withMultipliedAlpha (0.45f));
        g.setFont (12.0f);
        g.drawFittedText ("Curve pending", r.toNearestInt(), juce::Justification::centred, 1);
        return;
    }

    auto topArea    = r.removeFromTop (20.0f);
    auto bottomArea = r.removeFromBottom (28.0f);
    auto rightArea  = r.removeFromRight (32.0f);
    auto plot       = r;

    if (plot.getWidth() < 60.0f || plot.getHeight() < 40.0f)
        return;

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
        return plot.getX() + juce::jlimit (0.0f, 1.0f, n) * plot.getWidth();
    };

    auto mapY = [&] (float y) -> float
    {
        const float n = (y - yMin) / (yMax - yMin);
        return plot.getBottom() - juce::jlimit (0.0f, 1.0f, n) * plot.getHeight();
    };

    // Background grid
    g.setColour (juce::Colours::white.withMultipliedAlpha (0.08f));
    for (float xTick : { -60.0f, -48.0f, -36.0f, -24.0f, -12.0f, 0.0f })
        g.drawVerticalLine ((int) std::round (mapX (xTick)), plot.getY(), plot.getBottom());

    for (float yTick : { 0.0f, 6.0f, 12.0f, 18.0f, 24.0f })
        g.drawHorizontalLine ((int) std::round (mapY (yTick)), plot.getX(), plot.getRight());

    // Right-side GR ruler (displayed as 0 at top, 24 at bottom to match downward meter motion)
    g.setFont (14.0f);
    g.setColour (juce::Colours::white.withMultipliedAlpha (0.55f));
    for (float yTickDisplay : { 0.0f, 6.0f, 12.0f, 18.0f, 24.0f })
    {
        const float n = juce::jlimit (0.0f, 1.0f, yTickDisplay / yMax);
        const float y = plot.getY() + n * plot.getHeight();

        g.drawLine (plot.getRight(), y, plot.getRight() + 4.0f, y, 1.0f);
        g.drawText (juce::String ((int) yTickDisplay),
                    rightArea.toNearestInt().withY ((int) std::round (y - 10.0f)).withHeight (20),
                    juce::Justification::centredRight, false);
    }

    // Bottom LUFS ruler: top line = ticks/labels, bottom line = axis title
    const auto bottomTicksArea = bottomArea.removeFromTop (14.0f);
    const auto bottomTitleArea = bottomArea;

    for (float xTick : { -60.0f, -48.0f, -36.0f, -24.0f, -12.0f, 0.0f })
    {
        const float x = mapX (xTick);
        g.setColour (juce::Colours::white.withMultipliedAlpha (0.40f));
        g.drawLine (x, plot.getBottom(), x, plot.getBottom() + 4.0f, 1.0f);

        g.setColour (juce::Colours::white.withMultipliedAlpha (0.55f));
        g.drawText (juce::String ((int) xTick),
                    juce::Rectangle<int> ((int) std::round (x - 18.0f),
                                          (int) bottomTicksArea.getY(),
                                          36,
                                          (int) bottomTicksArea.getHeight()),
                    juce::Justification::centred, false);
    }

    // Shade the T2..T3 zone
    {
        const float xA = mapX (juce::jmin (t2, t3));
        const float xB = mapX (juce::jmax (t2, t3));

        g.setColour (juce::Colours::deepskyblue.withMultipliedAlpha (0.08f));
        g.fillRect (juce::Rectangle<float> (xA, plot.getY(), xB - xA, plot.getHeight()));
    }

    // Knee band around T2
    {
        const float xA = mapX (t2 - 0.5f * knee);
        const float xB = mapX (t2 + 0.5f * knee);

        // [BEGIN UI-CURVE-DOWN-KNEE-CONTRAST]
        g.setColour (juce::Colours::deepskyblue.withMultipliedAlpha (0.12f));
        // [END UI-CURVE-DOWN-KNEE-CONTRAST]
        g.fillRect (juce::Rectangle<float> (xA, plot.getY(), xB - xA, plot.getHeight()));
    }

    // Threshold markers
    const float t2X = mapX (t2);
    const float t3X = mapX (t3);

    g.setColour (juce::Colours::white.withMultipliedAlpha (0.35f));
    g.drawVerticalLine ((int) std::round (t2X), plot.getY(), plot.getBottom());
    g.drawVerticalLine ((int) std::round (t3X), plot.getY(), plot.getBottom());

    auto drawHandleCue = [&] (float x, int idx, float y)
    {
        const bool active = (thresholdDragging && activeThresholdIndex == idx);
        const bool hover  = (! thresholdDragging && hoverThresholdIndex == idx);

        g.setColour (juce::Colours::white.withMultipliedAlpha (active ? 0.95f : (hover ? 0.82f : 0.62f)));
        g.fillRoundedRectangle (juce::Rectangle<float> (x - 3.0f, y, 6.0f, 12.0f), 2.0f);
    };

    drawHandleCue (t2X, 2, plot.getY() + 16.0f);
    drawHandleCue (t3X, 3, plot.getY() + 30.0f);

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

    // Unified downward display mapping:
    // 0 at top, more GR lower on screen.
    auto mapMeterY = [&] (float gr) -> float
    {
        const float n = juce::jlimit (0.0f, 1.0f, gr / yMax);
        return plot.getY() + n * plot.getHeight();
    };

    const int numSteps = juce::jlimit (80, 240, (int) std::round (plot.getWidth()));
    for (int i = 0; i <= numSteps; ++i)
    {
        const float a = (float) i / (float) numSteps;
        const float x = xMin + a * (xMax - xMin);
        const float gr = gainReductionForInput (x);

        const float px = mapX (x);
        const float py = mapMeterY (gr);

        if (! started) { p.startNewSubPath (px, py); started = true; }
        else           { p.lineTo (px, py); }
    }

    g.setColour (juce::Colours::deepskyblue.withMultipliedAlpha (0.95f));
    g.strokePath (p, juce::PathStrokeType (2.0f));

    // Existing GR snapshot -> side indicator + true current x/y dot
    const auto met = processor.getDownwardMeteringSnapshot();
    const float grCur  = juce::jmax (0.0f, met.grDbCurrent);
    const float grHold = juce::jmax (0.0f, met.grDbHold);

    const float yCur  = mapMeterY (grCur);
    const float yHold = mapMeterY (grHold);

    g.setColour (juce::Colours::deepskyblue.withMultipliedAlpha (0.85f));
    g.drawLine (plot.getRight() - 18.0f, yCur, plot.getRight() - 2.0f, yCur, 1.6f);

    g.setColour (juce::Colours::white.withMultipliedAlpha (0.85f));
    g.drawLine (plot.getRight() - 18.0f, yHold, plot.getRight() - 2.0f, yHold, 1.2f);

    // True current x/y operating point:
    // x = detector loudness currently feeding the downward stage
    // y = current GR magnitude in the same display mapping as the curve/ruler
    if (std::isfinite (met.detectorLufsCurrent) && std::isfinite (met.grDbCurrent)
        && met.detectorLufsCurrent > -199.0f)
    {
        const float xDot = juce::jlimit (plot.getX(), plot.getRight(),
                                         mapX (juce::jlimit (xMin, xMax, met.detectorLufsCurrent)));
        const float yDot = juce::jlimit (plot.getY(), plot.getBottom(),
                                         mapMeterY (juce::jlimit (yMin, yMax, met.grDbCurrent)));

        g.setColour (juce::Colours::black.withMultipliedAlpha (0.55f));
        g.fillEllipse (xDot - 4.5f, yDot - 4.5f, 9.0f, 9.0f);

        g.setColour (juce::Colours::deepskyblue.withMultipliedAlpha (0.95f));
        g.fillEllipse (xDot - 3.0f, yDot - 3.0f, 6.0f, 6.0f);

        g.setColour (juce::Colours::white.withMultipliedAlpha (0.95f));
        g.drawEllipse (xDot - 4.0f, yDot - 4.0f, 8.0f, 8.0f, 1.0f);
    }

    // Top threshold labels
    g.setFont (14.0f);
    g.setColour (juce::Colours::white.withMultipliedAlpha (0.70f));
    g.drawText ("T2", (int) t2X - 12, (int) plot.getY() + 2, 24, 12, juce::Justification::centred);
    g.drawText ("T3", (int) t3X - 12, (int) plot.getY() + 16, 24, 12, juce::Justification::centred);

    // Footer info and labels
    g.setFont (14.0f);

    const int axisLabelW = 60;

    g.setColour (juce::Colours::white.withMultipliedAlpha (0.58f));
    g.drawText ("LUFS",
                juce::Rectangle<int> ((int) bottomTitleArea.getX(),
                                      (int) bottomTitleArea.getY(),
                                      axisLabelW,
                                      (int) bottomTitleArea.getHeight()),
                juce::Justification::centredLeft, false);

    const juce::String info = "Ratio " + juce::String (safeRatio, 2)
                            + "   Knee " + juce::String (safeKnee, 1) + " dB";
    g.drawFittedText (info,
                      juce::Rectangle<int> ((int) bottomTitleArea.getX() + axisLabelW,
                                            (int) bottomTitleArea.getY(),
                                            (int) bottomTitleArea.getWidth() - axisLabelW,
                                            (int) bottomTitleArea.getHeight()),
                      juce::Justification::centredLeft, 1);

    g.drawText ("GR",
                juce::Rectangle<int> ((int) rightArea.getX(), (int) topArea.getY(), (int) rightArea.getWidth(), 12),
                juce::Justification::centredRight, false);
}