#include "VolumeHistoryComponent.h"
#include "PluginProcessor.h"
// [BEGIN MTDM-THRESH-UI-INCLUDE-PARAMIDS]
#include "Core/Processing/Modules/MultiThresholdDynamicsParamIDs.h"
// [END MTDM-THRESH-UI-INCLUDE-PARAMIDS]

//==============================================================================
// VolumeHistoryComponent_Render.cpp
// - paint()
// - cached background drawing
// - playhead line
// - time ruler + dB ruler drawing
// - debug overlay rendering
//==============================================================================
//==============================================================================
// [EDIT-BLOCKS] (token-efficient patch anchors for future chats)
//   - [BEGIN VHC-RENDER-HELPER-FORMAT-TIME] ... [END VHC-RENDER-HELPER-FORMAT-TIME]
//   - [BEGIN VHC-RENDER-PAINT]             ... [END VHC-RENDER-PAINT]
//   - [BEGIN VHC-RENDER-CACHED-BG]         ... [END VHC-RENDER-CACHED-BG]
//   - [BEGIN VHC-RENDER-GEOMETRY]          ... [END VHC-RENDER-GEOMETRY]
//   - [BEGIN VHC-RENDER-RULER-HYST]        ... [END VHC-RENDER-RULER-HYST]
//   - [BEGIN VHC-RENDER-REP-CURVES]        ... [END VHC-RENDER-REP-CURVES]
//   - [BEGIN VHC-RENDER-PARSE-USER-OFFSET] ... [END VHC-RENDER-PARSE-USER-OFFSET]
//==============================================================================

// [BEGIN VHC-RENDER-HELPER-FORMAT-TIME]
// Helper: format seconds as HH:MM:SS (supports negative)
static juce::String formatTimeHMS (double seconds)
{
    const bool neg = (seconds < 0.0);
    seconds = std::abs (seconds);

    const int total = (int) std::floor (seconds + 0.5);
    const int h = total / 3600;
    const int m = (total % 3600) / 60;
    const int s = total % 60;

    const juce::String core = juce::String::formatted ("%02d:%02d:%02d", h, m, s);
    return neg ? "-" + core : core;
}
// [END VHC-RENDER-HELPER-FORMAT-TIME]

//==============================================================================
// DRAW
//==============================================================================

// [BEGIN VHC-RENDER-PAINT]
void VolumeHistoryComponent::paint (juce::Graphics& g)
{
    auto bounds = getLocalBounds().toFloat();
    const int width  = (int) bounds.getWidth();
    const int height = (int) bounds.getHeight();

    // [BEGIN ROLLING-LRA-SPLITTER-MAINPLOT-AREA]
    // Reserve bottom space for time ruler, and optionally for rolling LRA lane.
    const auto timeRuler = getTimeRulerArea();

    int plotBottom = timeRuler.getY();
    if (showRollingLra)
        plotBottom -= rollingLaneHeightPx;

    plotBottom = juce::jlimit (1, juce::jmax (1, height), plotBottom);

    mainPlotArea = juce::Rectangle<float> (0.0f, 0.0f, (float) width, (float) plotBottom);
    // [END ROLLING-LRA-SPLITTER-MAINPLOT-AREA]

    // [VIEW-NAV] Follow mode: follow playhead (so overwrite outside view jumps to that section)
    // Keep playhead at center while following.
    if (followRightEdge)
    {
        if (havePlayheadFrameIndex && zoomX > 1.0e-12)
        {
            const double visibleFrames = (double) width / zoomX;
            constexpr double playheadXFrac = 0.5; // 0.5 = center
            viewRightFrame = (double) playheadFrameIndex + visibleFrames * (1.0 - playheadXFrac);
        }
        else if (haveNowFrameIndex)
        {
            // fallback if playhead is unknown
            viewRightFrame = (double) nowFrameIndex;
        }
    }

    followButton.setToggleState (followRightEdge, juce::dontSendNotification);

    clampViewRightFrame (width);

    //==============================================================================
    // [AUTO-FOLLOW-HYST] Auto re-follow when playhead reaches right edge
    // - Only triggers when Follow is OFF
    // - Requires playhead to be sufficiently left first (arm), then near right edge (trigger)
    //==============================================================================
    if (! followRightEdge && dragMode == DragMode::none
       && havePlayheadFrameIndex && zoomX > 1.0e-12 && width > 1)
    {
        const double playheadX =
            (double) width - (viewRightFrame - (double) playheadFrameIndex) * zoomX;

        constexpr double enterPx = 10.0;  // trigger when within 10px of right edge
        constexpr double exitPx  = 80.0;  // re-arm only after playhead is at least 80px left

        if (playheadX <= (double) width - exitPx)
            autoRefollowArmed = true;

        // Only trigger if playhead is close to the right edge (not far offscreen)
        if (autoRefollowArmed
            && playheadX >= (double) width - enterPx
            && playheadX <= (double) width + enterPx)
        {
            followRightEdge = true;
            autoRefollowArmed = false;
            followButton.setToggleState (true, juce::dontSendNotification);

            // Immediately compute the followed view so it takes effect now
            const double visibleFrames = (double) width / zoomX;
            constexpr double playheadXFrac = 0.5; // keep playhead centered while following
            viewRightFrame = (double) playheadFrameIndex + visibleFrames * (1.0 - playheadXFrac);

            clampViewRightFrame (width);
        }
    }
    
    if (width <= 1 || height <= 0)
        return;

    rebuildStaticBackgroundIfNeeded();
    if (cachedStaticBackground.isValid())
        g.drawImageAt (cachedStaticBackground, 0, 0);
        // [BEGIN ROLLING-LRA-SPLITTER-HOIST-SELECTEDLEVEL]
        // Must live outside the plot-clip scope because we use it later for debug text.
        const int selectedLevel = selectBestLevelForCurrentZoom (width);
        // [END ROLLING-LRA-SPLITTER-HOIST-SELECTEDLEVEL]

        // [BEGIN ROLLING-LRA-SPLITTER-CLIP-PLOT]
        {
            juce::Graphics::ScopedSaveState plotClip (g);
            g.reduceClipRegion (mainPlotArea.toNearestInt());
        // [END ROLLING-LRA-SPLITTER-CLIP-PLOT]

    // [ROLLING-LRA-SPLITTER] selectedLevel moved above (outside clip scope)

    buildVisibleGroupsForLevel (selectedLevel, width, scratchVisibleGroups, scratchVisibleEndFrameIndex); // [TIMEBASE-FIX]

    const size_t n = scratchVisibleGroups.size();
    // [BEGIN ROLLING-LRA-SPLITTER-USE-PLOT-H]
    const float w = bounds.getWidth();
    const float h = mainPlotArea.getHeight(); // plot height changes when rolling lane is resized
    // [END ROLLING-LRA-SPLITTER-USE-PLOT-H]

    if (n >= 2)
    {
        computeRepresentativeCurves (scratchVisibleGroups, scratchRepMomentaryDb, scratchRepShortTermDb);

        scratchPathBandM.clear();
        scratchPathBandS.clear();

        if (showLines)
        {
            scratchPathRepM.clear();
            scratchPathRepS.clear();
        }

        if (showRollingLra)
        scratchPathRollingLra.clear(); // [ROLLING-LRA]

        if (showGate)
            scratchPathGate.clear(); // [LRAG]

        const float bandRangeThresholdDb = 3.0f;

        bool startedRepM = false, startedRepS = false;
        bool startedGate = false; // [LRAG]

        for (size_t i = 0; i < n; ++i)
        {
            const double xD = (double) w - (viewRightFrame - (double) scratchVisibleEndFrameIndex[i]) * zoomX;
            float x = (float) xD;
            if (x < -10.0f)
                continue;

            const auto& grp = scratchVisibleGroups[i];

            const float yMM = dbToY (grp.momentaryMaxDb, h);
            const float yMm = dbToY (grp.momentaryMinDb, h);
            const float ySM = dbToY (grp.shortTermMaxDb, h);
            const float ySm = dbToY (grp.shortTermMinDb, h);

            if (showBands && selectedLevel > 0)
            {
                const float rangeMM = grp.momentaryMaxDb - grp.momentaryMinDb;
                const float rangeSM = grp.shortTermMaxDb - grp.shortTermMinDb;

                if (rangeSM >= bandRangeThresholdDb)
                {
                    scratchPathBandS.startNewSubPath (x, ySM);
                    scratchPathBandS.lineTo (x, ySm);
                }

                if (rangeMM >= bandRangeThresholdDb)
                {
                    scratchPathBandM.startNewSubPath (x, yMM);
                    scratchPathBandM.lineTo (x, yMm);
                }
            }

            if (showLines)
            {
                const float yRepM = dbToY (scratchRepMomentaryDb[i], h);
                const float yRepS = dbToY (scratchRepShortTermDb[i], h);

                if (! startedRepS) { scratchPathRepS.startNewSubPath (x, yRepS); startedRepS = true; }
                else               { scratchPathRepS.lineTo         (x, yRepS); }

                if (! startedRepM) { scratchPathRepM.startNewSubPath (x, yRepM); startedRepM = true; }
                else               { scratchPathRepM.lineTo          (x, yRepM); }
            }

            // [LRAG] Gate curve path
            if (showGate)
            {
                const float gateCenter = 0.5f * (grp.gateMinDb + grp.gateMaxDb);
                const float yGate = dbToY (gateCenter, h);

                if (! startedGate)
                {
                    scratchPathGate.startNewSubPath (x, yGate);
                    startedGate = true;
                }
                else
                {
                    scratchPathGate.lineTo (x, yGate);
                }
            }
        }

       // Bands
       // Draw momentary first, then short-term on top.
       if (showBands && selectedLevel > 0)
       {
           g.setColour (juce::Colour::fromRGB (95, 117, 140).withMultipliedAlpha (0.7f)); // momentary
           g.strokePath (scratchPathBandM, juce::PathStrokeType (1.2f));

           g.setColour (juce::Colour::fromRGB (0, 80, 180).withMultipliedAlpha (0.6f));   // short-term (on top)
           g.strokePath (scratchPathBandS, juce::PathStrokeType (1.0f));
       } 

        // Lines (stroked paths only)
        // Draw momentary first, then short-term on top.
        if (showLines)
        {
            g.setColour (juce::Colour::fromRGB (95, 117, 140)); // momentary
            g.strokePath (scratchPathRepM, juce::PathStrokeType (2.0f));

            g.setColour (juce::Colour::fromRGB (0, 80, 180).withMultipliedAlpha (0.95f)); // short-term (on top)
            g.strokePath (scratchPathRepS, juce::PathStrokeType (1.5f));
        }

        // [LRAG] Gate line (drawn on top)
        if (showGate && startedGate)
        {
            g.setColour (juce::Colours::yellow.withMultipliedAlpha (0.9f));
            g.strokePath (scratchPathGate, juce::PathStrokeType (1.5f));
        }
    }

    // [BEGIN ROLLING-LRA-PLAYHEAD-MOVED]
    // [PLAYHEAD-LINE] moved ниже (after rLRA drawing) so it overlays the rLRA lane too.
    // [END ROLLING-LRA-PLAYHEAD-MOVED]

    // [BEGIN MTDM-THRESH-UI-PAINT-CALL]
    // Draw MTDM thresholds over the curves (O(1): 4 lines + 4 handles).
    // Done here so rulers/labels remain on top.
    drawMtdmThresholdOverlay (g);
    // [END MTDM-THRESH-UI-PAINT-CALL]
    // [BEGIN ROLLING-LRA-SPLITTER-UNCLIP-AFTER-PLOT]
    // End plot clip so rulers / rolling lane can draw in the bottom reserved area.
    // (plotClip ScopedSaveState will restore when it goes out of scope; we end scope here.)
    // NOTE: this relies on plotClip being created just above the curve drawing.
    // We close its scope by wrapping the plotted section in braces.
    // [END ROLLING-LRA-SPLITTER-UNCLIP-AFTER-PLOT]

    // [BEGIN ROLLING-LRA-SPLITTER-CLIP-PLOT-END]
    } // end plot clip scope
    // [END ROLLING-LRA-SPLITTER-CLIP-PLOT-END]

    //==============================================================================
    // [DBFS-SCALE] Right-side dBFS scale (overlay)
    // Shows the currently visible dB range (affected by zoomY).
    //==============================================================================
    {
        // Keep the scale above the time ruler area
        // [BEGIN ROLLING-LRA-SPLITTER-DBSCALE-TRIM]
        const float reservedBottomPx = bounds.getHeight() - mainPlotArea.getHeight();
        const auto scaleArea = bounds.withTrimmedBottom (reservedBottomPx);
        // [END ROLLING-LRA-SPLITTER-DBSCALE-TRIM]

        const float scaleH = scaleArea.getHeight();
        if (scaleH > 20.0f)
        {
            const float topDb = (float) viewTopDb;
            const float effectiveRange = (float) (baseDbRange / zoomY);
            const float bottomDb = topDb - effectiveRange;

            // Choose a tick step based on available pixel height
            static constexpr float candidates[] = { 1.0f, 2.0f, 3.0f, 6.0f, 10.0f, 12.0f, 20.0f, 30.0f };
            const float desiredTickSpacingPx = 32.0f;
            const float approxTicks = juce::jmax (1.0f, scaleH / desiredTickSpacingPx);
            const float approxStepDb = effectiveRange / approxTicks;

            float stepDb = candidates[(int) (sizeof (candidates) / sizeof (candidates[0])) - 1];
            for (float s : candidates)
            {
            if (s >= approxStepDb) { stepDb = s; break; }
            }

            const float rightX = scaleArea.getRight() - 1.0f;
            const float tickLen = 6.0f;
            const float labelWidth = 56.0f;
            const float labelHeight = 16.0f;

            g.setFont (14.0f); // [UI-FONTS]

            // Optional header
            g.setColour (juce::Colours::white.withMultipliedAlpha (0.6f));
            g.drawText ("LUFS",
                        (int) (rightX - labelWidth - 2.0f),
                        (int) (scaleArea.getY() + 2.0f),
                        (int) labelWidth,
                        12,
                        juce::Justification::right);

            // Major ticks: from first tick >= bottomDb up to topDb
            const float firstTick = std::ceil (bottomDb / stepDb) * stepDb;

            g.setColour (juce::Colours::white.withMultipliedAlpha (0.55f));

            for (float db = firstTick; db <= topDb + 0.001f; db += stepDb)
            {
                const float yLocal = dbToY (db, scaleH);
                const float y = scaleArea.getY() + yLocal;

                // Tick
                g.drawLine (rightX - tickLen, y, rightX, y, 1.0f);

                // Label (integer dB)
                const int dbInt = (int) std::round (db);
                g.drawText (juce::String (dbInt),
                            (int) (rightX - tickLen - labelWidth - 2.0f),
                            (int) (y - labelHeight * 0.5f),
                            (int) labelWidth,
                            (int) labelHeight,
                            juce::Justification::right);
            }

            // A faint vertical line to separate the scale
            g.setColour (juce::Colours::white.withMultipliedAlpha (0.18f));
            g.drawLine (rightX, scaleArea.getY(), rightX, scaleArea.getBottom(), 1.0f);
        }
    }

    //==============================================================================
    // [ROLLING-LRA] Rolling LRA curve strip (0..20 LU) drawn above time ruler
    //==============================================================================

    if (showRollingLra && haveNowFrameIndex && zoomX > 1.0e-12)
    {
        const auto timeRuler = getTimeRulerArea();
        const auto dbRuler   = getDbRulerArea();

        // [BEGIN ROLLING-LRA-SPLITTER-HEIGHT]
        const int rollingH = rollingLaneHeightPx;
        // [END ROLLING-LRA-SPLITTER-HEIGHT]
        int yTop = timeRuler.getY() - rollingH;

        if (yTop < 0)
            yTop = 0;

        const int graphW = juce::jmax (1, getWidth() - dbRuler.getWidth());
        juce::Rectangle<int> rollingArea (0, yTop, graphW, rollingH);

        // [BEGIN ROLLING-LRA-SPLITTER-DIVIDER-DRAW]
        {
            const int dividerY = rollingArea.getY();
            g.setColour (juce::Colours::white.withMultipliedAlpha (0.18f));
            g.drawHorizontalLine (dividerY, 0.0f, (float) rollingArea.getRight());

            // Slightly brighter "handle" segment on the left
            g.setColour (juce::Colours::white.withMultipliedAlpha (0.30f));
            g.drawLine (8.0f, (float) dividerY, 38.0f, (float) dividerY, 2.0f);
        }
        // [END ROLLING-LRA-SPLITTER-DIVIDER-DRAW]

        if (rollingArea.getWidth() > 20 && rollingArea.getHeight() > 12)
        {
            juce::Graphics::ScopedSaveState ss (g);
            g.reduceClipRegion (rollingArea);

            // Background strip
            g.setColour (juce::Colours::black.withMultipliedAlpha (0.20f));
            g.fillRect (rollingArea);

            // Grid (0..20 LU)
            const float vMin = 0.0f;
            const float vMax = 20.0f;

            auto valueToY = [&] (float v) -> float
            {
                v = juce::jlimit (vMin, vMax, v);
                const float norm = (v - vMin) / (vMax - vMin);
                return (float) rollingArea.getBottom() - norm * (float) rollingArea.getHeight();
            };

            g.setColour (juce::Colours::white.withMultipliedAlpha (0.10f));
            for (float v : { 0.0f, 5.0f, 10.0f, 15.0f, 20.0f })
                g.drawHorizontalLine ((int) std::round (valueToY (v)), 0.0f, (float) rollingArea.getRight());

            // Label
            g.setColour (juce::Colours::white.withMultipliedAlpha (0.55f));
            g.setFont (12.0f);
            g.drawText ("rLRA " + juce::String (rollingWindowSecondsCached) + "s (LU)",
                        rollingArea.getX() + 6,
                        rollingArea.getY() + 2,
                        rollingArea.getWidth() - 12,
                        14,
                        juce::Justification::left);

            // Visible time range in seconds
            const juce::int64 rightFrameI = (juce::int64) std::floor (viewRightFrame);
            const juce::int64 visibleFrames = (juce::int64) std::ceil ((double) getWidth() / zoomX);
            const juce::int64 leftFrameI = rightFrameI - visibleFrames;

            const juce::int64 leftSecondWanted  = floorDivInt64 (leftFrameI, 60) - 2;  // overscan
            const juce::int64 rightSecondWanted = floorDivInt64 (rightFrameI, 60) + 2;

            const juce::int64 latestSecond = floorDivInt64 (nowFrameIndex, 60);
            const juce::int64 earliestSecond = latestSecond - (juce::int64) (secondsCapacity - 1);

            const juce::int64 minSecond = juce::jmax (earliestSecond, leftSecondWanted);
            const juce::int64 maxSecond = juce::jmin (latestSecond, rightSecondWanted);

            const juce::int64 totalSeconds = (maxSecond - minSecond + 1);
            if (totalSeconds >= 2 && secondsCapacity > 0)
            {
                // Cap points for performance (simple chunking)
                const int maxPoints = juce::jlimit (256, 4096, (int) std::round ((double) rollingArea.getWidth() * 1.20));
                const juce::int64 step = (juce::int64) juce::jmax (1, (int) std::ceil ((double) totalSeconds / (double) maxPoints));

                // Negative-safe floor-to-grid
                auto floorDivLocal = [] (juce::int64 a, juce::int64 b) -> juce::int64
                {
                    if (b <= 0) return 0;
                    if (a >= 0) return a / b;
                    return - ((-a + b - 1) / b);
                };

                auto floorToGrid = [&] (juce::int64 v, juce::int64 grid) -> juce::int64
                {
                    if (grid <= 0) return v;
                    return floorDivLocal (v, grid) * grid;
                };

                juce::int64 firstSecond = floorToGrid (minSecond, step);
                if (firstSecond > minSecond)
                    firstSecond -= step;

                firstSecond = juce::jmax (earliestSecond, firstSecond - step);

                scratchPathRollingLra.clear();
                bool started = false;

                for (juce::int64 s0 = firstSecond; s0 <= maxSecond; s0 += step)
                {
                    const juce::int64 s1 = juce::jmin (maxSecond, s0 + step - 1);

                    // Aggregate chunk: take MAX (peak-preserving) rolling LRA in this chunk
                    float vmax = -1.0f;
                    bool any = false;

                    for (juce::int64 s = s0; s <= s1; ++s)
                    {
                        const int slot = wrapSecondSlot (s);
                        if (secAbsIndexTag[(size_t) slot] != s)
                            continue;

                        any = true;
                        vmax = juce::jmax (vmax, secRollingLraLu[(size_t) slot]);
                    }

                    if (! any)
                        continue;

                    const float y = valueToY (juce::jlimit (vMin, vMax, vmax));

                    // Map this point to the same X timeline as the main graph.
                    // Use chunk end time (s1) -> endFrame = (s1+1)*60
                    const juce::int64 endFrame = (s1 + 1) * 60;

                    const double xD = (double) getWidth() - (viewRightFrame - (double) endFrame) * zoomX;
                    float x = (float) xD;

                    if (x < -10.0f)
                        continue;

                    if (x > (float) getWidth() + 10.0f)
                        continue;

                    if (! started)
                    {
                        scratchPathRollingLra.startNewSubPath (x, y);
                        started = true;
                    }
                    else
                    {
                        scratchPathRollingLra.lineTo (x, y);
                    }
                }

                if (started)
                {
                    g.setColour (juce::Colours::limegreen.withMultipliedAlpha (0.90f));
                    g.strokePath (scratchPathRollingLra, juce::PathStrokeType (1.6f));
                }
            }

            // Border
            g.setColour (juce::Colours::white.withMultipliedAlpha (0.12f));
            g.drawRect (rollingArea);
        }
    }

    // [BEGIN ROLLING-LRA-PLAYHEAD-OVERLAY-ALL]
    // Draw playhead on top of plot + rLRA lane (but not into time ruler).
    if (haveNowFrameIndex && havePlayheadFrameIndex && zoomX > 1.0e-12)
    {
        const double framesFromRight = viewRightFrame - (double) playheadFrameIndex;
        float x = (float) ((double) width - framesFromRight * zoomX);

        x = std::floor (x) + 0.5f; // snap to pixel center

        const float yBottom = (float) getTimeRulerArea().getY(); // includes rLRA lane area above ruler
        if (x >= -2.0f && x <= (float) width + 2.0f && yBottom > 1.0f)
        {
            g.setColour (juce::Colours::white.withMultipliedAlpha (0.55f));
            g.drawLine (x, 0.0f, x, yBottom, 1.0f);
        }
    }
    // [END ROLLING-LRA-PLAYHEAD-OVERLAY-ALL]

    //==========================================================================
    // [RULER-FRAMES]  [VIEW-NAV]
    // Ruler must use the SAME x-mapping as the curves:
    // x = width - (viewRightFrame - tickFrame) * zoomX
    //==========================================================================
    const float rulerHeight   = 22.0f;   // [UI-FONTS] bigger
    const float rulerBaseY    = (float) height - 2.0f;
    const float tickTopY      = rulerBaseY - 6.0f;
    const float textTopY      = rulerBaseY - 14.0f;

    if (zoomX > 1.0e-12 && visualFrameRate > 0.0 && haveNowFrameIndex)
    {
        const juce::int64 rightFrameI = (juce::int64) std::floor (viewRightFrame);
        const double safeZoomX = zoomX;

        const double framesByWidth = (double) width / safeZoomX;
        const juce::int64 visibleFrames = (juce::int64) std::ceil (juce::jmax (1.0, framesByWidth));

        // Earliest frame we can possibly still have in our history window
        const juce::int64 earliestStoredFrame = nowFrameIndex - (juce::int64) rawCapacityFrames + 1;

        const juce::int64 leftFrame = juce::jmax (earliestStoredFrame, rightFrameI - visibleFrames);

        const double tickStepSec = getTickStepSecondsWithHysteresis (width);
        const juce::int64 tickStepFrames = (juce::int64) juce::jmax (1.0, std::round (tickStepSec * visualFrameRate));

        // Negative-safe floor-to-multiple
        auto floorDivInt64 = [] (juce::int64 a, juce::int64 b) -> juce::int64
        {
            if (b <= 0) return 0;
            if (a >= 0) return a / b;
            return - ((-a + b - 1) / b);
        };

        const juce::int64 lastTickFrame = floorDivInt64 (rightFrameI, tickStepFrames) * tickStepFrames;

        g.setColour (juce::Colours::white);
        g.setFont (14.0f);                   // [UI-FONTS] bigger

        for (juce::int64 tickFrame = lastTickFrame; tickFrame >= leftFrame; tickFrame -= tickStepFrames)
        {
            const float x = (float) ((double) width - (viewRightFrame - (double) tickFrame) * zoomX);

            if (x < -2.0f)
                break;

            if (x > (float) width + 2.0f)
                continue;

            g.drawLine (x, tickTopY, x, rulerBaseY, 1.0f);

            const double tSec = (double) tickFrame / visualFrameRate
                              + processor.getTimecodeOffsetSeconds()
                              + processor.getUserTimecodeOffsetSeconds();
            const float textWidth = 84.0f;

            g.drawText (formatTimeHMS (tSec),
                        x - textWidth * 0.5f,
                        textTopY,
                        textWidth,
                        rulerHeight,
                        juce::Justification::centred);
        }
    }

    // Overlay info
    g.setColour (juce::Colours::white);
    g.setFont (14.0f);

    const auto& Lsel         = levels[(size_t) selectedLevel];
    const double spanSeconds = (double) Lsel.spanFrames / visualFrameRate;

    const double tickStepNow = getTickStepSecondsWithHysteresis (width);

    juce::String info = "Level: " + juce::String (selectedLevel) +
                        " (span " + juce::String (Lsel.spanFrames) + " frames, " +
                        juce::String (spanSeconds, 3) + " s)" +
                        " | ZoomX: " + juce::String (zoomX, 4) +
                        " | ZoomY: " + juce::String (zoomY, 2) +
                        " | Tick: " + juce::String (tickStepNow, 3) + "s" +
                        " | Bands: " + juce::String (showBands ? "ON" : "OFF") +
                        " | Lines: " + juce::String (showLines ? "ON" : "OFF");
                        // [DEBUG-OVERLAY] Put debug fields on a second line so they don't get clipped
                        juce::String info2;

                        info2 += "Follow:" + juce::String (followRightEdge ? "ON" : "OFF")
                              +  " viewRight:" + juce::String (viewRightFrame, 1);

                        if (haveNowFrameIndex)
                            info2 += " now:" + juce::String (nowFrameIndex);
                        if (havePlayheadFrameIndex)
                            info2 += " play:" + juce::String (playheadFrameIndex);

                        info2 += " TCoff:" + juce::String (processor.getTimecodeOffsetSeconds(), 3) + "s"
                              +  " TCu:"   + juce::String (processor.getUserTimecodeOffsetSeconds(), 3) + "s"
                              +  " hostSamp:" + juce::String (processor.hostHasTimeInSamples() ? "Y" : "N")
                              +  " hostSec:"  + juce::String (processor.hostHasTimeInSeconds() ? "Y" : "N");

                        if (processor.hostHasTimeInSamples())
                            info2 += " S:" + juce::String (processor.getLastHostTimeInSamples());
                        if (processor.hostHasTimeInSeconds())
                            info2 += " s:" + juce::String (processor.getLastHostTimeInSeconds(), 3);

                        // Draw it
                        const int overlayX = 8;
                        const int overlayW = (int) std::min (bounds.getWidth() - 16.0f, 1600.0f);

                        // Main line
                        g.setColour (juce::Colours::white);
                        g.setFont (14.0f);
                        g.drawFittedText (info, overlayX, 6, overlayW, 18, juce::Justification::topLeft, 1);

                        // Debug line
                        g.setColour (juce::Colours::white.withMultipliedAlpha (0.85f));
                        g.setFont (12.0f);
                        g.drawFittedText (info2, overlayX, 24, overlayW, 16, juce::Justification::topLeft, 1);
}
// [END VHC-RENDER-PAINT]

//==============================================================================
// Cached background  [CACHE-STATIC]
//==============================================================================
// [BEGIN VHC-RENDER-CACHED-BG]
void VolumeHistoryComponent::markStaticBackgroundDirty() noexcept
{
    staticBackgroundDirty = true;
}

void VolumeHistoryComponent::rebuildStaticBackgroundIfNeeded()
{
    const int w = getWidth();
    const int h = getHeight();

    if (w <= 0 || h <= 0)
        return;

    const bool sizeChanged  = (w != cachedBgW || h != cachedBgH);
    const bool zoomYChanged = (std::abs (zoomY - cachedBgZoomY) > 1.0e-12);

    const bool topDbChanged = (std::abs (viewTopDb - cachedBgTopDb) > 1.0e-12); // [VIEW-NAV]

    if (! staticBackgroundDirty && ! sizeChanged && ! zoomYChanged && ! topDbChanged
        && cachedStaticBackground.isValid())
        return;

    cachedBgW = w;
    cachedBgH = h;
    cachedBgZoomY = zoomY;
    cachedBgTopDb = viewTopDb; // [VIEW-NAV]
    staticBackgroundDirty = false;

    cachedStaticBackground = juce::Image (juce::Image::RGB, w, h, true);
    juce::Graphics gg (cachedStaticBackground);

    auto backgroundColour = juce::Colour::fromRGB (16, 30, 50);
    gg.fillAll (backgroundColour);

    gg.setColour (juce::Colours::darkgrey.withMultipliedAlpha (0.4f));
    const int numLines = 4;
    const float effectiveRange = (float) (baseDbRange / zoomY);
    for (int i = 0; i <= numLines; ++i)
    {
        const float db = (float) viewTopDb - (effectiveRange / (float) numLines) * (float) i;
        const float y  = dbToY (db, (float) h);
        gg.drawHorizontalLine ((int) std::round (y), 0.0f, (float) w);
    }

    const float rulerBaseY = (float) h - 2.0f;
    gg.setColour (juce::Colours::darkgrey.withMultipliedAlpha (0.7f));
    gg.drawHorizontalLine ((int) std::round (rulerBaseY), 0.0f, (float) w);
}
// [END VHC-RENDER-CACHED-BG]

// [BEGIN MTDM-THRESH-UI-IMPL]
//==============================================================================
// [MTDM-THRESH-UI] Threshold overlay implementation (T0..T3)
//==============================================================================

bool VolumeHistoryComponent::mtdmParamsAvailable() const noexcept
{
    for (auto* a : mtdmThreshAtoms)
        if (a == nullptr)
            return false;

    for (auto* p : mtdmThreshParams)
        if (p == nullptr)
            return false;

    return true;
}

void VolumeHistoryComponent::initMtdmParamPointersIfNeeded()
{
    if (mtdmParamPtrsInitialised)
        return;

    auto& apvts = processor.getAPVTS();
    using namespace levelscope::mtdm::ParamIDs;

    mtdmThreshAtoms[0]  = apvts.getRawParameterValue (t0Lufs);
    mtdmThreshAtoms[1]  = apvts.getRawParameterValue (t1Lufs);
    mtdmThreshAtoms[2]  = apvts.getRawParameterValue (t2Lufs);
    mtdmThreshAtoms[3]  = apvts.getRawParameterValue (t3Lufs);

    mtdmThreshParams[0] = apvts.getParameter (t0Lufs);
    mtdmThreshParams[1] = apvts.getParameter (t1Lufs);
    mtdmThreshParams[2] = apvts.getParameter (t2Lufs);
    mtdmThreshParams[3] = apvts.getParameter (t3Lufs);

    mtdmParamPtrsInitialised = true;
}

float VolumeHistoryComponent::yToLufs (float y, float height) const noexcept
{
    if (height <= 0.0f)
        return 0.0f;

    // Must match dbToY() usage in the existing curve drawing:
    // dbToY (value, bounds.getHeight()).
    const float effectiveRange = (float) (baseDbRange / zoomY);
    const float top            = (float) viewTopDb;
    const float bottom         = top - effectiveRange;

    const float normY = juce::jlimit (0.0f, 1.0f, y / height); // y=0 top, y=height bottom
    return bottom + (1.0f - normY) * effectiveRange;
}

void VolumeHistoryComponent::updateMtdmThresholdOverlayGeometry()
{
    initMtdmParamPointersIfNeeded();
    if (! mtdmParamsAvailable())
        return;

    // [BEGIN ROLLING-LRA-SPLITTER-MTDM-PLOT-HEIGHT]
    const int w = getWidth();
    const int hTotal = getHeight();

    // Ensure mainPlotArea is valid even if paint hasn't run yet.
    const auto timeRuler = getTimeRulerArea();
    int plotBottom = timeRuler.getY();
    if (showRollingLra)
        plotBottom -= rollingLaneHeightPx;

    plotBottom = juce::jlimit (1, juce::jmax (1, hTotal), plotBottom);
    mainPlotArea = juce::Rectangle<float> (0.0f, 0.0f, (float) w, (float) plotBottom);

    const int h = (int) juce::jmax (1.0f, mainPlotArea.getHeight());
    // [END ROLLING-LRA-SPLITTER-MTDM-PLOT-HEIGHT]
    if (w <= 1 || h <= 1)
        return;

    // [BEGIN ROLLING-LRA-SPLITTER-MTDM-OVERLAY-AREA]
    // Exclude only the right dB ruler strip; height matches the main plot area.
    const auto dbRuler = getDbRulerArea();
    mtdmOverlayArea = juce::Rectangle<float> (0.0f,
                                             0.0f,
                                             (float) (w - dbRuler.getWidth()),
                                             (float) h);
    // [END ROLLING-LRA-SPLITTER-MTDM-OVERLAY-AREA]

    const float handleW = 34.0f;
    const float handleH = 18.0f;
    const float handleX = mtdmOverlayArea.getX() + 4.0f;

    for (int i = 0; i < 4; ++i)
    {
        const float v = mtdmThreshAtoms[(size_t) i]->load (std::memory_order_relaxed);

        const float yLocal = dbToY (v, (float) h);             // IMPORTANT: use full height (matches curves)
        const float yAbs   = mtdmOverlayArea.getY() + yLocal;   // area starts at y=0 anyway

        const float hy = juce::jlimit (mtdmOverlayArea.getY(),
                                       mtdmOverlayArea.getBottom() - handleH,
                                       yAbs - handleH * 0.5f);

        auto& th = thresholdHandles[(size_t) i];
        th.index      = i;
        th.drawBounds = juce::Rectangle<float> (handleX, hy, handleW, handleH);
        th.hitBounds  = th.drawBounds.expanded (5.0f, 6.0f);
    }
}

// [BEGIN MTDM-THRESH-UI-PUSH-FIX]
void VolumeHistoryComponent::computeOrderedThresholdsWithPush (int changedIndex,
                                                              float newValueLufs,
                                                              float outVals[4]) const noexcept
{
    if (changedIndex < 0 || changedIndex > 3)
        return;

    // Precondition: pointers valid (initMtdmParamPointersIfNeeded() already called)
    auto clampSnap = [&] (int idx, float v) -> float
    {
        if (auto* p = mtdmThreshParams[(size_t) idx])
        {
            const auto r = p->getNormalisableRange();
            v = juce::jlimit ((float) r.start, (float) r.end, v);
            v = r.snapToLegalValue (v);
        }
        return v;
    };

    // 1) Start from current values, clamped/snapped
    for (int i = 0; i < 4; ++i)
    {
        float v = mtdmThreshAtoms[(size_t) i]->load (std::memory_order_relaxed);
        outVals[i] = clampSnap (i, v);
    }

    // 2) Set the dragged threshold (ANCHOR) and do not move it again in this function
    outVals[changedIndex] = clampSnap (changedIndex, newValueLufs);

    // 3) Push to the RIGHT (higher index) upwards as needed
    for (int i = changedIndex + 1; i < 4; ++i)
    {
        const float minAllowed = outVals[i - 1] + kThreshMinGapLu;

        if (outVals[i] < minAllowed)
            outVals[i] = minAllowed;

        outVals[i] = clampSnap (i, outVals[i]);

        // If clamping at the max range prevents satisfying minAllowed, we accept the limit.
        // (At normal ranges this won't happen; this avoids "pulling back" the anchor.)
        if (outVals[i] < minAllowed)
            outVals[i] = outVals[i]; // keep clamped value
    }

    // 4) Push to the LEFT (lower index) downwards as needed
    for (int i = changedIndex - 1; i >= 0; --i)
    {
        const float maxAllowed = outVals[i + 1] - kThreshMinGapLu;

        if (outVals[i] > maxAllowed)
            outVals[i] = maxAllowed;

        outVals[i] = clampSnap (i, outVals[i]);

        // If clamping at the min range prevents satisfying maxAllowed, we accept the limit.
        if (outVals[i] > maxAllowed)
            outVals[i] = outVals[i]; // keep clamped value
    }
}
// [END MTDM-THRESH-UI-PUSH-FIX]

void VolumeHistoryComponent::applyThresholdValuesDuringDrag (const float newVals[4])
{
    // Lazy-start gestures for any param that gets pushed/changed during this drag.
    for (int i = 0; i < 4; ++i)
    {
        auto* p = mtdmThreshParams[(size_t) i];
        auto* a = mtdmThreshAtoms[(size_t) i];
        if (p == nullptr || a == nullptr)
            continue;

        const float oldV = a->load (std::memory_order_relaxed);
        const float newV = newVals[i];

        if (std::abs (newV - oldV) < 1.0e-6f)
            continue;

        if (! threshGestureActive[(size_t) i])
        {
            p->beginChangeGesture();
            threshGestureActive[(size_t) i] = true;
        }

        const float norm = juce::jlimit (0.0f, 1.0f, p->convertTo0to1 (newV));
        p->setValueNotifyingHost (norm);
    }
}

void VolumeHistoryComponent::endAllThresholdGestures()
{
    for (int i = 0; i < 4; ++i)
    {
        if (! threshGestureActive[(size_t) i])
            continue;

        if (auto* p = mtdmThreshParams[(size_t) i])
            p->endChangeGesture();

        threshGestureActive[(size_t) i] = false;
    }
}

void VolumeHistoryComponent::drawMtdmThresholdOverlay (juce::Graphics& g)
{
    initMtdmParamPointersIfNeeded();
    if (! mtdmParamsAvailable())
        return;

    updateMtdmThresholdOverlayGeometry();

    if (mtdmOverlayArea.getWidth() <= 1.0f || mtdmOverlayArea.getHeight() <= 1.0f)
        return;

    // [BEGIN ROLLING-LRA-SPLITTER-MTDM-DRAW-PLOT-H]
    const int h = (int) juce::jmax (1.0f, mainPlotArea.getHeight());
    // [END ROLLING-LRA-SPLITTER-MTDM-DRAW-PLOT-H]

    const juce::Colour lineColours[4] =
    {
        juce::Colours::orange.withMultipliedAlpha (0.85f),     // T0
        juce::Colours::cyan.withMultipliedAlpha   (0.80f),     // T1
        juce::Colours::magenta.withMultipliedAlpha(0.78f),     // T2
        juce::Colours::limegreen.withMultipliedAlpha (0.78f)   // T3
    };

    juce::Graphics::ScopedSaveState ss (g);
    g.reduceClipRegion (mtdmOverlayArea.toNearestInt());

    for (int i = 0; i < 4; ++i)
    {
        const float v = mtdmThreshAtoms[(size_t) i]->load (std::memory_order_relaxed);

        const float y = dbToY (v, (float) h); // IMPORTANT: full height (matches curves)
        const float yAbs = mtdmOverlayArea.getY() + y;

        // Line
        g.setColour (lineColours[i]);
        g.drawLine (mtdmOverlayArea.getX(),
                    yAbs,
                    mtdmOverlayArea.getRight(),
                    yAbs,
                    1.2f);

        // Handle
        const auto& th = thresholdHandles[(size_t) i];
        auto r = th.drawBounds;

        const bool isActive = (thresholdDragging && activeThresholdIndex == i);

        g.setColour (juce::Colours::black.withMultipliedAlpha (isActive ? 0.80f : 0.65f));
        g.fillRoundedRectangle (r, 3.5f);

        g.setColour (lineColours[i].withMultipliedAlpha (isActive ? 1.0f : 0.85f));
        g.drawRoundedRectangle (r, 3.5f, isActive ? 1.6f : 1.0f);

        g.setColour (juce::Colours::white.withMultipliedAlpha (0.95f));
        g.setFont (12.0f);
        g.drawFittedText ("T" + juce::String (i), r.toNearestInt(), juce::Justification::centred, 1);
    }
}
// [END MTDM-THRESH-UI-IMPL]

//==============================================================================
// Geometry helpers
//==============================================================================

// [BEGIN VHC-RENDER-GEOMETRY]
float VolumeHistoryComponent::dbToY (float db, float height) const noexcept
{
    if (height <= 0.0f)
        return 0.0f;

    const float effectiveRange = (float) (baseDbRange / zoomY);
    const float topDb          = (float) viewTopDb; // [VIEW-NAV]
    const float bottomDb       = topDb - effectiveRange;

    const float clamped = juce::jlimit (bottomDb, topDb, db);
    const float norm = (clamped - bottomDb) / effectiveRange;

    return height * (1.0f - norm);
}
// [END VHC-RENDER-GEOMETRY]

//==============================================================================
// Ruler tickStep hysteresis  [RULER-HYST-FIX]
//==============================================================================

// [BEGIN VHC-RENDER-RULER-HYST]
double VolumeHistoryComponent::getTickStepSecondsWithHysteresis (int widthPixels) noexcept
{
    const double pixelsPerSecond = zoomX * visualFrameRate;
    if (pixelsPerSecond <= 1.0e-12 || widthPixels <= 0)
        return 60.0;

    // Keep spacing in a band; add margins so we don't thrash near boundaries.
    const double minSpacingPx = (double) widthPixels / 20.0; // target <= 20 labels
    const double maxSpacingPx = (double) widthPixels / 10.0; // target >= 10 labels

    const double minSwitchPx = minSpacingPx * 0.95; // go coarser only if clearly too dense
    const double maxSwitchPx = maxSpacingPx * 1.05; // go finer only if clearly too sparse

    static constexpr double tickSteps[] = {
        0.5, 1.0, 2.0, 5.0,
        10.0, 15.0, 30.0,
        60.0, 120.0, 300.0, 600.0, 1200.0, 3600.0
    };
    static constexpr int numSteps = (int) (sizeof (tickSteps) / sizeof (tickSteps[0]));

    // Initialize once; DO NOT reset on every zoom tick, otherwise hysteresis can't work.
    if (tickStepIndex < 0 || tickStepIndex >= numSteps)
    {
        int idx = numSteps - 1;
        for (int i = 0; i < numSteps; ++i)
        {
            const double spacing = tickSteps[i] * pixelsPerSecond;
            if (spacing >= minSpacingPx)
            {
                idx = i;
                break;
            }
        }
        tickStepIndex = idx;
        return tickSteps[tickStepIndex];
    }

    // Adjust with guard (no infinite loops)
    for (int guard = 0; guard < numSteps + 2; ++guard)
    {
        const double spacing = tickSteps[tickStepIndex] * pixelsPerSecond;

        if (spacing < minSwitchPx && tickStepIndex < numSteps - 1)
        {
            ++tickStepIndex; // coarser
            continue;
        }

        if (spacing > maxSwitchPx && tickStepIndex > 0)
        {
            --tickStepIndex; // finer
            continue;
        }

        break;
    }

    tickStepIndex = juce::jlimit (0, numSteps - 1, tickStepIndex);
    return tickSteps[tickStepIndex];
}
// [END VHC-RENDER-RULER-HYST]

//==============================================================================
// Representative curves
//==============================================================================

// [BEGIN VHC-RENDER-REP-CURVES]
void VolumeHistoryComponent::computeRepresentativeCurves (const std::vector<FrameGroup>& groups,
                                                          std::vector<float>& repMomentary,
                                                          std::vector<float>& repShortTerm) const
{
    const size_t n = groups.size();

    repMomentary.clear();
    repShortTerm.clear();

    if (n == 0)
        return;

    repMomentary.resize (n);
    repShortTerm.resize (n);

    const float epsilonTrend = 0.1f;

    float prevCenterM = 0.0f;
    float prevCenterS = 0.0f;
    bool hasPrev = false;

    for (size_t i = 0; i < n; ++i)
    {
        const auto& grp = groups[i];

        const float minM    = grp.momentaryMinDb;
        const float maxM    = grp.momentaryMaxDb;
        const float centerM = 0.5f * (minM + maxM);

        const float minS    = grp.shortTermMinDb;
        const float maxS    = grp.shortTermMaxDb;
        const float centerS = 0.5f * (minS + maxS);

        float rawRepM = centerM;
        float rawRepS = centerS;

        if (! hasPrev)
        {
            hasPrev     = true;
            prevCenterM = centerM;
            prevCenterS = centerS;
        }
        else
        {
            const float trendM = centerM - prevCenterM;
            const float trendS = centerS - prevCenterS;

            float alphaM = 0.5f;
            if (trendM >  epsilonTrend)      alphaM = 1.0f;
            else if (trendM < -epsilonTrend) alphaM = 0.0f;

            float alphaS = 0.5f;
            if (trendS >  epsilonTrend)      alphaS = 1.0f;
            else if (trendS < -epsilonTrend) alphaS = 0.0f;

            rawRepM = minM + alphaM * (maxM - minM);
            rawRepS = minS + alphaS * (maxS - minS);

            rawRepM = juce::jlimit (minM, maxM, rawRepM);
            rawRepS = juce::jlimit (minS, maxS, rawRepS);

            prevCenterM = centerM;
            prevCenterS = centerS;
        }

        repMomentary[i] = rawRepM;
        repShortTerm[i] = rawRepS;
    }
}
// [END VHC-RENDER-REP-CURVES]

// [TIMECODE-USER] parse either seconds (e.g. -2.0) or HH:MM:SS (e.g. -00:00:02)

// [BEGIN VHC-RENDER-PARSE-USER-OFFSET]
static bool parseUserOffsetSeconds (juce::String text, double& outSeconds)
{
    text = text.trim();

    if (text.isEmpty())
        return false;

    // If it contains ':' parse as HH:MM:SS
    if (text.containsChar (':'))
    {
        bool neg = false;
        if (text.startsWithChar ('-'))
        {
            neg = true;
            text = text.substring (1).trim();
        }

        auto parts = juce::StringArray::fromTokens (text, ":", "");
        parts.removeEmptyStrings();

        if (parts.size() != 3)
            return false;

        const int hh = parts[0].getIntValue();
        const int mm = parts[1].getIntValue();
        const int ss = parts[2].getIntValue();

        double sec = (double) hh * 3600.0 + (double) mm * 60.0 + (double) ss;
        if (neg) sec = -sec;

        outSeconds = sec;
        return true;
    }

    // Otherwise treat as seconds
    outSeconds = text.getDoubleValue();
    return true;
}
// [END VHC-RENDER-PARSE-USER-OFFSET]