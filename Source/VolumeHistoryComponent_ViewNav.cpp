#include "VolumeHistoryComponent.h"
#include "PluginProcessor.h"

//==============================================================================
// VolumeHistoryComponent_ViewNav.cpp
// - view state (follow, viewRightFrame, viewTopDb)
// - wheel gestures (zoom/pan)
// - ruler hit test
// - mouse drag/double click interactions
//==============================================================================
//==============================================================================
// [EDIT-BLOCKS]
//   - [BEGIN VHC-VNAV-RULER-AREAS-RESETS] ... [END VHC-VNAV-RULER-AREAS-RESETS]
//   - [BEGIN VHC-VNAV-CLAMP-PAN-ZOOM]     ... [END VHC-VNAV-CLAMP-PAN-ZOOM]
//   - [BEGIN VHC-VNAV-RESIZED]            ... [END VHC-VNAV-RESIZED]
//   - [BEGIN VHC-VNAV-MOUSE]              ... [END VHC-VNAV-MOUSE]
//==============================================================================

//==============================================================================
// [UI-RULERS] Areas + resets
//==============================================================================

// [BEGIN VHC-VNAV-RULER-AREAS-RESETS]
juce::Rectangle<int> VolumeHistoryComponent::getTimeRulerArea() const
{
    // Keep consistent with paint() ruler height
    const int rulerH = 22; // [UI-FONTS] 150% bigger than before
    return { 0, getHeight() - rulerH, getWidth(), rulerH };
}

juce::Rectangle<int> VolumeHistoryComponent::getDbRulerArea() const
{
    // Right-side strip used for dB scale interaction
    const int rulerW = 72; // room for bigger labels
    const int timeRulerH = getTimeRulerArea().getHeight();

    return { getWidth() - rulerW, 0, rulerW, getHeight() - timeRulerH };
}

void VolumeHistoryComponent::resetXViewDefault()
{
    if (getWidth() <= 1 || visualFrameRate <= 0.0)
        return;

    const double desiredVisibleSeconds = 10.0;
    zoomX = (double) getWidth() / (desiredVisibleSeconds * visualFrameRate);
    zoomX = juce::jlimit (minZoomX, maxZoomX, zoomX);

    followRightEdge = true;
    if (haveNowFrameIndex)
        viewRightFrame = (double) nowFrameIndex;

    clampViewRightFrame (getWidth());
}

void VolumeHistoryComponent::fitXViewMaxZoomOut()
{
    if (getWidth() <= 1 || ! haveNowFrameIndex)
        return;

    // [FIT-FIRST-WRITTEN]
    // Find earliest written L0 frame index currently present in the GUI ring.
    // This avoids fitting the full 3h window when only a few minutes exist.
    const auto& tags = levels[0].absGroupIndexTag;

    juce::int64 minWritten = std::numeric_limits<juce::int64>::max();
    bool haveMin = false;

    for (const auto& t : tags)
    {
        if (t != (juce::int64) -1)
        {
            minWritten = juce::jmin (minWritten, t);
            haveMin = true;
        }
    }

    const juce::int64 right = nowFrameIndex;
    const juce::int64 left  = (haveMin ? minWritten : (right - (juce::int64) rawCapacityFrames + 1));

    const double visibleFrames = (double) (right - left + 1);
    zoomX = (double) getWidth() / juce::jmax (1.0, visibleFrames);
    zoomX = juce::jlimit (minZoomX, maxZoomX, zoomX);

    // Fit view to [left..right] and keep it there (do not force follow)
    followRightEdge = false;
    viewRightFrame = (double) right;
    clampViewRightFrame (getWidth());
}

void VolumeHistoryComponent::resetYViewDefault()
{
    zoomY = 1.0;
    viewTopDb = (double) maxDb; // default top at 0 dBFS
    markStaticBackgroundDirty();
}
// [END VHC-VNAV-RULER-AREAS-RESETS]

//==============================================================================
// [VIEW-NAV] clamp + pan/zoom
//==============================================================================

// [BEGIN VHC-VNAV-CLAMP-PAN-ZOOM]
void VolumeHistoryComponent::clampViewRightFrame (int widthPixels) noexcept
{
    if (! haveNowFrameIndex || zoomX <= 1.0e-12 || widthPixels <= 1)
        return;

    const double visibleFrames = (double) widthPixels / zoomX;

    const double maxRight = (double) nowFrameIndex;
    const double earliest = (double) nowFrameIndex - (double) rawCapacityFrames + 1.0;

    double minRight = earliest + visibleFrames;
    if (minRight > maxRight)
        minRight = maxRight;

    viewRightFrame = juce::jlimit (minRight, maxRight, viewRightFrame);
}

void VolumeHistoryComponent::panTime (float wheelDelta)
{
    const int w = getWidth();
    if (w <= 1 || wheelDelta == 0.0f || zoomX <= 1.0e-12)
        return;

    if (followRightEdge)
        viewRightFrame = (double) nowFrameIndex;

    const double visibleFrames = (double) w / zoomX;

    // Pan by ~15% of the visible range per wheel "unit"
    const double panFrames = (double) wheelDelta * visibleFrames * 0.15;

    // Wheel up (positive) moves to earlier time (left)
    viewRightFrame -= panFrames;

    followRightEdge = false;
    autoRefollowArmed = false; // [AUTO-FOLLOW-HYST] disarm until playhead moves away from edge
    clampViewRightFrame (w);
}

void VolumeHistoryComponent::panDb (float wheelDelta)
{
    const int h = getHeight();
    if (h <= 1 || wheelDelta == 0.0f)
        return;

    const double effectiveRange = (double) baseDbRange / zoomY;

    // Pan by ~10% of visible range per wheel "unit"
    const double panDbAmount = (double) wheelDelta * effectiveRange * 0.10;

    viewTopDb += panDbAmount;

    const double topMin = viewMinDbLimit + effectiveRange;
    const double topMax = viewMaxDbLimit;
    viewTopDb = juce::jlimit (topMin, topMax, viewTopDb);

    markStaticBackgroundDirty();
}

void VolumeHistoryComponent::applyHorizontalZoom (float wheelDelta, float anchorX)
{
    const int w = getWidth();
    if (w <= 1 || wheelDelta == 0.0f)
        return;

    if (! haveNowFrameIndex)
        return;

    const double oldZoom = zoomX;

    //--------------------------------------------------------------------------
    // [FOLLOW-ZOOM-FIX]
    // If Follow is ON, we must base the zoom anchor on the *currently followed*
    // viewRightFrame (playhead-centered), NOT on nowFrameIndex (furthest written).
    // Otherwise, when overwriting earlier material, zoom snaps to the far-right
    // prior written region.
    //--------------------------------------------------------------------------
    if (followRightEdge)
    {
        if (havePlayheadFrameIndex && oldZoom > 1.0e-12)
        {
            const double visibleFrames = (double) w / oldZoom;
            constexpr double playheadXFrac = 0.5; // keep playhead centered while following
            viewRightFrame = (double) playheadFrameIndex + visibleFrames * (1.0 - playheadXFrac);
        }
        else if (haveNowFrameIndex)
        {
            // Fallback only if playhead is unknown
            viewRightFrame = (double) nowFrameIndex;
        }
    }

    const double zoomBase   = 1.1;
    const double zoomFactor = std::pow (zoomBase, (double) wheelDelta);

    zoomX *= zoomFactor;
    zoomX  = juce::jlimit (minZoomX, maxZoomX, zoomX);

    // Keep the timeline frame under the mouse fixed
    const double ax = juce::jlimit (0.0, (double) w, (double) anchorX);

    const double safeOld = (oldZoom > 1.0e-12 ? oldZoom : 1.0e-12);
    const double safeNew = (zoomX   > 1.0e-12 ? zoomX   : 1.0e-12);

    const double frameAtCursor = viewRightFrame - ((double) w - ax) / safeOld;
    viewRightFrame = frameAtCursor + ((double) w - ax) / safeNew;

    // Manual zoom disables follow (consistent with existing pan behavior)
    followRightEdge = false;
    autoRefollowArmed = false; // [AUTO-FOLLOW-HYST]
    hasCustomZoomX = true;

    clampViewRightFrame (w);
}

void VolumeHistoryComponent::applyVerticalZoom (float wheelDelta, float anchorY)
{
    const int h = getHeight();
    if (h <= 1 || wheelDelta == 0.0f)
        return;

    const double oldZoomY = zoomY;

    const double zoomBase   = 1.1;
    const double zoomFactor = std::pow (zoomBase, (double) wheelDelta);

    zoomY *= zoomFactor;
    zoomY  = juce::jlimit (minZoomY, maxZoomY, zoomY);

    // [VIEW-NAV-Y-LIMITS]
    // Prevent impossible states where the visible range exceeds the allowed view limits.
    // effectiveRange = baseDbRange / zoomY must be <= (viewMaxDbLimit - viewMinDbLimit)
    const double maxVisibleRange = viewMaxDbLimit - viewMinDbLimit;
    const double minZoomYByLimits = (double) baseDbRange / maxVisibleRange;
    zoomY = juce::jmax (zoomY, minZoomYByLimits);

    // Keep the dB under the mouse fixed while zooming
    const double ay = juce::jlimit (0.0, (double) h, (double) anchorY);
    const double norm = 1.0 - (ay / (double) h);

    const double effectiveOld = (double) baseDbRange / oldZoomY;
    const double bottomOld = viewTopDb - effectiveOld;
    const double dbAtCursor = bottomOld + norm * effectiveOld;

    const double effectiveNew = (double) baseDbRange / zoomY;
    const double bottomNew = dbAtCursor - norm * effectiveNew;
    double topNew = bottomNew + effectiveNew;

    const double topMin = viewMinDbLimit + effectiveNew; // ensures bottom >= viewMinDbLimit
    const double topMax = viewMaxDbLimit;                // allows positive dBFS
    topNew = juce::jlimit (topMin, topMax, topNew);
    
    viewTopDb = topNew;

    markStaticBackgroundDirty();
}
// [END VHC-VNAV-CLAMP-PAN-ZOOM]

//==============================================================================
// Component lifecycle
//==============================================================================

// [BEGIN VHC-VNAV-RESIZED]
void VolumeHistoryComponent::resized()
{
    if (! hasCustomZoomX)
    {
        const int w = getWidth();
        if (w > 0 && visualFrameRate > 0.0)
        {
            const double desiredVisibleSeconds = 10.0;
            zoomX = (double) w / (desiredVisibleSeconds * visualFrameRate);
            zoomX = juce::jlimit (minZoomX, maxZoomX, zoomX);
        }
    }

    // Reserve scratch buffers
    {
        const int w = juce::jmax (1, getWidth());
        const size_t reserveCount = (size_t) juce::jlimit (256, 8192, w + 512);

        scratchVisibleGroups.reserve (reserveCount);
        scratchVisibleEndFrameIndex.reserve (reserveCount); // [TIMEBASE-FIX]
        scratchRepMomentaryDb.reserve (reserveCount);
        scratchRepShortTermDb.reserve (reserveCount);
    }

    // Do NOT reset tickStepIndex here; hysteresis should smoothly adapt.
    markStaticBackgroundDirty();

    // [FOLLOW-BUTTON] small toggle in the top-right
    followButton.setBounds (getWidth() - 88, 6, 80, 22);
    gateButton.setBounds   (getWidth() - 172, 6, 76, 22);
    rollingLraButton.setBounds (getWidth() - 256, 6, 76, 22);
}
// [END VHC-VNAV-RESIZED]

//==============================================================================
// Mouse
//==============================================================================

// [BEGIN VHC-VNAV-MOUSE]
void VolumeHistoryComponent::mouseWheelMove (const juce::MouseEvent& event,
                                             const juce::MouseWheelDetails& wheel)
{
    // Use deltaY when available, otherwise fall back to deltaX
    const float wheelDelta = (wheel.deltaY != 0.0f ? wheel.deltaY : wheel.deltaX);
    if (wheelDelta == 0.0f)
        return;

    // Gesture map:
    // - Shift + wheel: pan time
    // - Alt   + wheel: zoom Y at mouse
    // - Ctrl  + wheel: pan Y
    // - wheel: zoom X at mouse
    if (event.mods.isShiftDown())
        panTime (wheelDelta);
    else if (event.mods.isAltDown())
        applyVerticalZoom (wheelDelta, event.position.y);
    else if (
       #if JUCE_MAC
            event.mods.isCommandDown()
   #else
        event.mods.isCtrlDown()
   #endif
    )
        panDb (wheelDelta);
    else
        applyHorizontalZoom (wheelDelta, event.position.x);

    repaint();
}

void VolumeHistoryComponent::mouseDoubleClick (const juce::MouseEvent& event)
{
    const auto p = event.getPosition();

    // Double-click dB ruler resets Y view
    if (getDbRulerArea().contains (p))
    {
        resetYViewDefault();
        repaint();
        return;
    }

    // Double-click time ruler resets X view
    // Shift + double-click: fit/max zoom out
    if (getTimeRulerArea().contains (p))
    {
        if (event.mods.isShiftDown())
            fitXViewMaxZoomOut();
        else
            resetXViewDefault();

        followButton.setToggleState (followRightEdge, juce::dontSendNotification);
        repaint();
        return;
    }
}

void VolumeHistoryComponent::mouseDown (const juce::MouseEvent& event)
{
    const auto p = event.getPosition();

    // [TIMECODE-USER] Right-click on time ruler to set/reset user timecode offset
    if (getTimeRulerArea().contains (p) && event.mods.isPopupMenu())
    {
        juce::PopupMenu m;
        m.addItem (1, "Timecode offset: -2.0 s (Cubase common preroll)");
        m.addItem (2, "Timecode offset: 0.0 s (reset)");
        m.addSeparator();
        m.addItem (3, "Offset -1.0 s");
        m.addItem (4, "Offset +1.0 s");
        m.addItem (5, "Offset -0.1 s");
        m.addItem (6, "Offset +0.1 s");

        m.showMenuAsync (juce::PopupMenu::Options(),
                         [this] (int r)
                         {
                             auto setOff = [this] (double s)
                             {
                                 processor.setUserTimecodeOffsetSeconds (s);
                             };

                             auto addOff = [this] (double ds)
                             {
                                 processor.setUserTimecodeOffsetSeconds (processor.getUserTimecodeOffsetSeconds() + ds);
                             };

                             if (r == 1) { setOff (-2.0); repaint(); return; }
                             if (r == 2) { setOff ( 0.0); repaint(); return; }
                                     if (r == 3) { addOff (-1.0); repaint(); return; }
                             if (r == 4) { addOff (+1.0); repaint(); return; }
                             if (r == 5) { addOff (-0.1); repaint(); return; }
                             if (r == 6) { addOff (+0.1); repaint(); return; }
                         });
        return;
    }

    // [UI-RULERS] Drag time ruler to pan time
    if (getTimeRulerArea().contains (p))
    {
        dragMode = DragMode::timeRuler;
        dragStartPos = p;
        dragStartViewRightFrame = viewRightFrame;
        followRightEdge = false;
        autoRefollowArmed = false; // [AUTO-FOLLOW-HYST] disarm until playhead moves away from edge
        followButton.setToggleState (false, juce::dontSendNotification);
        return;
    }

    // [UI-RULERS] Drag dB ruler to pan dB
    if (getDbRulerArea().contains (p))
    {
        dragMode = DragMode::dbRuler;
        dragStartPos = p;
        dragStartViewTopDb = viewTopDb;
        return;
    }

    // Existing behavior: Shift-click toggles bands; plain click toggles lines.
    if (event.mods.isShiftDown())
        showBands = ! showBands;
    else
        showLines = ! showLines;

    repaint();
}

void VolumeHistoryComponent::mouseDrag (const juce::MouseEvent& event)
{
    if (dragMode == DragMode::none)
        return;

    const auto p = event.getPosition();

    if (dragMode == DragMode::timeRuler)
    {
        const int dx = p.x - dragStartPos.x;

        // "Grab" behavior: drag right -> earlier time (content moves right)
        viewRightFrame = dragStartViewRightFrame - (double) dx / (zoomX > 1.0e-12 ? zoomX : 1.0e-12);
        clampViewRightFrame (getWidth());
        repaint();
        return;
    }

    if (dragMode == DragMode::dbRuler)
    {
        const int dy = p.y - dragStartPos.y;

        const double effectiveRange = (double) baseDbRange / zoomY;
        const double dbPerPixel = effectiveRange / (double) juce::jmax (1, getHeight());

        // [DB-DRAG-DIRECTION] drag down -> curves go down (topDb increases)
        viewTopDb = dragStartViewTopDb + (double) dy * dbPerPixel;

        const double topMin = viewMinDbLimit + effectiveRange;
        const double topMax = viewMaxDbLimit;
        viewTopDb = juce::jlimit (topMin, topMax, viewTopDb);

        markStaticBackgroundDirty();
        repaint();
        return;
    }
}

void VolumeHistoryComponent::mouseUp (const juce::MouseEvent&)
{
    dragMode = DragMode::none;
}
// [END VHC-VNAV-MOUSE]