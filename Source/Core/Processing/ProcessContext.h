#pragma once

#include <juce_audio_basics/juce_audio_basics.h>
#include <cstdint>

namespace levelscope
{
    // ProcessContext is created by the wrapper/core per audio block and passed through the module chain.
    // RT-safety: must remain a trivial, stack-friendly struct; no owning allocations.
    struct ProcessContext
    {
        // Audio
        juce::AudioBuffer<float>& audio;
        juce::MidiBuffer* midi = nullptr; // may be null if caller does not use MIDI

        // Format
        double sampleRate = 48000.0;
        int numSamples = 0;
        juce::AudioChannelSet channelSet; // channel roles/layout for role-aware processing

        // Transport / timebase (populated by wrapper)
        bool isRealtime = true;
        bool isPlaying = false;
        bool isRecording = false;

        // Discontinuity = seek/loop/jump detected by wrapper timebase logic.
        // Modules may use this to reset detectors, envelopes, lookahead buffers, etc.
        bool isDiscontinuity = false;

        // Analysis guard: when true, analysis should not update (e.g. stop-time freeze behavior).
        bool freezeAnalysis = false;

        // Absolute timeline indices at start of this block (provided by wrapper if known)
        int64_t absoluteSampleIndex = 0;

        // 60 Hz "timeline truth" frame index at block start (optional; wrapper decides availability)
        bool hasFrameIndex60Hz = false;
        int64_t absoluteFrameIndex60Hz = 0;

        ProcessContext (juce::AudioBuffer<float>& bufferIn,
                        juce::MidiBuffer* midiIn,
                        double sampleRateIn,
                        int numSamplesIn,
                        const juce::AudioChannelSet& channelSetIn) noexcept
            : audio (bufferIn),
              midi (midiIn),
              sampleRate (sampleRateIn),
              numSamples (numSamplesIn),
              channelSet (channelSetIn)
        {
        }
    };
} // namespace levelscope
