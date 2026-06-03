cmake_minimum_required(VERSION 3.15)

project(LevelScope VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default to Release if not specified
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

# macOS universal binary + deployment target
if (APPLE)
    if (NOT CMAKE_OSX_DEPLOYMENT_TARGET)
        set(CMAKE_OSX_DEPLOYMENT_TARGET "10.13" CACHE STRING "Minimum macOS deployment version" FORCE)
    endif()

    if (NOT CMAKE_OSX_ARCHITECTURES)
        set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64" CACHE STRING "Build architectures for macOS" FORCE)
    endif()
endif()

include(FetchContent)

# Fetch JUCE (adjust tag if you want a different JUCE version)
FetchContent_Declare(
    juce
    GIT_REPOSITORY https://github.com/juce-framework/JUCE.git
    GIT_TAG 7.0.12
)

FetchContent_MakeAvailable(juce)

#==============================
# Core static library (no GUI)
#==============================
add_library(LevelScopeCore STATIC
    Source/Core/LevelScopeHistoryModel.cpp
    Source/Core/LevelScopeHistoryModel.h
    Source/Core/BS1770KWeighting.h
    Source/Core/RunningLoudnessStats.cpp
    Source/Core/RunningLoudnessStats.h

    # Core context-pack umbrella (public-facing)
    Source/Core/PublicAPI.h

    # [BEGIN MTDM-CMAKE-CORESOURCES]
    # Processing / module-chain scaffolding
    Source/Core/Processing/ProcessContext.h
    Source/Core/Processing/IAudioModule.h
    Source/Core/Processing/ProcessorCore.h

    # Processing / modules
    # [BEGIN LVLR-CMAKE-CORESOURCES]
    Source/Core/Processing/Modules/LevelerModule.h
    Source/Core/Processing/Modules/LevelerModule.cpp
    # [END LVLR-CMAKE-CORESOURCES]

    Source/Core/Processing/Modules/MultiThresholdDynamicsModule.h
    Source/Core/Processing/Modules/MultiThresholdDynamicsModule.cpp
    # [END MTDM-CMAKE-CORESOURCES]

    # [BEGIN MTDM-CMAKE-SUC-DSP]
    # Processing / DSP building blocks (Stage D1a)
    Source/Core/Processing/DSP/SpectralUpwardCompressor.h
    Source/Core/Processing/DSP/SpectralUpwardCompressor.cpp
    # [END MTDM-CMAKE-SUC-DSP]

    # [BEGIN LVLR-CMAKE-PARAMIDS]
    Source/Core/Processing/Modules/LevelerParamIDs.h
    # [END LVLR-CMAKE-PARAMIDS]

    # [BEGIN MTDM-CMAKE-PARAMIDS]
    Source/Core/Processing/Modules/MultiThresholdDynamicsParamIDs.h
    # [END MTDM-CMAKE-PARAMIDS]

    # [BEGIN MTDM-CMAKE-PARAMGROUPS]
    Source/Core/Processing/Modules/MultiThresholdDynamicsParamGroups.h
    # [END MTDM-CMAKE-PARAMGROUPS]

    # [BEGIN MTDM-CMAKE-BUC-DSP]
    Source/Core/Processing/DSP/BroadbandUpwardCompressor.h
    Source/Core/Processing/DSP/BroadbandUpwardCompressor.cpp
    # [END MTDM-CMAKE-BUC-DSP]

    # [BEGIN MTDM-CMAKE-BDC-DSP]
    Source/Core/Processing/DSP/BroadbandDownwardCompressor.h
    Source/Core/Processing/DSP/BroadbandDownwardCompressor.cpp
    # [END MTDM-CMAKE-BDC-DSP]

    # [BEGIN MTDM-CMAKE-LIMITER-DSP]
    Source/Core/Processing/DSP/LookaheadLimiter.h
    Source/Core/Processing/DSP/LookaheadLimiter.cpp
    # [END MTDM-CMAKE-LIMITER-DSP]
)

target_include_directories(LevelScopeCore
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/Source
)

# [BEGIN MTDM-CMAKE-JUCE-DSP]
target_link_libraries(LevelScopeCore
    PUBLIC
        juce::juce_core
        juce::juce_audio_basics       # ProcessContext uses AudioBuffer/AudioChannelSet
        juce::juce_data_structures    # IAudioModule uses ValueTree
        juce::juce_dsp                # Stage D1a spectral processor (FFT/windowing)
)
# [END MTDM-CMAKE-JUCE-DSP]

# The plugin target
juce_add_plugin(LevelScope
    PRODUCT_NAME "LevelScope"
    COMPANY_NAME "Studio Paradiso"
    BUNDLE_ID "com.StudioParadiso.levelscope"

    IS_SYNTH FALSE
    NEEDS_MIDI_INPUT FALSE
    NEEDS_MIDI_OUTPUT FALSE
    IS_MIDI_EFFECT FALSE
    EDITOR_WANTS_KEYBOARD_FOCUS FALSE
    COPY_PLUGIN_AFTER_BUILD TRUE

    PLUGIN_MANUFACTURER_CODE StPa   # 4 chars
    PLUGIN_CODE LvSc                # 4 chars

    FORMATS VST3
)

juce_generate_juce_header(LevelScope)

target_sources(LevelScope
    PRIVATE
        Source/PluginProcessor.cpp
        Source/PluginProcessor.h
        Source/PluginEditor.cpp
        Source/PluginEditor.h
        Source/DynamicsCurveComponent.cpp
        Source/DynamicsCurveComponent.h
        Source/LevelerCurveComponent.cpp
        Source/LevelerCurveComponent.h
        Source/VolumeHistoryComponent.cpp
        Source/VolumeHistoryComponent.h
        Source/VolumeHistoryComponent_History.cpp
        Source/VolumeHistoryComponent_ViewNav.cpp
        Source/VolumeHistoryComponent_Render.cpp
        Source/LoudnessStatsComponent.cpp
        Source/LoudnessStatsComponent.h
)

target_compile_definitions(LevelScope
    PRIVATE
        JUCE_WEB_BROWSER=0
        JUCE_USE_CURL=0
        JUCE_VST3_CAN_REPLACE_VST2=0
        JUCE_REPORT_APP_USAGE=0
)

target_link_libraries(LevelScope
    PRIVATE
        LevelScopeCore
        juce::juce_audio_utils
        juce::juce_audio_processors
        juce::juce_audio_formats
        juce::juce_audio_devices
        juce::juce_audio_basics
        juce::juce_core
        juce::juce_data_structures
        juce::juce_events
        juce::juce_graphics
        juce::juce_gui_basics
        juce::juce_gui_extra
)

if (MSVC)
    target_compile_options(LevelScope PRIVATE /W4)
else()
    target_compile_options(LevelScope PRIVATE -Wall -Wextra -Wpedantic)
endif()

#==============================
# Standalone app scaffold
#==============================

juce_add_gui_app(LevelScopeStandalone
    PRODUCT_NAME "LevelScope Standalone"
    COMPANY_NAME "Studio Paradiso"
    BUNDLE_ID "com.studioparadiso.levelscope.standalone"
)

target_sources(LevelScopeStandalone
    PRIVATE
        Source/Standalone/StandaloneApp.cpp
        Source/Standalone/StandaloneMainComponent.cpp
        Source/Standalone/StandaloneMainComponent.h
        Source/Standalone/StandaloneSessionModel.h
)

target_compile_definitions(LevelScopeStandalone
    PRIVATE
        JUCE_WEB_BROWSER=0
        JUCE_USE_CURL=0
        JUCE_REPORT_APP_USAGE=0
        JUCE_APPLICATION_NAME_STRING="LevelScope Standalone"
        JUCE_APPLICATION_VERSION_STRING="0.1.0"
)

target_link_libraries(LevelScopeStandalone
    PRIVATE
        LevelScopeCore
        juce::juce_gui_basics
)

if (MSVC)
    target_compile_options(LevelScopeStandalone PRIVATE /W4)
else()
    target_compile_options(LevelScopeStandalone PRIVATE -Wall -Wextra -Wpedantic)
endif()