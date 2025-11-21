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
        Source/VolumeHistoryComponent.cpp
        Source/VolumeHistoryComponent.h
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