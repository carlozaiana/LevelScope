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
    GIT_SHALLOW TRUE
)

# JUCE 7 bootstraps juceaide in a nested CMake process. Let that helper inherit
# the configured minimum too, rather than the host SDK's 15.x default (which
# makes JUCE 7's CGWindowListCreateImage call unavailable at compile time).
if (APPLE)
    set(_levelscope_previous_deployment_env "$ENV{MACOSX_DEPLOYMENT_TARGET}")
    set(ENV{MACOSX_DEPLOYMENT_TARGET} "${CMAKE_OSX_DEPLOYMENT_TARGET}")
endif()
FetchContent_MakeAvailable(juce)
if (APPLE)
    set(ENV{MACOSX_DEPLOYMENT_TARGET} "${_levelscope_previous_deployment_env}")
endif()

#==============================
# Core static library (no GUI)
#==============================
add_library(LevelScopeMeasurementCore STATIC
    Source/Core/Measurement/OfflineMeasurement.cpp)
target_include_directories(LevelScopeMeasurementCore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Source)
if (NOT MSVC)
    target_compile_options(LevelScopeMeasurementCore PRIVATE -Wall -Wextra -Wpedantic -ffp-contract=off)
endif()
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
        LevelScopeMeasurementCore
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

add_library(LevelScopeStandaloneSupport STATIC
    Source/Standalone/SessionDocument.cpp
    Source/Standalone/VerifiedSource.cpp
    Source/Core/Proposal/SectionAnalysis.cpp
    Source/Standalone/StandaloneSections.cpp
    Source/Core/Proposal/AnalysisEngine.cpp
    Source/Core/Proposal/AttributionEngine.cpp
    Source/Standalone/StandaloneDiagnostics.cpp
    Source/Core/Proposal/Methods/RuleBasedMethod.cpp
    Source/Standalone/StandaloneProposal.cpp
    Source/Core/Proposal/Model/TargetProfile.cpp
    Source/Core/Proposal/VerificationEngine.cpp
    Source/Standalone/StandaloneCompliance.cpp
    Source/Standalone/StandaloneMedia.cpp
    Source/Standalone/StandaloneAudition.cpp
    Source/Standalone/StandaloneMeasurementService.cpp
    Source/Standalone/StandaloneVideo.cpp
)
if (APPLE)
    enable_language(OBJCXX)
    target_sources(LevelScopeStandaloneSupport PRIVATE Source/Standalone/StandaloneVideo_mac.mm)
    set_source_files_properties(Source/Standalone/StandaloneVideo_mac.mm PROPERTIES COMPILE_FLAGS "-fobjc-arc -Wno-deprecated-declarations")
    target_link_libraries(LevelScopeStandaloneSupport PUBLIC "-framework AVFoundation" "-framework CoreMedia" "-framework CoreVideo")
endif()
target_include_directories(LevelScopeStandaloneSupport PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Source)
file(GLOB_RECURSE _measurement_identity_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/Source/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/Source/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/Source/*.mm")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${_measurement_identity_sources})
set(_measurement_identity_text "${CMAKE_CXX_COMPILER_ID};${CMAKE_CXX_COMPILER_VERSION};${CMAKE_BUILD_TYPE};${CMAKE_CXX_FLAGS};${CMAKE_CXX_FLAGS_RELEASE};no-fast-math;fp-contract-off")
foreach(_measurement_source IN LISTS _measurement_identity_sources)
    file(SHA256 "${_measurement_source}" _measurement_source_hash)
    string(APPEND _measurement_identity_text "${_measurement_source_hash};")
endforeach()
string(SHA256 _measurement_build_identity "${_measurement_identity_text}")
target_compile_definitions(LevelScopeStandaloneSupport PUBLIC LEVELSCOPE_MEASUREMENT_BUILD="${_measurement_build_identity}")
if (NOT MSVC)
    target_compile_options(LevelScopeStandaloneSupport PRIVATE -fno-fast-math -ffp-contract=off)
    target_compile_options(LevelScopeMeasurementCore PRIVATE -fno-fast-math)
endif()
target_compile_definitions(LevelScopeStandaloneSupport PUBLIC JUCE_WEB_BROWSER=0 JUCE_USE_CURL=0)
target_link_libraries(LevelScopeStandaloneSupport PUBLIC
    LevelScopeMeasurementCore juce::juce_audio_formats juce::juce_audio_devices juce::juce_cryptography juce::juce_graphics)

add_library(LevelScopeDiagnosticsUI STATIC
    Source/Standalone/StandaloneSectionsPanel.cpp
    Source/Standalone/MeasurementTimeline.cpp
    Source/Standalone/StandaloneDiagnosticsPanel.cpp)
target_link_libraries(LevelScopeDiagnosticsUI PUBLIC LevelScopeStandaloneSupport juce::juce_gui_basics)

juce_add_gui_app(LevelScopeStandalone
    PRODUCT_NAME "LevelScope Standalone"
    COMPANY_NAME "Studio Paradiso"
    BUNDLE_ID "com.studioparadiso.levelscope.standalone"
)

target_link_libraries(LevelScopeStandalone PRIVATE LevelScopeDiagnosticsUI)

target_sources(LevelScopeStandalone
    PRIVATE
        Source/Standalone/StandaloneSessionUI.cpp
        Source/Standalone/StandaloneApp.cpp
        Source/Standalone/StandaloneMainComponent.cpp
        Source/Standalone/StandaloneMainComponent.h
        Source/Standalone/StandaloneSessionModel.h
        Source/Standalone/StandaloneTargetProfiles.h
        Source/Standalone/StandaloneWorkflowReadiness.h
        Source/Standalone/StandaloneAudioPanel.cpp
        Source/Standalone/StandaloneAudioPanel.h
        Source/Standalone/StandaloneServices.h
        Source/Standalone/StandaloneVideoPanel.cpp
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
        LevelScopeStandaloneSupport
        juce::juce_gui_basics
)

if (MSVC)
    target_compile_options(LevelScopeStandalone PRIVATE /W4)
    target_compile_options(LevelScopeStandaloneSupport PRIVATE /W4)
else()
    target_compile_options(LevelScopeStandalone PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(LevelScopeStandaloneSupport PRIVATE -Wall -Wextra -Wpedantic)
endif()

option(LEVELSCOPE_BUILD_VIDEO_TESTS "Build hardware-independent video/transport tests" OFF)
if (LEVELSCOPE_BUILD_VIDEO_TESTS)
    enable_testing()
    juce_add_console_app(LevelScopeVideoTests PRODUCT_NAME "LevelScope Video Tests")
    target_sources(LevelScopeVideoTests PRIVATE Tests/VideoTests.cpp)
    if (APPLE)
        target_sources(LevelScopeVideoTests PRIVATE Tests/VideoReaderAudit_mac.mm)
        target_link_libraries(LevelScopeVideoTests PRIVATE "-framework AVFoundation")
    endif()
    target_link_libraries(LevelScopeVideoTests PRIVATE LevelScopeStandaloneSupport)
    add_test(NAME video-timing COMMAND LevelScopeVideoTests)
    set_tests_properties(video-timing PROPERTIES TIMEOUT 60)
    # Manual opt-in only: never register a physical audio test with CTest.
    juce_add_gui_app(LevelScopeVideoQualification PRODUCT_NAME "LevelScope Video Qualification"
        BUNDLE_ID "com.studioparadiso.levelscope.videoqualification")
    target_sources(LevelScopeVideoQualification PRIVATE Tests/VideoPlaybackQualification.cpp)
    target_link_libraries(LevelScopeVideoQualification PRIVATE LevelScopeStandaloneSupport juce::juce_gui_basics)
endif()

# Explicit diagnostic target, never a replacement for the frozen qualification harness.
option(LEVELSCOPE_BUILD_VIDEO_PAINT_DIAGNOSTIC "Build bounded video paint tracing tests" OFF)
if (LEVELSCOPE_BUILD_VIDEO_PAINT_DIAGNOSTIC)
    if (NOT APPLE)
        message(FATAL_ERROR "Video paint diagnostic currently requires macOS")
    endif()
    juce_add_gui_app(LevelScopeVideoPaintDiagnostic PRODUCT_NAME "LevelScope Video Paint Diagnostic"
        BUNDLE_ID "com.studioparadiso.levelscope.videopaintdiagnostic")
    target_sources(LevelScopeVideoPaintDiagnostic PRIVATE Tests/VideoPlaybackQualification.cpp Tests/VideoPaint/Trace_mac.mm)
    target_compile_definitions(LevelScopeVideoPaintDiagnostic PRIVATE LEVELSCOPE_VIDEO_PAINT_DIAGNOSTICS=1)
    target_link_libraries(LevelScopeVideoPaintDiagnostic PRIVATE LevelScopeStandaloneSupport juce::juce_gui_basics)
    add_executable(LevelScopeVideoPaintTraceTests Tests/VideoPaint/TraceBufferTests.cpp)
endif()

option(LEVELSCOPE_BUILD_MEASUREMENT_TESTS "Build hardware-independent offline measurement checks" OFF)
if (LEVELSCOPE_BUILD_MEASUREMENT_TESTS)
    enable_testing()
    juce_add_console_app(LevelScopeMeasurementTests PRODUCT_NAME "LevelScope Measurement Tests")
    target_sources(LevelScopeMeasurementTests PRIVATE Tests/OfflineMeasurementTests.cpp)
    target_link_libraries(LevelScopeMeasurementTests PRIVATE LevelScopeStandaloneSupport)
    add_test(NAME offline-measurement COMMAND LevelScopeMeasurementTests)
    set_tests_properties(offline-measurement PROPERTIES TIMEOUT 300)
endif()

option(LEVELSCOPE_BUILD_COMPLIANCE_TESTS "Build M2 profile/evaluation checks" OFF)
if (LEVELSCOPE_BUILD_COMPLIANCE_TESTS)
    enable_testing()
    juce_add_console_app(LevelScopeComplianceTests PRODUCT_NAME "LevelScope Compliance Tests")
    target_sources(LevelScopeComplianceTests PRIVATE Tests/ComplianceTests.cpp)
    target_link_libraries(LevelScopeComplianceTests PRIVATE LevelScopeStandaloneSupport)
    add_test(NAME measured-compliance COMMAND LevelScopeComplianceTests)
    set_tests_properties(measured-compliance PROPERTIES TIMEOUT 60)
endif()

option(LEVELSCOPE_BUILD_PROPOSAL_TESTS "Build M3 global-trim checks" OFF)
if (LEVELSCOPE_BUILD_PROPOSAL_TESTS)
    enable_testing()
    juce_add_console_app(LevelScopeProposalTests PRODUCT_NAME "LevelScope Proposal Tests")
    target_sources(LevelScopeProposalTests PRIVATE Tests/ProposalTests.cpp)
    target_link_libraries(LevelScopeProposalTests PRIVATE LevelScopeStandaloneSupport)
    add_test(NAME global-trim-proposals COMMAND LevelScopeProposalTests)
    set_tests_properties(global-trim-proposals PROPERTIES TIMEOUT 120)
endif()

option(LEVELSCOPE_BUILD_DIAGNOSTICS_TESTS "Build M4 measured diagnostics checks" OFF)
if (LEVELSCOPE_BUILD_DIAGNOSTICS_TESTS)
    enable_testing()
    juce_add_console_app(LevelScopeDiagnosticsTests PRODUCT_NAME "LevelScope Diagnostics Tests")
    target_sources(LevelScopeDiagnosticsTests PRIVATE Tests/DiagnosticsTests.cpp)
    target_link_libraries(LevelScopeDiagnosticsTests PRIVATE LevelScopeDiagnosticsUI)
    add_test(NAME measured-diagnostics COMMAND LevelScopeDiagnosticsTests)
    set_tests_properties(measured-diagnostics PROPERTIES TIMEOUT 120)
endif()

option(LEVELSCOPE_BUILD_SECTIONS_TESTS "Build M5 programme section checks" OFF)
if (LEVELSCOPE_BUILD_SECTIONS_TESTS)
    enable_testing()
    juce_add_console_app(LevelScopeSectionsTests PRODUCT_NAME "LevelScope Sections Tests")
    target_sources(LevelScopeSectionsTests PRIVATE Tests/SectionsTests.cpp)
    target_link_libraries(LevelScopeSectionsTests PRIVATE LevelScopeDiagnosticsUI)
    add_test(NAME programme-sections COMMAND LevelScopeSectionsTests)
    set_tests_properties(programme-sections PROPERTIES TIMEOUT 120)
endif()

option(LEVELSCOPE_BUILD_AUDITION_TESTS "Build M6 verified audition checks" OFF)
if (LEVELSCOPE_BUILD_AUDITION_TESTS)
    enable_testing()
    juce_add_console_app(LevelScopeAuditionTests PRODUCT_NAME "LevelScope Audition Tests")
    target_sources(LevelScopeAuditionTests PRIVATE Tests/AuditionTests.cpp Source/Standalone/StandaloneAudioPanel.cpp)
    target_link_libraries(LevelScopeAuditionTests PRIVATE LevelScopeStandaloneSupport juce::juce_gui_basics)
    add_test(NAME verified-audition COMMAND LevelScopeAuditionTests)
    set_tests_properties(verified-audition PROPERTIES TIMEOUT 300)
endif()

option(LEVELSCOPE_BUILD_MULTICHANNEL_TESTS "Build M7 System I PCM 7.1 checks" OFF)
if (LEVELSCOPE_BUILD_MULTICHANNEL_TESTS)
    enable_testing()
    juce_add_console_app(LevelScopeMultichannelTests PRODUCT_NAME "LevelScope Multichannel Tests")
    target_sources(LevelScopeMultichannelTests PRIVATE Tests/MultichannelTests.cpp Source/Standalone/StandaloneAudioPanel.cpp)
    if (APPLE)
        target_sources(LevelScopeMultichannelTests PRIVATE Tests/VideoReaderAudit_mac.mm)
    endif()
    target_link_libraries(LevelScopeMultichannelTests PRIVATE LevelScopeDiagnosticsUI)
    add_test(NAME multichannel-pcm COMMAND LevelScopeMultichannelTests)
    set_tests_properties(multichannel-pcm PROPERTIES TIMEOUT 600)
endif()

option(LEVELSCOPE_BUILD_SCAFFOLD_TESTS "Build hardware-independent standalone checks" OFF)
if (LEVELSCOPE_BUILD_SCAFFOLD_TESTS)
    enable_testing()
    juce_add_console_app(LevelScopeScaffoldTests PRODUCT_NAME "LevelScope Scaffold Tests")
    target_sources(LevelScopeScaffoldTests PRIVATE Tests/StandaloneScaffoldTests.cpp)
    target_link_libraries(LevelScopeScaffoldTests PRIVATE LevelScopeStandaloneSupport)
    add_test(NAME standalone-scaffold COMMAND LevelScopeScaffoldTests)
    set_tests_properties(standalone-scaffold PROPERTIES TIMEOUT 60)
endif()

option(LEVELSCOPE_BUILD_SESSION_TESTS "Build M8 reproducible-session checks" OFF)
if (LEVELSCOPE_BUILD_SESSION_TESTS)
    enable_testing()
    juce_add_console_app(LevelScopeSessionTests PRODUCT_NAME "LevelScope Session Tests")
    target_sources(LevelScopeSessionTests PRIVATE Tests/SessionTests.cpp)
    target_compile_definitions(LevelScopeSessionTests PRIVATE JUCE_MODAL_LOOPS_PERMITTED=1)
    target_link_libraries(LevelScopeSessionTests PRIVATE LevelScopeStandaloneSupport)
    add_test(NAME reproducible-sessions COMMAND LevelScopeSessionTests)
    set_tests_properties(reproducible-sessions PROPERTIES TIMEOUT 300)
    juce_add_console_app(LevelScopeSessionPlaybackTests PRODUCT_NAME "LevelScope Session Playback Tests")
    target_sources(LevelScopeSessionPlaybackTests PRIVATE Tests/SessionPlaybackTests.cpp Source/Standalone/StandaloneAudioPanel.cpp)
    if (APPLE)
        target_sources(LevelScopeSessionPlaybackTests PRIVATE Tests/VideoReaderAudit_mac.mm)
    endif()
    target_link_libraries(LevelScopeSessionPlaybackTests PRIVATE LevelScopeStandaloneSupport juce::juce_gui_basics)
    add_test(NAME restored-session-audition COMMAND LevelScopeSessionPlaybackTests)
    set_tests_properties(restored-session-audition PROPERTIES TIMEOUT 300)
endif()
