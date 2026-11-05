// Supplemental end-to-end population grouping + envelope-selection stress.
#include "Core/Proposal/LraPassages.h"
#include <chrono>
#include <iostream>
#include <mach/mach.h>
int main(int argc,char** argv)
{
    if(argc!=2||juce::File(argv[1]).exists())return 2;
    using namespace levelscope;using Clock=std::chrono::steady_clock;
    auto a=std::make_shared<proposal::ProgramAnalysis>();a->sampleRate=192000;a->frames=int64_t(a->sampleRate)*10800;a->identity="three-hour-sized-population";a->p10=-35;a->p95=-20;a->lra={15,measurement::Status::available};a->lraStable=true;
    for(int n=0;n<107985;++n){proposal::DiagnosticWindow w;w.begin=int64_t(n)*19200;w.end=w.begin+576000;w.included=true;w.lufs=n%100<10?-20.:-35.;a->lraWindows.push_back(w);}
    processing::SectionGainPlan before;before.inputIdentity=juce::String::repeatedString("a",64);before.sampleRate=a->sampleRate;before.frames=a->frames;before.intervals={{0,a->frames,0}};
    const auto start=Clock::now();auto population=proposal::findLraPassages(a,10);const auto grouped=Clock::now();auto selected=proposal::selectLraTrims(*population,before,{});const auto end=Clock::now();
    const double grouping=std::chrono::duration<double>(grouped-start).count(),planning=std::chrono::duration<double>(end-grouped).count();
    mach_task_basic_info_data_t info{};mach_msg_type_number_t count=MACH_TASK_BASIC_INFO_COUNT;task_info(mach_task_self(),MACH_TASK_BASIC_INFO,reinterpret_cast<task_info_t>(&info),&count);
    auto out=juce::var(new juce::DynamicObject());auto* o=out.getDynamicObject();o->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);o->setProperty("grouping_seconds",grouping);o->setProperty("planning_seconds",planning);o->setProperty("resident_bytes",juce::int64(info.resident_size));o->setProperty("selected_regions",int(selected.selected.size()));o->setProperty("available_regions",int(population->high.size()));o->setProperty("affected_frames",juce::int64(selected.affectedFrames));o->setProperty("cancelled_cleanly",!proposal::findLraPassages(a,10,[]{return true;}));
#if defined(__arm64__) || defined(__aarch64__)
    const double budget=2;
#else
    const double budget=5;
#endif
    const bool passed=grouping+planning<=budget&&info.resident_size<64*1024*1024&&selected.available;
    o->setProperty("passed",passed);o->setProperty("budget_seconds",budget);o->setProperty("full_duration_audio_decoded",false);
    std::cout<<juce::JSON::toString(out)<<'\n';return juce::File(argv[1]).replaceWithText(juce::JSON::toString(out))?(passed?0:1):2;
}
