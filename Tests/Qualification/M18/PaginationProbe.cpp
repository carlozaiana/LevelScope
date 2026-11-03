#include "AuditionTestSupport.h"
#include "Standalone/StandaloneProposalPreferences.h"
#include <algorithm>
int main(int argc,char** argv){
    if(argc!=2)return 2;juce::ScopedJuceInitialiser_GUI gui;juce::Array<juce::var> rows;
    for(int trial=0;trial<3;++trial){
        juce::MemoryOutputStream text;for(int n=0;n<10000;++n)text<<"Changed annotation "<<n<<"; explicit proposal details; no omitted rows.\n";const auto all=text.toUTF8();
        struct Pulse:juce::Timer{std::vector<double> ticks;double last=juce::Time::getMillisecondCounterHiRes();void timerCallback()override{const auto now=juce::Time::getMillisecondCounterHiRes();ticks.push_back(now-last);last=now;}} pulse;
        pulse.startTimer(5);juce::MessageManager::getInstance()->runDispatchLoopUntil(40);
        std::vector<double> pages;ProposalInspectionPanel panel(all);
        for(int n=0;n<panel.pageCount();++n){const auto t=juce::Time::getMillisecondCounterHiRes();panel.nextPage();pages.push_back(juce::Time::getMillisecondCounterHiRes()-t);juce::MessageManager::getInstance()->runDispatchLoopUntil(5);}
        juce::MessageManager::getInstance()->runDispatchLoopUntil(100);pulse.stopTimer();std::sort(pages.begin(),pages.end());std::sort(pulse.ticks.begin(),pulse.ticks.end());
        auto row=juce::var(new juce::DynamicObject());auto* r=row.getDynamicObject();r->setProperty("trial",trial+1);r->setProperty("pages",panel.pageCount());r->setProperty("complete",panel.allText()==all);r->setProperty("page_p95_ms",pages[size_t(.95*(pages.size()-1))]);r->setProperty("page_max_ms",pages.back());r->setProperty("heartbeat_p95_ms",pulse.ticks[size_t(.95*(pulse.ticks.size()-1))]);r->setProperty("heartbeat_max_ms",pulse.ticks.back());rows.add(row);std::cout<<juce::JSON::toString(row,true)<<std::endl;
    }
    auto v=juce::var(new juce::DynamicObject());v.getDynamicObject()->setProperty("trials",rows);v.getDynamicObject()->setProperty("base_build",LEVELSCOPE_MEASUREMENT_BUILD);v.getDynamicObject()->setProperty("os",juce::SystemStats::getOperatingSystemName());juce::File(argv[1]).replaceWithText(juce::JSON::toString(v,true));return 0;
}
