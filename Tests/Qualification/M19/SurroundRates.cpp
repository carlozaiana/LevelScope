// Complete the remaining 5.1 rate cells with the frozen workflow assertions.
#define main m19_original_test_main
#include "../../LraProposalTests.cpp"
#undef main
int main(int argc,char** argv)
{
    if(argc!=3)return 2;juce::ScopedJuceInitialiser_GUI gui;juce::File root(argv[1]);bool passed=false;juce::String failure;
    try
    {
        expect(!root.exists()&&root.createDirectory().wasOk(),"fresh supplemental 5.1 rate root");
        for(int rate:{44100,96000,192000})
        {workflow(root.getChildFile(juce::String(rate)),6,rate);cases.getReference(cases.size()-1).getDynamicObject()->setProperty("rate",rate);}
        passed=true;
    }
    catch(const std::exception& e){failure=e.what();}
    auto out=juce::var(new juce::DynamicObject());auto* o=out.getDynamicObject();o->setProperty("passed",passed);o->setProperty("failure",failure);o->setProperty("checks",checks);o->setProperty("cases",cases);o->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);o->setProperty("physical_devices_opened",false);o->setProperty("video_loaded",false);
    return juce::File(argv[2]).replaceWithText(juce::JSON::toString(out,true))?(passed?0:1):2;
}
