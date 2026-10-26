#pragma once
namespace
{
[[noreturn]] void crashCase(const juce::File& root,const juce::String& stage)
{
    expect(root.createDirectory().wasOk(),"crash fixture directory");auto s=measured(root,2,48000);const auto request=ex::capture(s);JobControl c;pcm::Hooks hook;
    hook.reached=[&](const char* name){if(stage==name)::_exit(86);};
    auto prepared=prepareVerifiedSource(request.current,c);expect(prepared.value.has_value(),"crash owned copy prepared");
    juce::StringArray owned;for(const auto& file:(*prepared.value)->files())owned.add(file.getFullPathName());
    expect(root.getChildFile("owned-copies.txt").replaceWithText(owned.joinIntoString("\n")),"record exact crash-owned copies for test cleanup");
    auto result=ex::render(request,root,c,hook,nullptr,*prepared.value);expect(result.value.has_value(),"crash fixture rendered");
    const auto final=root.getChildFile("final");expect(ex::publish(**result.value,request,final,ex::approvalIdentity(**result.value,request),c,hook).value.has_value(),"crash fixture committed");
    if(stage=="committed")::_exit(86);::_exit(87);
}
void crashTests()
{
    Temp temp;const auto executable=juce::File::getSpecialLocation(juce::File::currentExecutableFile);
    for(const char* stage:{"audio-write","publish-rename","committed"})
    {
        const auto folder=temp.directory.getChildFile(stage);juce::ChildProcess child;
        expect(child.start(juce::StringArray{executable.getFullPathName(),"--crash-folder",folder.getFullPathName(),"--crash-stage",stage}),"start isolated crash simulation");
        expect(child.waitForProcessToFinish(30000)&&child.getExitCode()==86,"abrupt worker process exit without destructors");
        const auto final=folder.getChildFile("final");const bool committed=juce::String(stage)=="committed";
        expect(final.exists()==committed,"crash exposes final bundle only after atomic commit");
        if(committed)
        {
            JobControl c;const auto report=juce::JSON::parse(final.getChildFile("report.json"));
            expect(report["artifact"]["sha256"].toString()==pcm::hashFile(final.getChildFile("programme.wav"),c),"post-commit crash retains exact complete published bytes");
        }
        else expect(!folder.findChildFiles(juce::File::findDirectories,false,".levelscope-export-*").isEmpty(),"crash staging is identifiable and never auto-published");
        const auto survivor=folder.getChildFile("unowned.txt");expect(survivor.replaceWithText("preserve"),"unowned recovery sentinel");
        juce::StringArray owned;owned.addLines(folder.getChildFile("owned-copies.txt").loadFileAsString());
        for(const auto& path:owned)if(path.isNotEmpty())expect(juce::File(path).deleteFile(),"remove only child-recorded owned snapshot after crash evidence");
        StandaloneSessionModel reopened;expect(!reopened.sourceDocument.hasSource&&!reopened.currentStateDocument.isInitialized&&survivor.loadFileAsString()=="preserve","new runtime inherits no export job/approval and deletes no unowned files");
    }
}
}
