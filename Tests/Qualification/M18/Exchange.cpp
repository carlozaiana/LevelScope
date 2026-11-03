#define main unusedM13QualificationMain
#include "SectionProposalTests.cpp"
#undef main
namespace {
void saveExchange(StandaloneSessionModel& s,const juce::File& file){JobControl job;expect(sf::save(sf::capture(s),file,job).succeeded(),"save exact exchange fixture");}
void prepareExchange(const juce::File& from,const juce::File& root){
    for(auto name:{"stereo","surround"}){
        JobControl job;auto document=sf::read(from.getChildFile(juce::String("native-")+name+"/job.levelscope"),job);expect(bool(document.document),"read native accepted session");
        auto verified=sf::verify(*document.document,job);StandaloneSessionModel s;juce::String notice;
        expect(verified.prepared&&sf::adopt(*verified.document,*verified.prepared,s,notice)&&s.source.hasMeasurement&&s.currentState.hasMeasurement&&s.lastAppliedExplanation,"prepare independently verified schema-8 exchange");
        const auto folder=root.getChildFile(name);folder.createDirectory();expect(s.sourceDocument.file.copyFileTo(folder.getChildFile(s.sourceDocument.file.getFileName())),"copy only owned short exchange PCM");
        saveExchange(s,folder.getChildFile("verified.levelscope"));
        auto edited=s;auto req=proposal::JobRequirements::editable(s.getSelectedTargetProfile());
        req.maximumLraTenths=180;req.maximumMomentaryTenths=-150;req.maximumShortTermTenths=-180;
        expect(edited.setJobRequirements(req,notice),"explicit full target fields");saveExchange(edited,folder.getChildFile("target-lra-m-s.levelscope"));
        req.integratedTenths.reset();req.truePeakTenths.reset();req.maximumLraTenths.reset();req.maximumMomentaryTenths.reset();req.maximumShortTermTenths.reset();
        expect(edited.setJobRequirements(req,notice),"explicit No Limit for every target");saveExchange(edited,folder.getChildFile("target-no-limit.levelscope"));
        expect(edited.setCurrentGainDb(edited.currentStateDocument.staticGainDb-.1)&&!edited.currentState.hasMeasurement,"later manual processing invalidates Current");saveExchange(edited,folder.getChildFile("processing-unverified.levelscope"));
    }
}
void readExchange(const juce::File& from,const juce::File& root){
    for(auto name:{"stereo","surround"})for(auto variant:{"verified","target-lra-m-s","target-no-limit","processing-unverified"}){
        const auto folder=from.getChildFile(name),file=folder.getChildFile(juce::String(variant)+".levelscope");JobControl job;auto read=sf::read(file,job);
        expect(read.document&&read.document->loadedSchema==8&&read.document->lastAppliedExplanation,"foreign schema-8 history is present");
        auto before=*read.document;const auto history=encodeAppliedExplanation(*before.lastAppliedExplanation);const auto plan=proposalPlanIdentity(before.current.processingPlan());
        const auto requirements=before.jobRequirements.digest(),preferences=before.proposalPolicy.digest();
        const auto geometry=before.sections.digest();const auto roles=before.source.layout;
        const bool previouslyVerified=juce::String(variant)!="processing-unverified";
        expect(bool(before.currentCache)==previouslyVerified,"sender valid versus stale cache matches fixture");
        const auto origin=before.sourceCache->context.platformAndFloatingPointMode;
        expect(sf::relinkAudio(*read.document,folder.getChildFile(before.source.file.getFileName())),"foreign audio relink is explicit");
        auto verified=sf::verify(*read.document,job);StandaloneSessionModel s;juce::String notice;
        expect(verified.prepared&&sf::adopt(*verified.document,*verified.prepared,s,notice),"foreign source verified and session adopted");
        const auto receiver=s.measurementContext(MeasurementDomain::source).platformAndFloatingPointMode;
        expect(origin!=receiver,"exchange actually crosses floating-point platform identity");
        expect(!s.source.hasMeasurement&&!s.currentState.hasMeasurement&&!s.sectionProposal,"foreign caches and private candidate never become verification");
        expect(s.lastAppliedExplanation&&encodeAppliedExplanation(*s.lastAppliedExplanation)==history,"complete historical explanation byte-exact across Macs");
        expect(proposalPlanIdentity(s.currentStateDocument.processingPlan())==plan&&s.sectionDocument.map().digest()==geometry,"section geometry and exact applied processing retained");
        expect(s.jobRequirements.digest()==requirements&&s.proposalPolicy.digest()==preferences,"all target fields and preferences exact");
        expect(s.sourceDocument.layout.count==roles.count&&s.sourceDocument.layout.roles==roles.roles&&s.sourceDocument.layout.confirmed,"saved channel roles retained");
        for(int route:s.sourceDocument.routing.fileToDevice)expect(route==-1,"reopened physical routes remain muted");
        if(juce::String(variant)=="verified"){
            measuredComposition(s,MeasurementDomain::source);measuredComposition(s,MeasurementDomain::currentState);
            expect(s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"foreign applied plan independently remeasured on receiving Mac");
        }else if(juce::String(variant)=="processing-unverified")expect(s.appliedExplanationStatus().contains("processing edited"),"changed processing visibly historical");
        else expect(s.appliedExplanationStatus().contains("Requirements changed"),"changed target visibly historical");
        const auto destination=root.getChildFile(juce::String(name)+"-"+variant+".levelscope");saveExchange(s,destination);
        auto reread=sf::read(destination,job);auto reopened=sf::verify(*reread.document,job);
        expect(reopened.prepared&&reopened.prepared->model.currentState.hasMeasurement==(juce::String(variant)=="verified"),"receiving-platform save/reopen preserves independent verification state");
        auto row=juce::var(new juce::DynamicObject());auto* r=row.getDynamicObject();r->setProperty("layout",name);r->setProperty("variant",variant);r->setProperty("origin",origin);r->setProperty("receiver",receiver);r->setProperty("history_sha256",juce::SHA256(history.data(),history.size()).toHexString());r->setProperty("current_remeasured",s.currentState.hasMeasurement);exchangeEvidence.add(row);
    }
}
}
int main(int argc,char** argv){
    if(argc!=4)return 2;juce::ScopedJuceInitialiser_GUI gui;juce::File root(argv[3]);bool passed=false;juce::String error;
    try{expect(!root.exists()&&root.createDirectory().wasOk(),"fresh exchange evidence directory");if(juce::String(argv[1])=="prepare")prepareExchange(juce::File(argv[2]),root);else if(juce::String(argv[1])=="read")readExchange(juce::File(argv[2]),root);else throw std::runtime_error("Unknown exchange action");passed=true;}
    catch(const std::exception& e){error=e.what();std::cerr<<error<<std::endl;}
    auto v=juce::var(new juce::DynamicObject());auto* r=v.getDynamicObject();r->setProperty("passed",passed);r->setProperty("failure",error);r->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);r->setProperty("checks",checks);r->setProperty("exchange_cases",exchangeEvidence);root.getChildFile("result.json").replaceWithText(juce::JSON::toString(v,true));return passed?0:1;
}
