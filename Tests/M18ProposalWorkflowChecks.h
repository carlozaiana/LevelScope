#pragma once
#include "Standalone/ProposalWorkflowPanel.h"
#include "Standalone/TimelineGainLane.h"
#include "Standalone/MeasurementTimeline.h"
namespace
{
juce::TextButton& proposalButton(juce::Component& ui,const juce::String& name)
{
    for(auto* c:ui.getChildren())if(auto* b=dynamic_cast<juce::TextButton*>(c);b&&b->getButtonText()==name)return *b;
    throw std::runtime_error("Missing proposal button "+name.toStdString());
}
void proposalImage(juce::Component& ui,const juce::File& path)
{auto out=path.createOutputStream();expect(out&&juce::PNGImageFormat().writeImageToStream(ui.createComponentSnapshot(ui.getLocalBounds()),*out),"save proposal UI capture");}
std::shared_ptr<const SectionProposalResult> proposeFor(StandaloneSessionModel& s)
{
    JobControl c;auto result=generateSectionProposal(s.sectionProposalRequest(),c);
    expect(result.value&&s.acceptSectionProposal(*result.value),"complete existing global-first section-led service");return *result.value;
}
void limitsWindowAfterPropose(StandaloneSessionModel& s, ProposalWorkflowPanel& panel)
{
    auto content = std::make_unique<juce::Component>(); content->setSize(900,650);
    ProposalWindow owner("Test timeline",std::move(content),nullptr); owner.present();
    std::unique_ptr<ProposalWindow> limits;
    panel.onPreferences = [&]{ ProposalPreferencesPanel::showWindow(limits,s,&owner,[]{}); };
    auto& open = proposalButton(panel,"Proposal limits...");
    const auto pump = []{ juce::MessageManager::getInstance()->runDispatchLoopUntil(60); };
    open.onClick(); pump();
    expect(limits && limits->getPeer()->isFocused(), "opening limits after Propose brings the native editor to the front");
    auto* form = limits->getContentComponent();
    juce::TextEditor* duration = nullptr;
    for (auto* child : form->getChildren())
        if (child->getName().startsWith("Newly affected duration")) duration = dynamic_cast<juce::TextEditor*>(child);
    expect(duration!=nullptr, "proposal limits exposes editable duration");
    duration->setText("32.1",false);
    owner.present(); pump();
    expect(owner.getPeer()->isFocused(), "timeline can move in front of the modeless limits editor");
    open.onClick(); pump();
    expect(limits->getContentComponent()==form && limits->getPeer()->isFocused() && duration->getText()=="32.1"
           && s.proposalPolicy.durationThousandths==200, "repeat click raises the same native editor with unconfirmed text intact");
    limits->closeButtonPressed();
    open.onClick(); pump();
    expect(limits->getContentComponent()!=form && limits->getPeer()->isFocused()
           && proposalButton(panel,"Apply proposal").isEnabled(), "closing and reopening starts a fresh editor and leaves the proposal actionable");
    panel.onPreferences = {};
}
// Exercise the same repeated control/lane/dirty checks found in the live hang.
// Inflate only the in-memory window count: this is a UI fixture, not a
// full-duration decode, delivery measurement or qualification run.
void longWindowRefresh(StandaloneSessionModel s)
{
    for (auto* summary : { &s.source, &s.currentState })
    {
        auto value = *summary->analysis;
        auto fill = [&] (std::vector<measurement::Window>& windows, int count)
        {
            windows.resize(count);
            for (int n = 0; n < count; ++n)
                windows[n] = { value.detail.frames * (n + 1) / count, .0123 + n * 1.e-9 };
        };
        fill(value.detail.momentary, 36000); fill(value.detail.shortTerm, 36000); fill(value.detail.lraWindows, 3600);
        summary->analysis = std::make_shared<const AnalysisSnapshot>(std::move(value));
    }
    const auto result = proposeFor(s);
    expect(result->status()==SectionProposalStatus::globalOnly, "large window payload has an actionable prediction-only proposal");
    const auto request = s.sectionProposalRequest();
    const auto start = juce::Time::getMillisecondCounterHiRes();
    const auto uncached = trimRequest(request.source.get(),request.current.get(),request.sourceContext,request.currentContext,
                                      request.currentContext.profile,request.applicability);
    const auto hashMs = juce::Time::getMillisecondCounterHiRes() - start;
    const auto cached = request.scalar();
    expect(cached.source.snapshotIdentity==uncached.source.snapshotIdentity && cached.current.snapshotIdentity==uncached.current.snapshotIdentity,
           "published Source and Current preserve the complete canonical measurement digests");
    ProposalWorkflowPanel panel(s); TimelineGainLane lane(s); lane.synchronize("long-windows");
    std::vector<double> times;
    for (int n = 0; n < 20; ++n)
    {
        const auto t = juce::Time::getMillisecondCounterHiRes();
        panel.refresh(); lane.synchronize("long-windows");
        const bool actionable = s.hasActionableSectionProposal();
        times.push_back(juce::Time::getMillisecondCounterHiRes() - t);
        expect(actionable && proposalButton(panel,"Apply proposal").isEnabled(), "live refresh retains the actionable proposal");
    }
    std::sort(times.begin(),times.end());
    expect(times[18]<25, "long-window proposal UI refresh p95 stays below a quarter of its 100 ms tick");
    auto metrics = juce::var(new juce::DynamicObject());
    metrics.getDynamicObject()->setProperty("long_window_proposal_refresh_p95_ms",times[18]);
    metrics.getDynamicObject()->setProperty("uncached_selected_measurement_pair_ms",hashMs); cases.add(metrics);
    limitsWindowAfterPropose(s,panel);

    auto changed = *s.currentState.analysis;
    changed.detail.shortTerm.back().energy = std::nextafter(changed.detail.shortTerm.back().energy, 1.);
    s.currentState.analysis = std::make_shared<const AnalysisSnapshot>(std::move(changed));
    panel.refresh(); lane.synchronize("long-windows");
    expect(s.sectionProposalRequest().scalar().current.snapshotIdentity!=cached.current.snapshotIdentity
           && panel.outcomeText().contains("Stale inputs") && !proposalButton(panel,"Apply proposal").isEnabled()
           && !s.applySectionProposal(), "a new measurement with just one last-window bit changed invalidates the preview and Apply");
    auto wrong = request.sourceContext; wrong.sessionId = "foreign-session";
    expect(!measuredProgrammeForSnapshot(request.source, wrong, MeasurementDomain::source, request.applicability).currentIdentityMatches,
           "digest reuse does not reuse a context match for a different session");
    std::weak_ptr<const AnalysisSnapshot> old = s.currentState.analysis;
    s.currentState.analysis.reset();
    expect(old.expired(), "derived digests do not retain old measurement arrays");
}
void historicalRoundTrip(StandaloneSessionModel& s,const juce::File& path)
{
    const auto history=s.lastAppliedExplanation;expect(bool(history),"successful Apply captures a complete historical explanation");
    JobControl io;const auto doc=sf::capture(s);expect(doc.lastAppliedExplanation==history,"save captures immutable explanation with one document revision");
    expect(sf::save(doc,path,io).succeeded(),"save applied explanation attachment");const auto read=sf::read(path,io);
    expect(read.document&&read.document->lastAppliedExplanation&&read.document->lastAppliedExplanation->text()==history->text(),"complete historical explanation survives binary session round trip");
    auto verified=sf::verify(*read.document,io);StandaloneSessionModel reopened;juce::String notice;
    expect(verified.prepared&&sf::adopt(*verified.document,*verified.prepared,reopened,notice),"adopt saved session with historical explanation after media verification");
    expect(reopened.lastAppliedExplanation&&reopened.lastAppliedExplanation->text()==history->text()&&!reopened.sectionProposal&&!reopened.currentState.hasMeasurement&&reopened.currentStateNeedsRemeasurement(),"private candidate history never becomes reopened Current verification");
    expect(reopened.appliedExplanationStatus().contains("for this plan")&&reopened.appliedExplanationStatus().contains("unverified"),"reopen identifies plan association and independent unverified status");
    for(int route:reopened.sourceDocument.routing.fileToDevice)expect(route==-1,"history reopening keeps physical outputs muted");
    const auto same=proposalPlanIdentity(reopened.currentStateDocument.processingPlan());
    expect(reopened.setCurrentGainDb(reopened.currentStateDocument.staticGainDb-.1)&&reopened.appliedExplanationStatus().contains("processing edited"),"manual gain retains explanation and marks its processing association historical");
    expect(reopened.setCurrentGainDb(history->after.staticGainDb)&&proposalPlanIdentity(reopened.currentStateDocument.processingPlan())==same&&reopened.currentStateNeedsRemeasurement(),"restoring exact historical processing does not restore verification");
    auto target=reopened.jobRequirements;target.maximumLraTenths=180;expect(reopened.setJobRequirements(target,notice)&&reopened.appliedExplanationStatus().contains("Requirements changed"),"requirements changes independently label historical explanation");
    auto prefs=reopened.proposalPolicy;prefs.durationThousandths=900;expect(reopened.setProposalPolicy(prefs,notice)&&reopened.appliedExplanationStatus().contains("preferences changed"),"preferences changes independently label historical explanation");
    const auto before=history->text();expect(reopened.sectionDocument.rename(1,"Renamed after application",notice)&&reopened.lastAppliedExplanation->text()==before,"historical names and intervals survive subsequent annotation edits");
    auto report=ex::report(ex::capture(reopened));expect(ex::reportJson(report).contains("Historical explanation")&&ex::reportHtml(report).contains("Historical explanation"),"reports expose saved reasons as escaped historical evidence");
    auto bytes=encodeAppliedExplanation(*history);JobControl decode;expect(decodeAppliedExplanation(bytes,decode)->text()==before,"independent explanation codec retains exact scalar values and structured reasons");
    const auto exact=decodeAppliedExplanation(bytes,decode);
    for(size_t n=0;n<history->baseline.metrics.size();++n)expect(exact->baseline.metrics[n].value==history->baseline.metrics[n].value&&exact->baseline.metrics[n].status==history->baseline.metrics[n].status,"historical baseline scalars and statuses round-trip exactly");
    expect(exact->beforeIdentity==history->beforeIdentity&&exact->afterIdentity==history->afterIdentity&&exact->requirements.digest()==history->requirements.digest()&&exact->target.digest()==history->target.digest(),"historical full target and exact plan identities survive serialization");
    auto replaced=reopened;expect(replaced.setValidatedSource(s.sourceDocument.file,s.sourceDocument.media)&&!replaced.lastAppliedExplanation,"replacing Source clears the prior programme explanation");
    auto bad=path.getSiblingFile(path.getFileNameWithoutExtension()+"-damaged.levelscope");juce::MemoryBlock file;expect(path.loadFileAsData(file),"load owned history corruption fixture");
    auto* raw=static_cast<uint8_t*>(file.getData());raw[file.getSize()-1]^=1;expect(bad.replaceWithData(raw,file.getSize()),"damage independent explanation bytes only");auto recovered=sf::read(bad,io);
    expect(recovered.document&&!recovered.document->lastAppliedExplanation&&recovered.document->sourceCache&&recovered.document->current.staticGainDb==doc.current.staticGainDb&&recovered.document->notices.joinIntoString(" ").contains("explanation discarded"),"damaged history recovers core processing and independent Source cache with visible notice");
    auto old=path.getSiblingFile(path.getFileNameWithoutExtension()+"-old.levelscope");rewritePolicyCore(path,old,[](juce::var& v){v.getDynamicObject()->removeProperty("appliedExplanation");});auto legacy=sf::read(old,io);
    expect(legacy.document&&!legacy.document->lastAppliedExplanation&&legacy.document->notices.joinIntoString(" ").contains("No reasons were reconstructed"),"older session without complete history migrates without invented reasons");
    auto foreign=*history;foreign.programmeId=juce::Uuid().toString();auto mismatched=doc;mismatched.lastAppliedExplanation=std::make_shared<const AppliedProposalExplanation>(foreign);
    auto wrong=path.getSiblingFile("wrong-history.levelscope");expect(sf::save(mismatched,wrong,io).succeeded(),"write intentionally mismatched documentary fixture");auto wrongRead=sf::read(wrong,io);
    expect(wrongRead.document&&!wrongRead.document->lastAppliedExplanation,"unrelated historical programme is discarded without attaching reasons to current audio");
    sf::FileHooks late;late.reached=[&](const char* step){if(juce::String(step)=="replace")s.documentChanged();};
    expect(sf::save(doc,path,io,late).succeeded()&&sf::saved(s,doc,path)&&s.documentDirty(),"late edit during save stays unsaved while coherent history revision is written");
}
void m18ProposalWorkflow(const juce::File& root)
{
    auto s=proposalFixture(root.getChildFile("sections"),48000,2,0);const auto baseline=s.currentState.analysis;const auto revision=s.revision;
    ProposalWorkflowPanel panel(s);panel.setSize(1200,ProposalWorkflowPanel::height);bool reviewed=false,historyOpened=false;
    panel.onPropose=[&]{proposeFor(s);panel.refresh();};panel.onReview=[&]{reviewed=true;};panel.onHistory=[&]{historyOpened=true;};
    expect(proposalButton(panel,"Propose").isEnabled()&&!proposalButton(panel,"Apply proposal").isEnabled()&&panel.actionText().contains("No ready proposal"),"main proposal actions expose disabled Apply reason");
    proposalButton(panel,"Propose").onClick();expect(s.sectionProposal&&s.sectionProposal->status()==SectionProposalStatus::ready&&s.currentState.analysis==baseline&&s.revision==revision,"main Propose retains measured baseline and applied settings");
    const auto candidate=s.sectionProposal;expect(panel.outcomeText().contains("Ready to apply")&&proposalButton(panel,"Review proposal").isEnabled()&&proposalButton(panel,"Apply proposal").isEnabled(),"main outcome prominently offers review and explicit Apply");
    proposalButton(panel,"Review proposal").onClick();expect(reviewed,"main Review proposal callback accessible");
    TimelineGainLane lane(s);MeasurementTimeline timeline;timeline.setSize(1200,520);timeline.attachGainLane(lane);timeline.setProgramme(s.sourceDocument.media.frames,48000,"proposal-lane");lane.synchronize("proposal-lane");
    const auto& correction=candidate->selection().corrections.front();lane.selectFrame(correction.begin,true);
    expect(lane.inspectionText().contains("Proposed local")&&lane.inspectionText().contains("global")&&lane.inspectionText().contains("combined")&&lane.inspectionText().contains("Reserve achieved")&&s.currentState.analysis==baseline,"automatic lane inspection distinguishes local/global/combined gain and complete reason without applying");
    proposalImage(panel,root.getChildFile("ready-proposal.png"));proposalImage(timeline,root.getChildFile("proposed-gain-lane.png"));
    auto preferences=s.proposalPolicy;preferences.durationThousandths=300;juce::String error;expect(s.setProposalPolicy(preferences,error),"change preference while proposal is displayed");panel.refresh();lane.synchronize("proposal-lane");
    expect(panel.outcomeText().contains("Stale inputs")&&!proposalButton(panel,"Apply proposal").isEnabled()&&lane.inspectionText().contains("STALE")&&!s.applySectionProposal()&&s.currentState.analysis==baseline,"stale proposal retained for inspection with Apply rejected and baseline intact");
    proposalButton(panel,"Discard proposal").onClick();expect(!s.sectionProposal&&s.currentState.analysis==baseline,"main Discard preserves processing and measurements");
    expect(s.setProposalPolicy({},error),"restore explicit preferences");panel.refresh();proposalButton(panel,"Propose").onClick();
    const auto accepted=s.sectionProposal;proposalButton(panel,"Apply proposal").onClick();
    expect(s.lastAppliedExplanation&&s.lastAppliedExplanation->baseline.identity==baseline->identity&&s.lastAppliedExplanation->candidate&&s.lastAppliedExplanation->candidate->identity==accepted->assessment()->identity,"Apply snapshots baseline and private candidate scalars before invalidation");
    expect(s.currentStateNeedsRemeasurement()&&!s.currentState.hasMeasurement&&panel.outcomeText().contains("unverified")&&proposalButton(panel,"Last applied explanation").isEnabled(),"Apply leaves Current unverified and history directly inspectable");
    proposalButton(panel,"Last applied explanation").onClick();expect(historyOpened,"main saved explanation callback accessible");historicalRoundTrip(s,root.getChildFile("section-applied.levelscope"));
    const auto history=s.lastAppliedExplanation;measuredComposition(s,MeasurementDomain::currentState);auto noop=proposeFor(s);panel.refresh();
    expect(noop->status()==SectionProposalStatus::noChange&&panel.outcomeText().contains("Already meets requirements")&&!proposalButton(panel,"Apply proposal").isEnabled()&&s.applySectionProposal()&&s.lastAppliedExplanation==history,"already passes no-op preserves historical explanation and independent measurement");
    s.discardSectionProposal();expect(s.beginSectionTrimDraft(),"manual draft blocks new proposal explicitly");panel.refresh();expect(!proposalButton(panel,"Propose").isEnabled()&&panel.actionText().contains("manual gain draft"),"disabled Propose names unresolved manual draft");s.discardSectionTrimDraft();
    for(auto code:{ErrorCode::cancelled,ErrorCode::staleRevision,ErrorCode::decodeFailed})
    {
        s.sectionProposalError=code;s.sectionProposalNotice="Injected worker diagnostic";panel.refresh();
        expect(panel.outcomeText().contains(code==ErrorCode::cancelled?"Proposal cancelled":code==ErrorCode::staleRevision?"Stale inputs":"Proposal failed")&&!proposalButton(panel,"Apply proposal").isEnabled(),"cancellation stale completion and failure are distinct visible outcomes");
    }
    s.sectionProposalError=ErrorCode::none;expect(s.job.begin(JobStage::proposal,s.revision),"begin progress UI fixture");panel.refresh(.42f);
    expect(panel.outcomeText().contains("42%")&&!proposalButton(panel,"Apply proposal").isEnabled(),"proposal progress visible without private result authority");s.job.requestCancel();panel.refresh(.42f);expect(panel.outcomeText().contains("Cancelling"),"cancelling progress explicit");s.job.finish(s.job.ticket,s.revision,ErrorCode::cancelled,"cancelled");

    auto global=proposalFixture(root.getChildFile("global"),48000,2,0,.04f);auto g=proposeFor(global);
    longWindowRefresh(global);
    expect(g->status()==SectionProposalStatus::globalOnly&&g->passes()==0&&!g->assessment()&&g->previewPlan()->staticGainDb==*g->global().totalGainDb,"global-first preview is prediction only with no unnecessary candidate passes");
    expect(global.applySectionProposal()&&global.lastAppliedExplanation&&!global.lastAppliedExplanation->candidate&&global.lastAppliedExplanation->global.predictedIntegrated,"global-only Apply records predictions without inventing measured candidate");historicalRoundTrip(global,root.getChildFile("global-applied.levelscope"));
    auto direct=proposalFixture(root.getChildFile("direct"),48000,2,0,.04f);const auto trim=proposal::draftGlobalTrim(direct.currentTrimRequest());
    expect(direct.acceptTrimProposal(trim)&&direct.applyTrim()&&direct.lastAppliedExplanation&&direct.lastAppliedExplanation->method==proposal::trimMethodVersion,"legacy direct global Apply also captures complete explanation");historicalRoundTrip(direct,root.getChildFile("direct-applied.levelscope"));

    auto blocked=proposalFixture(root.getChildFile("blocked"),48000,2,0);auto request=blocked.sectionProposalRequest();
    auto map=mapFor(48000,480000,4);map.inputIdentity=request.annotations.inputIdentity;map.sections={{1,0,120000,"Long A"},{2,120000,240000,"Long B"},{3,240000,363024,"Long C"},{4,363024,480000,"Remaining"}};
    expect(blocked.sectionDocument.restore(map,map.inputIdentity),"synthetic long-section preference fixture");request=blocked.sectionProposalRequest();auto peaks=peakMap(map,-.5,3);peaks.channels=2;peaks.snapshotIdentity=request.scalar().current.snapshotIdentity;
    for(size_t n=0;n<peaks.sections.size();++n){auto& row=peaks.sections[n];row.truePeak.resize(2);row.sample.resize(2);if(n<3)row.truePeak[0]={std::pow(10.,-.5/20.),map.sections[n].begin*8,false};}request.attribution=std::make_shared<const proposal::SectionPeakMap>(peaks);JobControl limit;auto limited=generateSectionProposal(request,limit);
    expect(limited.value&&(*limited.value)->status()==SectionProposalStatus::noSolution&&(*limited.value)->issue()==SectionProposalIssue::preferences&&(*limited.value)->passes()==0&&!(*limited.value)->assessment(),"duration preference rejection does not claim candidate measurement");
    expect(blocked.acceptSectionProposal(*limited.value),"adopt rejected selection for inspection only");ProposalWorkflowPanel limitsPanel(blocked);limitsPanel.setSize(1200,ProposalWorkflowPanel::height);limitsPanel.refresh();
    std::cout << limitsPanel.outcomeText() << "\n";
    expect(limitsPanel.outcomeText().contains("76.13%")&&limitsPanel.outcomeText().contains("20.0%")&&limitsPanel.outcomeText().contains("Long A")&&limitsPanel.outcomeText().contains("Long B")&&limitsPanel.outcomeText().contains("Long C")&&!proposalButton(limitsPanel,"Apply proposal").isEnabled(),"76.13 percent rejection leads with 20 percent allowance and three long sections");
    TimelineGainLane rejectedLane(blocked);rejectedLane.setSize(1200,TimelineGainLane::height);rejectedLane.setViewport(0,480000,60,1130);rejectedLane.synchronize("rejected");rejectedLane.selectFrame(0,false);
    expect(rejectedLane.inspectionText().contains("REJECTED")&&rejectedLane.inspectionText().contains("baseline global"),"rejected local selection clearly uses baseline global, not final combined solution");
    for(auto begin:{120000,240000}){proposalButton(rejectedLane,">").onClick();expect(rejectedLane.inspectionText().contains("["+juce::String(begin)+","),"next change reaches every proposed section");}
    proposalImage(limitsPanel,root.getChildFile("duration-rejection.png"));
    auto requirements=blocked.jobRequirements;requirements.maximumLraTenths=180;expect(blocked.setJobRequirements(requirements,error),"set explicit 18 LU LRA requirement");auto unsupported=proposeFor(blocked);limitsPanel.refresh();
    expect(unsupported->passes()==0&&!unsupported->assessment()&&limitsPanel.outcomeText().contains("Unsupported correction")&&limitsPanel.outcomeText().contains("18.0 LU")&&limitsPanel.outcomeText().contains("No candidate was measured"),"LRA capability rejection names 18 LU requirement without claiming a measured candidate failed it");
    proposalImage(limitsPanel,root.getChildFile("unsupported-lra.png"));

    auto rejected=proposalFixture(root.getChildFile("measured-rejection"),48000,6,0,2.f);const auto rejectedResult=proposeFor(rejected);ProposalWorkflowPanel rejectedPanel(rejected);rejectedPanel.setSize(1200,ProposalWorkflowPanel::height);
    expect(rejectedResult->assessment()&&rejectedResult->issue()==SectionProposalIssue::candidateRejected&&rejectedPanel.outcomeText().contains("Measured private")&&rejectedPanel.outcomeText().contains("against")&&!proposalButton(rejectedPanel,"Apply proposal").isEnabled(),"measured candidate rejection names actual IL and TP values and selected limits without claiming LRA failure");
    expect(rejectedResult->text().contains("Current baseline before proposal")&&rejectedResult->text().contains("Measured private candidate"),"full review clearly separates baseline and candidate measurements");
    proposalImage(rejectedPanel,root.getChildFile("measured-candidate-rejection.png"));
    auto keyMap=mapFor(48000,480000,10000);std::ostringstream keyOracle;keyOracle.imbue(std::locale::classic());keyOracle << keyMap.inputIdentity.toStdString() << ':' << keyMap.sampleRate << ':' << keyMap.frames;
    for(const auto& row:keyMap.sections)keyOracle << ';' << row.id << ':' << row.begin << ':' << row.end;
    expect(proposal::sectionGeometryKey(keyMap).toStdString()==keyOracle.str()&&StandaloneSessionModel::sectionGeometry(keyMap).toStdString()==keyOracle.str(),"linear section geometry encoding preserves every established identity byte");
    auto denseRequest=rejected.sectionProposalRequest();auto denseMap=mapFor(48000,480000,10000);denseMap.inputIdentity=denseRequest.annotations.inputIdentity;
    expect(rejected.sectionDocument.restore(denseMap,denseMap.inputIdentity),"dense proposal uses 10000 original annotations");denseRequest=rejected.sectionProposalRequest();
    auto densePeaks=peakMap(denseMap,1.2,10000);densePeaks.snapshotIdentity=denseRequest.scalar().current.snapshotIdentity;denseRequest.attribution=std::make_shared<const proposal::SectionPeakMap>(densePeaks);
    JobControl denseJob;auto dense=generateSectionProposal(denseRequest,denseJob);expect(dense.value&&rejected.acceptSectionProposal(*dense.value)&&(*dense.value)->selection().corrections.size()==10000,"dense rejected proposal retains every change for review");
    TimelineGainLane denseLane(rejected);denseLane.setSize(1200,TimelineGainLane::height);denseLane.setViewport(0,480000,60,1130);denseLane.synchronize("dense-proposal");denseLane.selectFrame(479999,false);
    expect(denseLane.inspectionText().contains("Section 10000")&&denseLane.inspectionText().contains("REJECTED"),"last automatic correction remains inspectable in the lane");
    std::vector<double> paints;
    for(int n=0;n<12;++n){const auto start=juce::Time::getMillisecondCounterHiRes();rejectedPanel.refresh();denseLane.synchronize("dense-proposal");auto rendered=denseLane.createComponentSnapshot(denseLane.getLocalBounds());paints.push_back(juce::Time::getMillisecondCounterHiRes()-start);}
    std::sort(paints.begin(),paints.end());std::cout << "Dense proposal lane p95 ms " << paints[10] << "\n";expect(paints[10]<100,"dense automatic proposal controls and gain-lane paint p95 below 100 ms");
    auto densePerf=juce::var(new juce::DynamicObject());densePerf.getDynamicObject()->setProperty("dense_proposal_controls_lane_p95_ms",paints[10]);cases.add(densePerf);
    auto large=*history;large.annotations=mapFor(48000,4800000,10000);large.annotations.inputIdentity=large.inputIdentity;large.frames=large.annotations.frames;large.before.sectionGains.reset();large.before.schemaVersion=1;large.before.engineVersion="static-gain-v1";large.beforeIdentity=proposalPlanIdentity(large.before);
    for(auto& a:large.annotations.sections)a.name=juce::String::repeatedString(juce::String::fromUTF8("🙂"),128);
    auto scalarRequest=scalar(-24,1.2);auto largePeaks=peakMap(large.annotations,1.2,10000);auto policy=large.preferences;policy.durationThousandths=1000;
    large.preferences=policy;large.selection=proposal::selectSectionTrims(scalarRequest,large.annotations,{},largePeaks,policy);large.after.sectionGains=large.selection.plan.canonical();large.after.schemaVersion=2;large.after.engineVersion=processing::SectionGainPlan::version;large.afterIdentity=proposalPlanIdentity(large.after);
    const auto started=juce::Time::getMillisecondCounterHiRes();auto data=encodeAppliedExplanation(large);JobControl largeIo;auto decoded=decodeAppliedExplanation(data,largeIo);
    expect(data.size()<=AppliedProposalExplanation::maximumBytes&&decoded->selection.corrections.size()==10000&&decoded->annotations.sections.back().name==large.annotations.sections.back().name,"10000 complete corrections with maximum UTF-8 names fit 8 MiB history budget");
    auto all=decoded->text();ProposalInspectionPanel complete(all);
    expect(complete.pageCount()>1&&complete.allText()==all&&all.contains("Section 10000"),"all supported historical correction rows remain accessible without truncation");
    auto perf=juce::var(new juce::DynamicObject());perf.getDynamicObject()->setProperty("history_bytes_10000_max_names",juce::int64(data.size()));perf.getDynamicObject()->setProperty("history_roundtrip_and_inspection_ms",juce::Time::getMillisecondCounterHiRes()-started);cases.add(perf);
    for(auto mode:{0,1,2})
    {
        bool rejected=false;auto invalid=data;if(mode==0)invalid.resize(20);else if(mode==1)invalid.push_back(0);else invalid[0]=2;
        try{decodeAppliedExplanation(invalid,largeIo);}catch(...){rejected=true;}expect(rejected,"truncated trailing or future historical payload safely rejected");
    }
}
}
