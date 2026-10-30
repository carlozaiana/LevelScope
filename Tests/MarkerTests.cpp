#include "CompositionTestFixtures.h"
#include "Standalone/StandaloneMarkerImport.h"
#include "Standalone/StandaloneMarkerImportPanel.h"
#include <chrono>
#include <sys/resource.h>
namespace mk=levelscope::proposal::markers;
namespace
{
void requireReject(const std::function<void()>& fn,const char* name){bool rejected=false;try{fn();}catch(const std::exception&){rejected=true;}expect(rejected,name);}
proposal::SectionMap mapFor(int rate=48000,int64_t frames=269875200)
{proposal::SectionDocument d;d.bind(juce::String::repeatedString("a",64),rate,frames);return d.map();}
mk::Options pasted(mk::Unit unit=mk::Unit::videoFrames){mk::Options o;o.format=mk::Format::pastedPositions;o.unit=unit;o.confirmedPointMarkers=o.confirmedUnits=true;return o;}
void actualNuendo(const juce::File& file)
{
    if(file==juce::File())return;
    juce::MemoryBlock block;
    expect(file.loadFileAsData(block),"actual Nuendo 15.0.30 fixture read");
    const std::string bytes(static_cast<const char*>(block.getData()),block.getSize());
    expect(mk::contentHash(bytes)=="24941511d781d103f801a348b022001f3eda2b33887d694b37a87393f1b835b7","actual Nuendo fixture exact content identity");
    auto o=pasted(mk::Unit::nonDropTimecode);
    o.format=mk::Format::nuendoCsv;o.positionColumn=1;o.nameColumn=3;o.origin="00:00:00:00";
    auto p=mk::prepare(bytes,file.getFileName(),o,mapFor(),"nuendo-fixture");
    expect(p.receipt.rows.size()==6&&p.receipt.imported.sections.size()==7,"Nuendo actual point CSV produces six points without treating metadata as markers");
    const int64_t rawFrames[]={1209600,1495680,1633920,1660800,28823040,30205440};
    for(int i=0;i<6;++i)
        expect(p.receipt.rows[size_t(i)].frame==rawFrames[i],"raw exported labels under explicit zero origin; not claimed as DAW programme positions");

    // User reports 25 fps, project start two seconds before programme zero, and
    // first two visible markers 23:05 / 29:04 versus CSV 25:05 / 31:04.
    // Test an explicit export-coordinate origin, never an automatic Nuendo fix.
    o.origin="00:00:02:00";
    const auto mapped=mk::prepare(bytes,file.getFileName(),o,mapFor(),"nuendo-fixture");
    const int64_t programmeFrames[]={1113600,1399680,1537920,1564800,28727040,30109440};
    for(int i=0;i<6;++i)
    {
        expect(mapped.receipt.rows[size_t(i)].frame==programmeFrames[i]
            &&mapped.receipt.rows[size_t(i)].quantization=="0/1","explicit 50-frame origin maps exactly at 48 kHz");
        expect(mapped.receipt.rows[size_t(i)].position==p.receipt.rows[size_t(i)].position,"origin mapping preserves original exported labels");
    }
    expect(mapped.receipt.rows[0].seconds=="116/5"&&mapped.receipt.rows[1].seconds=="729/25","explicit origin agrees with first two user-reported Nuendo ruler positions");
    expect(mapped.receipt.identity()!=p.receipt.identity(),"origin participates in marker provenance identity");
    expect(mk::decode(mk::encode(mapped.receipt)).identity()==mapped.receipt.identity(),"nonzero Nuendo origin survives receipt validation and round trip");

    // A second untouched export from the same project adds a named marker at
    // displayed programme zero. It exports at 02:00, confirming this mapping.
    const auto anchorFile=file.getSiblingFile("nuendo-15.0.30-programme-zero.csv");
    juce::MemoryBlock anchorBlock;
    expect(anchorFile.loadFileAsData(anchorBlock),"actual Nuendo programme-zero anchor fixture read");
    const std::string anchorBytes(static_cast<const char*>(anchorBlock.getData()),anchorBlock.getSize());
    expect(mk::contentHash(anchorBytes)=="f4709037d4cbe30436be995b4fa3cf5cd9ce1795575ee97594fca60216e09cfd","programme-zero export exact content identity");
    const auto anchored=mk::prepare(anchorBytes,anchorFile.getFileName(),o,mapFor(),"nuendo-fixture");
    expect(anchored.receipt.rows.size()==7&&anchored.receipt.imported.sections.size()==7,"programme-zero marker does not create an empty section");
    expect(anchored.receipt.rows[0].position=="00:00:02:00"&&anchored.receipt.rows[0].frame==0
        &&anchored.receipt.imported.sections[0].name=="AudioVideoStart at 00:00:00","actual named programme-zero anchor maps to sample zero");
    for(int i=0;i<6;++i)
        expect(anchored.receipt.rows[size_t(i+1)].position==mapped.receipt.rows[size_t(i)].position
            &&anchored.receipt.rows[size_t(i+1)].frame==programmeFrames[i],"all original marker rows retain their verified mapping in the anchor export");
    expect(mk::decode(mk::encode(anchored.receipt)).identity()==anchored.receipt.identity(),"actual anchor name and origin survive receipt validation");
    cases.add("nuendo-actual-structure");
}
void actualProTools(const juce::File& file,const juce::File& original)
{
    if(file==juce::File())return;
    auto read=[](const juce::File& f){juce::MemoryBlock b;expect(f.loadFileAsData(b),"actual Pro Tools export read");return std::string(static_cast<const char*>(b.getData()),b.getSize());};
    const auto bytes=read(file);
    expect(mk::contentHash(bytes)=="f9f5bf49f5ddda6707e9031b587d78a383797fb84eb2af5b64b29c2ba8aaa460","Pro Tools 10.3.10 verbatim timing headers and marker table identity");
    auto o=pasted(mk::Unit::audioFrames);o.format=mk::Format::proToolsText;o.rateNumerator=24;
    auto m=mapFor(48000,48808914);const auto p=mk::prepare(bytes,file.getFileName(),o,m,"pro-tools-fixture");
    expect(p.receipt.rows.size()==45&&p.receipt.imported.sections.size()==46,"actual MARKERS LISTING at end imports all 45 points only");
    const int64_t samples[]={587888,1610400,3432560,3544000,4704912,5656912,6203840,6926912,8518000,8822000,9439590,9790000,10070000,13020000,18392000,23951856,24155856,24242528,24438000,25012912,26658000,28769904,29271152,29345152,29411152,29441152,29544304,29687152,29923152,30059984,30091152,30100000,30234000,30444000,30622912,31391856,31446128,31586000,32668000,36202752,37128304,38727792,38849792,39750000,41343984};
    for(size_t i=0;i<45;++i)expect(p.receipt.rows[i].frame==samples[i]&&p.receipt.rows[i].quantization=="0/1","TIME REFERENCE Samples preserved exactly without timecode-frame rounding");
    expect(p.receipt.rows.front().name.trim()=="Start Audio"&&p.receipt.rows.back().name.trim()=="END","actual fixed-width tab fields retain marker names");
    expect(mk::decode(mk::encode(p.receipt)).identity()==p.receipt.identity(),"actual sample mapping and provenance round trip");
    o.origin="587888";const auto shifted=mk::prepare(bytes,file.getFileName(),o,m,"pro-tools-fixture");
    expect(shifted.receipt.rows[0].frame==0&&shifted.receipt.rows[1].frame==1022512&&shifted.receipt.imported.sections.size()==45,"explicit sample-origin subtraction is exact; this test does not infer programme start from marker name");
    o.origin="0";
    requireReject([&]{mk::prepare(bytes,file.getFileName(),o,mapFor(44100,48808914),"r");},"different programme rate rejects sample-reference mapping without resampling");
    auto replace=[&](const juce::String& from,const juce::String& to){return juce::String::fromUTF8(bytes.data(),int(bytes.size())).replaceFirstOccurrenceOf(from,to).toStdString();};
    requireReject([&]{mk::prepare(replace("SAMPLE RATE:","UNUSED RATE:"),"x",o,m,"r");},"missing sample-rate header rejected");
    requireReject([&]{mk::prepare("SAMPLE RATE:\t48000\n"+bytes,"x",o,m,"r");},"ambiguous duplicate sample-rate headers rejected");
    requireReject([&]{mk::prepare(replace("Samples","Ticks"),"x",o,m,"r");},"mixed musical/sample reference units reject entire table");
    requireReject([&]{mk::prepare(replace("TIME REFERENCE","OTHER"),"x",o,m,"r");},"missing sample reference never falls back to LOCATION");
    requireReject([&]{mk::prepare(replace("587888","587888.5"),"x",o,m,"r");},"fractional sample references rejected");
    requireReject([&]{mk::prepare(replace("587888","9223372036854775808"),"x",o,m,"r");},"sample-reference overflow rejected");
    requireReject([&]{mk::prepare(bytes+"M A R K E R S\n","x",o,m,"r");},"both recognized marker headings count toward ambiguity rejection");
    auto tc=o;tc.unit=mk::Unit::nonDropTimecode;tc.origin="00:00:00:00";
    const auto timecode=mk::prepare(bytes,file.getFileName(),tc,m,"r");
    expect(timecode.receipt.rows[0].frame==586000&&timecode.receipt.rows[1].frame==1610000,"explicit 24 fps LOCATION maps display labels independently of sample references");
    for(size_t i=0;i<45;++i)expect(samples[i]>=timecode.receipt.rows[i].frame&&samples[i]<timecode.receipt.rows[i].frame+2000,"actual sample references fall within their displayed 24 fps frame; lost subframe precision is disclosed");
    expect(timecode.receipt.notices.joinIntoString("\n").contains("subframe"),"timecode-only precision limitation visible in preview/provenance");
    if(original!=juce::File())
    {
        const auto full=read(original);expect(mk::contentHash(full)=="10aecddf80f7dcb668315f68cfd4782a10fe0b218a3c2d487529267fbc60dc25","full user export exact identity, kept outside repository");
        const auto whole=mk::prepare(full,original.getFileName(),o,m,"r");
        expect(whole.receipt.rows.size()==45&&whole.receipt.rows.front().line==6347,"large unrelated track/file listings are excluded from actual full export");
        for(size_t i=0;i<45;++i)expect(whole.receipt.rows[i].frame==p.receipt.rows[i].frame&&whole.receipt.rows[i].name==p.receipt.rows[i].name,"full original and minimal verbatim fixture have identical selected marker content");
    }
    cases.add("pro-tools-actual-structure");
}
void mappingCases()
{
    auto o=pasted();auto m=mapFor();auto p=mk::prepare("42397,42575,45099,47283,130243,130547","pasted",o,m,"runtime");
    const int64_t expected[]={81402240,81744000,86590080,90783360,250066560,250650240};
    expect(p.receipt.rows.size()==6&&p.receipt.imported.sections.size()==7,"six cinema frame points produce seven sections without media I/O");
    for(int i=0;i<6;++i)expect(p.receipt.rows[size_t(i)].frame==expected[i]&&p.receipt.rows[size_t(i)].quantization=="0/1","independent cinema boundary equals exact 25fps audio sample");
    expect(p.receipt.imported.sections.back().end==269875200,"EOF unchanged");
    for(int rate:{44100,48000,96000,192000})
    {
        auto m0=mapFor(rate,int64_t(rate)*10);auto seconds=pasted(mk::Unit::seconds);seconds.origin="3600";
        auto q=mk::prepare("3600.1\n7203/2","paste",seconds,m0,"r");
        expect(q.receipt.rows[0].frame==rate/10&&q.receipt.rows[1].frame==rate*3/2,"exact decimal/rational seconds with nonzero origin");
        auto frames=pasted(mk::Unit::audioFrames);frames.origin="17";q=mk::prepare("18","paste",frames,m0,"r");expect(q.receipt.rows[0].frame==1,"audio frame origin subtraction");
        auto tc=pasted(mk::Unit::nonDropTimecode);tc.origin="01:00:00:00";tc.rateNumerator=24000;tc.rateDenominator=1001;
        q=mk::prepare("01:00:00:01","paste",tc,m0,"r");expect(q.receipt.rows[0].frame==int64_t(rate)*1001/24000,"fractional-rate NDF uses nominal 24 labels and rational 24000/1001");
        tc.rateNumerator=30000;q=mk::prepare("01:00:00:01","paste",tc,m0,"r");expect(q.receipt.rows[0].frame==int64_t(rate)*1001/30000,"fractional-rate 30-base NDF origin");
    }
    o=pasted(mk::Unit::seconds);auto q=mk::prepare("1/96000","p",o,m,"r");expect(q.receipt.rows[0].frame==0&&q.receipt.rows[0].quantization=="1/2","floor reports exact half-sample earlier adjustment");
    for(auto bad:{"nan","inf","-1","1/0","9223372036854775808","1e2","5622.400000001","1/2/3"})requireReject([&]{mk::prepare(bad,"p",o,m,"r");},"invalid seconds/overflow/position beyond EOF rejected before rounding");
    o=pasted(mk::Unit::nonDropTimecode);o.origin="01:00:00:00";for(auto bad:{"01:00:00;01","01:60:00:00","00:59:59:24","24:00:00:00","01:00:00:25"})requireReject([&]{mk::prepare(bad,"p",o,m,"r");},"drop-frame invalid field negative-relative and wrap rejected");
    o=pasted();o.rateNumerator=2997;o.rateDenominator=100;requireReject([&]{mk::prepare("1","p",o,m,"r");},"approximate or unsupported rate rejected");
    o=pasted();o.confirmedUnits=false;requireReject([&]{mk::prepare("1","p",o,m,"r");},"explicit actual-unit confirmation required");
    cases.add("exact-mapping");
}
void parserCases()
{
    auto o=pasted();o.format=mk::Format::nuendoCsv;auto m=mapFor();
    auto p=mk::prepare("\xef\xbb\xbfPosition,Name,Unused\r\n25,\"A, \"\"quote\"\"\",x\r\n50,é,x\r\n","test.csv",o,m,"r");
    expect(p.receipt.rows[0].name=="A, \"quote\""&&p.receipt.rows[1].frame==96000,"UTF8 BOM CRLF quoted delimiters escaped quotes and Unicode");expect(p.receipt.notices.size()==1,"unused CSV columns disclosed");
    auto spaced=mk::prepare("Position,Name\n 25 , padded name ","x",o,m,"r");expect(spaced.receipt.rows[0].position==" 25 "&&spaced.receipt.rows[0].name==" padded name ","original selected position/name text retained exactly");
    auto restored=mk::decode(mk::encode(p.receipt));expect(restored.identity()==p.receipt.identity(),"receipt exact serialization round trip");
    for(char separator:{';','\t'}){o.separator=separator;auto s=std::string("Position")+separator+"Name\n25"+separator+"A\n";expect(mk::prepare(s,"test.csv",o,m,"r").receipt.rows[0].frame==48000,"explicit alternate separator");}o.separator=',';
    for(auto bad:{"Position,Name\n25,A\n\n50,B","Position,Name\n25,\"A","Position,Name,Type\n25,A,Cycle","Position,Name,End\n25,A,50","Timecode In,Description,Timecode Out\n25,A,50","Position,Name,Duration\n25,A,1","Position,Name,Track\n25,A,one\n50,B,two","Position,Name\n25,A\n25,B"})requireReject([&]{mk::prepare(bad,"test.csv",o,m,"r");},"malformed rows and unsupported or conflicting markers block whole import");
    o.typeColumn=2;requireReject([&]{mk::prepare("Position,Name,Type\n25,A,Cycle","x",o,m,"r");},"explicit cycle marker rejected");
    p=mk::prepare("Position,Name,Type\n50,B,Point\n25,A,Point\n25,A,Point\n0,Start,Point","x",o,m,"r");
    expect(p.receipt.rows[2].coalesced&&p.receipt.imported.sections[0].name=="Start"&&p.receipt.imported.sections.size()==3,"visible sort duplicate coalescing and zero endpoint");
    o=pasted();auto bytes=std::string("\xff");requireReject([&]{mk::prepare(bytes,"x",o,m,"r");},"unsupported encoding rejected");
    bytes=std::string("25\0",3);requireReject([&]{mk::prepare(bytes,"x",o,m,"r");},"embedded NUL rejected");
    o.format=mk::Format::proToolsText;o.unit=mk::Unit::nonDropTimecode;o.origin="01:00:00:00";
    const std::string pt="SESSION NAME:\tSynthetic test\nTIME FORMAT:\tNot authoritative for LOCATION\n\nM A R K E R S\n#\tLOCATION\tTIME REFERENCE\tUNITS\tNAME\tCOMMENTS\n1\t01:00:01:00\tAbsolute\tTimecode\tIntro\tignored\n2\t01:00:02:00\tAbsolute\tTimecode\tNext\tignored\n";
    p=mk::prepare(pt,"session.txt",o,m,"r");expect(p.receipt.rows.size()==2&&p.receipt.rows[0].frame==48000&&p.receipt.rows[0].name=="Intro","dedicated Pro Tools marker table reader ignores session headers and keeps selected field");
    requireReject([&]{mk::prepare("Position,Name\n25,A","x",o,m,"r");},"generic CSV is not Pro Tools text");
    requireReject([&]{mk::prepare(pt+"unsupported structure","x",o,m,"r");},"unrecognized Pro Tools structure never silently discarded");
    requireReject([&]{mk::prepare(pt+"M A R K E R S\n","x",o,m,"r");},"ambiguous multiple Pro Tools sections rejected");
    auto v=mk::encode(p.receipt);v.getDynamicObject()->setProperty("inputSha256","bad");requireReject([&]{mk::decode(v);},"damaged receipt identity rejected");
    o=pasted();requireReject([&]{mk::prepare("25","x",o,m,"r",[]{return true;});},"cancelled parsing leaves no partial map");cases.add("bounded-readers");
}

StandaloneSessionModel smallModel(const juce::File& root,int channels)
{
    const int rate=48000,frames=rate*3;juce::AudioBuffer<float> raw(channels,frames);
    for(int c=0;c<channels;++c)for(int n=0;n<frames;++n)raw.setSample(c,n,float(.018*std::sin((n+c)*.031)));
    const auto f=root.getChildFile("source.wav");floatFile(f,rate,raw);JobControl j;auto source=validateMedia(f,j);expect(source.succeeded(),"short source validated");
    StandaloneSessionModel s;expect(s.setValidatedSource(f,source.media),"adopt short source");s.sourceDocument.layout=layoutFor(channels);s.sourceDocument.layout.confirmed=true;
    expect(s.initializeCurrentStateFromSource(),"initialize Current");s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);s.synchronizeSections();return s;
}
void transactionCases(const juce::File& root)
{
    for(int channels:{1,2,6})
    {
        auto directory=root.getChildFile("transaction-"+juce::String(channels));directory.createDirectory();auto s=smallModel(directory,channels);juce::String error;
        expect(s.sectionDocument.split(s.sectionDocument.map().sections[0].id,48000,error),"existing annotation boundary");s.documentChanged();
        expect(s.beginSectionTrimDraft()&&s.setSectionDraftGain(0,-1.23456789012345)&&s.applySectionTrimDraft(),"apply non-unity absolute section plan explicitly");
        measuredComposition(s,MeasurementDomain::source);Capture before(channels,144000);measuredComposition(s,MeasurementDomain::currentState,&before);
        auto draft=proposal::draftGlobalTrim(s.currentTrimRequest());expect(draft.status==proposal::TrimStatus::draft&&s.acceptTrimProposal(draft)&&s.applyTrim(),"explicit passing global trim above applied sections");measuredComposition(s,MeasurementDomain::currentState);
        Capture baseline(channels,144000);measuredComposition(s,MeasurementDomain::currentState,&baseline);
        expect(s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"measured Current passes before annotations");
        auto source=s.source.analysis,current=s.currentState.analysis;auto oldRequest=ex::capture(s);JobControl control;
        auto artifact=ex::render(oldRequest,directory,control);expect(artifact.value.has_value(),"short verified artifact before annotation change");const auto hash=(*artifact.value)->hash(),oldApproval=ex::approvalIdentity(**artifact.value,oldRequest);
        const auto plan=s.currentStateDocument.processingPlan();const auto gainKey=s.currentStateDocument.sectionGains->audioKey(),applied=s.currentStateDocument.sectionGainAnnotations.digest();const auto audioRevision=s.revision,documentRevision=s.documentRevision;const auto proposalBinding=s.sectionProposalRequest().binding();
        auto o=pasted(mk::Unit::audioFrames);const std::string bytes="24000,72000,120000";auto p=mk::prepare(bytes,"Pasted positions",o,s.sectionDocument.map(),s.sessionId);
        expect(s.beginSectionTrimDraft(),"pre-existing manual draft for stale-geometry check");expect(commitMarkerImport(s,p,mk::contentHash(bytes),error),"atomic annotation import while applied processing and draft exist");
        expect(s.revision==audioRevision&&s.documentRevision==documentRevision+1&&s.source.analysis==source&&s.currentState.analysis==current,"import preserves audio revision and measured snapshots; one document revision");
        expect(s.sectionProposalRequest().binding()!=proposalBinding,"section-led proposal binding becomes stale after annotation geometry change");
        expect(!s.sectionDraftFresh()&&s.sectionTrimDraft.has_value(),"geometry-dependent draft retained visibly stale");
        expect(s.currentStateDocument.sectionGains->audioKey()==gainKey&&s.currentStateDocument.sectionGainAnnotations.digest()==applied&&s.currentStateDocument.staticGainDb==plan.staticGainDb,"applied absolute section intervals and full-precision global gain unchanged");
        expect(s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"unchanged measured compliance remains current");s.discardSectionTrimDraft();
        auto newRequest=ex::capture(s);expect(ex::audioIdentity(oldRequest)==ex::audioIdentity(newRequest)&&ex::requestIdentity(oldRequest)!=ex::requestIdentity(newRequest),"report identity changes but rendered audio identity does not");
        expect((*artifact.value)->matchesAudio(newRequest)&&ex::approvalIdentity(**artifact.value,newRequest)!=oldApproval,"existing artifact reusable; obsolete report approval rejected");
        auto stalePublish=ex::publish(**artifact.value,newRequest,directory.getChildFile("stale-export"),oldApproval,control);expect(!stalePublish.value,"stale report publication authorization rejected");
        auto refreshed=ex::artifactReport(**artifact.value,newRequest,directory.getChildFile("refreshed-report"),control);expect(refreshed.value.has_value()&&(*artifact.value)->hash()==hash,"report refreshed without rerendering unchanged artifact");
        auto report=ex::report(newRequest,artifact.value->get());expect(report["marker_import_provenance"].isObject()&&bool(report["annotations_match_last_import"]),"report carries current receipt separately from measurements");
        auto hostile=newRequest;hostile.sections.sections[0].name="<script>alert('marker')</script> & \"quote\"";const auto html=ex::reportHtml(ex::report(hostile));expect(!html.contains("<script>")&&html.contains("&lt;script&gt;"),"marker/annotation labels are escaped report text, never executable HTML");
        const auto revision=s.documentRevision;auto again=mk::prepare(bytes,"Pasted positions",o,s.sectionDocument.map(),s.sessionId);expect(commitMarkerImport(s,again,mk::contentHash(bytes),error)&&s.documentRevision==revision,"identical map and receipt are genuine no-op");
        auto replaced=s;replaced.sourceDocument.media.sha256=juce::String::repeatedString("b",64);const auto replacedReport=ex::capture(replaced);expect(replacedReport.sections.sections.empty()&&!replacedReport.markerImport,"report never carries a previous programme annotation receipt after source replacement");
        auto wrong=again;wrong.runtime="old-runtime";expect(!commitMarkerImport(s,wrong,mk::contentHash(bytes),error),"stale runtime cannot commit");expect(!commitMarkerImport(s,again,"changed bytes",error),"changed marker bytes cannot commit");
        bool busy=true;s.auditionBusy=[&]{return busy;};expect(!commitMarkerImport(s,again,mk::contentHash(bytes),error),"active audition blocks confirmation without stopping it");busy=false;s.auditionBusy={};
        expect(s.job.begin(JobStage::analysis,s.revision),"start isolated busy job");expect(!commitMarkerImport(s,again,mk::contentHash(bytes),error),"active measurement job blocks confirmation");s.job=SessionJob{};
        auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));OutputRouting muted;
        expect(s.selectAuditionDomain(MeasurementDomain::currentState),"explicit Current audition selection");auto request=s.auditionRequest();audition.synchronize(request,muted);expect(audition.connect(48000).isEmpty(),"injected-only output opened");Capture played(channels,144000);audition.setCaptureObserver(&played);expect(audition.play(request,muted,true,true).isEmpty(),"explicit injected audition after annotation import");awaitPlaying(audition,*injected);drain(audition,*injected,997);
        expect(played.valid&&played.next==144000&&sameSamples(played.data,baseline.data)&&audition.underruns()==0,"all-channel pre-monitor capture including LFE matches unchanged offline processing at EOF");audition.setCaptureObserver(nullptr);const auto* verifiedCopy=audition.preparedCopy();expect(audition.seek(47003,request,muted).isEmpty(),"annotation import does not change audio-only seek coordinates");Capture sought(channels,144000-47003,47003);audition.setCaptureObserver(&sought);expect(audition.play(request,muted,false,true).isEmpty(),"explicit restart from audio coordinate");awaitPlaying(audition,*injected);drain(audition,*injected,127);bool same=sought.valid&&sought.next==144000;for(int ch=0;ch<channels;++ch)same&=std::memcmp(sought.data.getReadPointer(ch),baseline.data.getReadPointer(ch,47003),size_t(144000-47003)*4)==0;expect(same&&audition.preparedCopy()==verifiedCopy,"seek samples and source-copy reuse invariant");audition.setCaptureObserver(nullptr);audition.disconnect();
        auto doc=sf::capture(s);const auto file=directory.getChildFile("markers.levelscope");expect(sf::save(doc,file,control).succeeded(),"current schema coherent session saved");auto read=sf::read(file,control);expect(read.document&&read.document->loadedSchema==6&&read.document->markerImport&&read.document->markerImport->identity()==s.markerImport->identity(),"current schema exact receipt recovered without marker file");
        auto verified=sf::verify(*read.document,control);StandaloneSessionModel restored;expect(verified.document&&verified.prepared&&sf::adopt(*verified.document,*verified.prepared,restored,error),"current schema media verification and atomic adoption");expect(restored.source.hasMeasurement&&restored.currentState.hasMeasurement&&restored.sourceDocument.routing.fileToDevice[0]==-1,"same-build caches reused and routes muted on reopen");
        expect(restored.sessionId!=s.sessionId&&!restored.acceptMeasurement(*current),"late measurement cannot affect reopened runtime");
        const auto coherent=sf::capture(s);expect(s.sectionDocument.rename(s.sectionDocument.map().sections[0].id,"Changed after save capture",error),"later annotation edit");s.documentChanged();expect(sf::saved(s,coherent,file)&&s.documentDirty(),"later edits remain unsaved after coherent snapshot completion");expect(!commitMarkerImport(s,again,mk::contentHash(bytes),error),"competing annotation edit rejects old preview");
        Capture after(channels,144000);measuredComposition(s,MeasurementDomain::currentState,&after);expect(sameSamples(baseline.data,after.data),"all-channel processed samples bit-identical after marker import and name edit");
        sf::FileHooks failure;failure.fail=[](const char* stage){return juce::String(stage)=="replace";};auto original=juce::SHA256(file).toHexString();expect(!sf::save(sf::capture(s),file,control,failure).succeeded()&&juce::SHA256(file).toHexString()==original,"failed session replacement retains original complete document");
        cases.add("transaction/session/report "+juce::String(channels));
    }
}

void rewrittenCore(const juce::File& original,const juce::File& destination,const std::function<void(juce::var&)>& edit)
{
    juce::MemoryBlock block;expect(original.loadFileAsData(block),"load owned legacy-shape fixture");const auto* data=static_cast<const uint8_t*>(block.getData());uint64_t length=0;for(int i=0;i<8;++i)length|=uint64_t(data[8+i])<<(8*i);
    auto value=juce::JSON::parse(juce::String::fromUTF8(reinterpret_cast<const char*>(data+80),int(length)));edit(value);const auto json=juce::JSON::toString(value,true).toStdString();auto out=destination.createOutputStream();expect(out!=nullptr,"write owned session-shape fixture");out->write("LSSESS01",8);out->writeInt64(int64_t(json.size()));const auto hash=mk::contentHash(json).toStdString();out->write(hash.data(),hash.size());out->write(json.data(),json.size());out->write(data+80+length,block.getSize()-80-size_t(length));out->flush();expect(out->getStatus().wasOk(),"session-shape fixture finalized");
}
void legacyCases(const juce::File& root)
{
    root.createDirectory();auto s=smallModel(root,2);measuredComposition(s,MeasurementDomain::source);measuredComposition(s,MeasurementDomain::currentState);JobControl job;const auto original=root.getChildFile("original-schema6.levelscope");expect(sf::save(sf::capture(s),original,job).succeeded(),"current schema baseline for legacy core-shape tests");const auto hash=juce::SHA256(original).toHexString();
    for(int schema:{1,2,3,4})
    {
        auto file=root.getChildFile("schema-"+juce::String(schema)+".levelscope");rewrittenCore(original,file,[schema](juce::var& v){v.getDynamicObject()->setProperty("schema",schema);v["programme"].getDynamicObject()->setProperty("schema",schema);v.getDynamicObject()->removeProperty("markerImport");});const auto oldHash=juce::SHA256(file).toHexString();auto read=sf::read(file,job);expect(read.document&&read.document->loadedSchema==schema&&!read.document->markerImport,"legacy schema 1-4 recovers settings without invented import provenance");expect(read.document->notices.joinIntoString(" ").contains("M16 reader"),"older session upgrade and Save As warning explicit");
        auto verified=sf::verify(*read.document,job);StandaloneSessionModel opened;juce::String error;expect(verified.document&&verified.prepared&&sf::adopt(*verified.document,*verified.prepared,opened,error),"legacy core settings install without rewriting original");expect(juce::SHA256(file).toHexString()==oldHash,"opening leaves legacy file unchanged");auto copy=root.getChildFile("schema-"+juce::String(schema)+"-SaveAs.levelscope");expect(sf::save(sf::capture(opened),copy,job).succeeded()&&juce::SHA256(file).toHexString()==oldHash,"Save As writes new format and preserves legacy readable original");
    }
    expect(juce::SHA256(original).toHexString()==hash,"baseline original unchanged");
    auto damaged=root.getChildFile("damaged-cache.levelscope");juce::MemoryBlock b;original.loadFileAsData(b);static_cast<uint8_t*>(b.getData())[b.getSize()-1]^=1;damaged.replaceWithData(b.getData(),b.getSize());auto read=sf::read(damaged,job);expect(read.document&&!read.document->currentCache&&read.document->sourceCache&&!read.document->notices.isEmpty(),"replaceable damaged Current cache discarded; source/settings retained");
    cases.add("legacy schemas / cache recovery");
    auto sevenRoot=root.getChildFile("seven");sevenRoot.createDirectory();auto seven=smallModel(sevenRoot,8);auto o=pasted(mk::Unit::seconds);auto p=mk::prepare("1\n2","7.1 markers",o,seven.sectionDocument.map(),seven.sessionId);const auto revision=seven.revision;juce::String error;expect(commitMarkerImport(seven,p,mk::contentHash("1\n2"),error)&&seven.revision==revision&&seven.sectionDocument.map().sections.size()==3,"7.1 annotations supported without enabling unsupported section processing/export");expect(!seven.beginSectionTrimDraft(),"7.1 processing boundary preserved");
}

juce::Component* named(juce::Component& root,const juce::String& name)
{if(root.getName()==name)return &root;for(auto* child:root.getChildren())if(auto* result=named(*child,name))return result;return nullptr;}
void uiBudgetCase(const juce::File& root)
{
    root.createDirectory();auto s=smallModel(root,2);MarkerImportWindow window(s,true,[]{return true;},[]{});
    auto* units=dynamic_cast<juce::ComboBox*>(named(window,"Marker position units"));auto* text=dynamic_cast<juce::TextEditor*>(named(window,"Pasted marker positions"));auto* status=dynamic_cast<juce::Label*>(named(window,"Marker import status"));
    auto* points=dynamic_cast<juce::ToggleButton*>(named(window,"I confirm this export contains point markers only"));auto* confirmed=dynamic_cast<juce::ToggleButton*>(named(window,"I confirm the actual position-field units, rate and origin; no inference from filenames or general DAW settings"));auto* preview=dynamic_cast<juce::TextButton*>(named(window,"Preview mapping"));auto* commit=dynamic_cast<juce::TextButton*>(named(window,"Confirm: replace annotations"));
    expect(units&&text&&status&&points&&confirmed&&preview&&commit,"real marker panel controls found for resource heartbeat harness");units->setSelectedId(1);std::string input;for(int i=0;i<10000;++i)input+=std::to_string(i*10)+"\n";text->setText(input);points->setToggleState(true,juce::sendNotification);confirmed->setToggleState(true,juce::sendNotification);
    struct Heartbeat:juce::Timer{std::vector<double> samples;double last=juce::Time::getMillisecondCounterHiRes();void timerCallback()override{auto now=juce::Time::getMillisecondCounterHiRes();samples.push_back(now-last);last=now;}} heartbeat;heartbeat.startTimer(5);
    auto pump=[&](const std::function<bool()>& done){const auto deadline=juce::Time::getMillisecondCounterHiRes()+5000;while(!done()&&juce::Time::getMillisecondCounterHiRes()<deadline)juce::MessageManager::getInstance()->runDispatchLoopUntil(5);expect(done(),"bounded panel worker completes within five seconds");};
    juce::MessageManager::getInstance()->runDispatchLoopUntil(30);preview->triggerClick();pump([&]{return status->getText().startsWith("Preview only:");});
    expect(s.sectionDocument.map().sections.size()==1,"10000-row preview has no annotation side effects");commit->triggerClick();pump([&]{return status->getText().startsWith("Annotations replaced.");});expect(s.sectionDocument.map().sections.size()==10000,"real panel confirms complete annotation map atomically");
    juce::MessageManager::getInstance()->runDispatchLoopUntil(100);heartbeat.stopTimer();expect(!heartbeat.samples.empty(),"UI heartbeat recorded");auto samples=heartbeat.samples;std::sort(samples.begin(),samples.end());const auto p95=samples[size_t(.95*double(samples.size()-1))],maximum=samples.back();expect(p95<=100.,"real panel 10000-row preview/commit heartbeat p95 within 100 ms");
    auto evidence=juce::var(new juce::DynamicObject());evidence.getDynamicObject()->setProperty("group","ui-resource");evidence.getDynamicObject()->setProperty("ui_heartbeat_p95_ms",p95);evidence.getDynamicObject()->setProperty("ui_heartbeat_max_ms",maximum);evidence.getDynamicObject()->setProperty("ticks",int(samples.size()));evidence.getDynamicObject()->setProperty("rows",10000);cases.add(evidence);window.closeButtonPressed();
}
void resourceCases()
{
    auto o=pasted(mk::Unit::audioFrames);auto m=mapFor(48000,48000LL*10800);std::string text;for(int i=0;i<10000;++i)text+=std::to_string(int64_t(i)*48000)+"\n";
    struct rusage before{},after{};getrusage(RUSAGE_SELF,&before);auto start=std::chrono::steady_clock::now();auto p=mk::prepare(text,"resource",o,m,"r");auto duration=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();getrusage(RUSAGE_SELF,&after);
    expect(p.receipt.rows.size()==10000&&p.receipt.imported.sections.size()==10000,"maximum representative rows produce full map without truncation");expect(duration<=2.,"maximum representative parse/map within two seconds");expect(after.ru_maxrss-before.ru_maxrss<=64LL*1024*1024,"incremental importer peak RSS within 64 MiB");
    auto evidence=juce::var(new juce::DynamicObject());evidence.getDynamicObject()->setProperty("group","resource-limits");evidence.getDynamicObject()->setProperty("parse_map_seconds",duration);evidence.getDynamicObject()->setProperty("incremental_peak_rss_bytes",double(after.ru_maxrss-before.ru_maxrss));evidence.getDynamicObject()->setProperty("rows",10000);evidence.getDynamicObject()->setProperty("input_bytes",int(text.size()));evidence.getDynamicObject()->setProperty("three_hour_audio_executed",false);cases.add(evidence);
    auto csvOptions=o;csvOptions.format=mk::Format::nuendoCsv;std::string large="Position,Name,Unused\n";for(int n=0;n<10000;++n)large+=std::to_string(int64_t(n)*48000)+",M,"+std::string(810,'x')+"\n";
    const auto padding=mk::maximumBytes-large.size();large.insert(large.size()-1,padding,'x');expect(large.size()==mk::maximumBytes,"representative maximum 8 MiB input constructed without media");
    start=std::chrono::steady_clock::now();auto largePreview=mk::prepare(large,"large.csv",csvOptions,m,"r");duration=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();getrusage(RUSAGE_SELF,&after);
    expect(largePreview.receipt.rows.size()==10000&&duration<=2.,"8 MiB 10000-row CSV parse/map within two seconds");expect(after.ru_maxrss-before.ru_maxrss<=64LL*1024*1024,"maximum byte and row workload added peak RSS within 64 MiB");
    evidence.getDynamicObject()->setProperty("eight_mib_seconds",duration);evidence.getDynamicObject()->setProperty("eight_mib_peak_rss_delta",double(after.ru_maxrss-before.ru_maxrss));
    requireReject([&]{mk::prepare(text+"480000001\n","resource",o,m,"r");},"more than 10000 input markers rejected");
    requireReject([&]{mk::prepare(std::string(mk::maximumBytes+1,'1'),"resource",o,m,"r");},"input above 8 MiB rejected");
    int polls=0;start=std::chrono::steady_clock::now();requireReject([&]{mk::prepare(text,"resource",o,m,"r",[&]{return ++polls>10;});},"in-progress parser cancellation rejects all annotations");duration=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();expect(duration<=1.,"cancellation completes within one second");
}
}
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI init;juce::File root(argc>1?argv[1]:"/tmp/levelscope-marker-tests");root.createDirectory();bool passed=true;juce::String error;
    auto begin=juce::Time::getMillisecondCounterHiRes();try{mappingCases();parserCases();actualNuendo(argc>2?juce::File(argv[2]):juce::File());actualProTools(argc>2?juce::File(argv[2]).getSiblingFile("pro-tools-10.3.10-markers.txt"):juce::File(),argc>3?juce::File(argv[3]):juce::File());resourceCases();transactionCases(root);legacyCases(root.getChildFile("legacy"));uiBudgetCase(root.getChildFile("ui-resource"));}catch(const std::exception& e){passed=false;error=e.what();std::cerr<<error<<'\n';}
    auto report=juce::var(new juce::DynamicObject());auto* r=report.getDynamicObject();r->setProperty("passed",passed);r->setProperty("error",error);r->setProperty("checks",checks);r->setProperty("cases",cases);juce::Array<juce::var> completed;for(const auto& item:cases)completed.add(item.isString()?item:item["group"]);r->setProperty("completed_groups",completed);r->setProperty("source_build",LEVELSCOPE_MEASUREMENT_BUILD);r->setProperty("runtime",juce::SystemStats::getOperatingSystemName());r->setProperty("wall_ms",juce::Time::getMillisecondCounterHiRes()-begin);r->setProperty("actual_daw_interoperability",argc>2?"Nuendo 15.0.30 actual CSV at user-confirmed 25 fps: raw labels and explicit two-second origin tested; first two mapped positions agree with user observations. Second actual export confirms the named programme-zero anchor; support is limited to this tested export setup. Pro Tools 10.3.10 actual table: explicit 24 fps LOCATION and matching-rate TIME REFERENCE / Samples tested; programme-origin and visible-position confirmation recorded separately.":"Development run without actual export fixture; synthetic reader checks only.");root.getChildFile("marker-tests.json").replaceWithText(juce::JSON::toString(report,true));std::cout<<checks.size()<<" checks, "<<cases.size()<<" groups, "<<(passed?"PASS":"FAIL")<<'\n';return passed?0:1;
}
