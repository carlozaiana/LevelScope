#include "StandaloneExport.h"
#include "StandaloneProcessing.h"
#include <juce_cryptography/juce_cryptography.h>
#include <sstream>
#include <iomanip>
#include <locale>
#include <stdexcept>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mount.h>
#include <stdio.h>

namespace levelscope::standalone::exporting
{
namespace
{
struct Failure : std::runtime_error
{ Failure (ErrorCode c, const juce::String& s) : std::runtime_error (s.toStdString()), code (c) {} ErrorCode code; };
void demand (bool ok, ErrorCode code, const juce::String& message) { if (!ok) throw Failure (code, message); }
juce::String digest (const juce::String& s) { return juce::SHA256 (s.toRawUTF8(), s.getNumBytesAsUTF8()).toHexString(); }
juce::String exact (double v) { std::ostringstream s; s.imbue (std::locale::classic()); s << std::hexfloat << v; return s.str(); }
bool sameFolder (const juce::File& f, uint64_t device, uint64_t inode)
{ struct stat s {}; return ::lstat (f.getFullPathName().toRawUTF8(), &s) == 0 && S_ISDIR (s.st_mode)
      && uint64_t (s.st_dev) == device && uint64_t (s.st_ino) == inode; }
void removeOwnedFolder (const juce::File& f, uint64_t device, uint64_t inode)
{
    if (!sameFolder (f, device, inode)) return;
    for (const auto& child : f.findChildFiles (juce::File::findFilesAndDirectories, false))
        ::unlink (child.getFullPathName().toRawUTF8()); // Never follow links or recursively delete foreign directories.
    ::rmdir (f.getFullPathName().toRawUTF8());
}
struct Folder
{
    juce::File file; uint64_t device = 0, inode = 0; bool released = false;
    explicit Folder (const juce::File& parent)
    {
        struct statfs fs {};
        demand (parent.isDirectory() && ::statfs (parent.getFullPathName().toRawUTF8(), &fs) == 0
                && (fs.f_flags & MNT_LOCAL), ErrorCode::invalidRequest, "Choose an existing local destination folder");
        file = parent.getChildFile (".levelscope-export-" + juce::Uuid().toString());
        demand (::mkdir (file.getFullPathName().toRawUTF8(), 0700) == 0, ErrorCode::writeFailed, "Cannot create owned export staging directory");
        struct stat s {}; demand (::lstat (file.getFullPathName().toRawUTF8(), &s) == 0, ErrorCode::writeFailed, "Cannot identify staging directory");
        device = static_cast<uint64_t> (s.st_dev); inode = static_cast<uint64_t> (s.st_ino);
    }
    ~Folder() { if (!released) removeOwnedFolder (file, device, inode); }
};
void space (const juce::File& parent, uint64_t payload, const pcm24::Hooks& hooks)
{
    const auto available = hooks.availableBytes ? hooks.availableBytes (parent) : parent.getBytesFreeOnVolume();
    demand (available >= 0 && uint64_t (available) >= payload + 80u * 1024u * 1024u,
            ErrorCode::insufficientDiskSpace, "Insufficient destination space for PCM, reports and 80 MiB reserve");
}
std::vector<std::pair<juce::File, juce::String>> originals (const Request& r)
{
    std::vector<std::pair<juce::File, juce::String>> out;
    if (r.current.media.sha256.isNotEmpty())for (const auto& m : r.current.sourceMembers()) out.emplace_back (m.file, m.media.sha256);
    return out;
}
Request checkedReportMedia(const Request& request,JobControl& control)
{
    auto checked=request;
    for(const auto& item:originals(request))
    {
        try{if(pcm24::hashFile(item.first,control)!=item.second)checked.mediaIssue="Programme bytes changed; recorded hashes are historical identities. Remeasure after verification.";}
        catch(const std::exception&){if(control.isCancelled())throw;checked.mediaIssue="Programme member unavailable; recorded hashes are historical identities. Relink and verify before remeasurement.";}
    }
    return checked;
}
void destinationCheck (const Request& r, const juce::File& dest)
{
    demand (dest != juce::File() && !dest.exists() && !dest.isSymbolicLink(), ErrorCode::invalidRequest, "Destination must be a new export/report folder; existing files are never overwritten");
    for (const auto& item : originals (r))
        demand (!samePhysicalFile (item.first, dest), ErrorCode::invalidRequest, "Destination aliases programme media");
    for (const auto& file : {r.sessionFile, r.video.file}) if (file != juce::File())
        demand (!samePhysicalFile (file, dest), ErrorCode::invalidRequest, "Destination aliases session/reference media");
}
void writeText (const juce::File& f, const juce::String& text, JobControl& control, const pcm24::Hooks& hooks, const char* step)
{
    pcm24::checkpoint (hooks, step, control);
    const int fd = ::open (f.getFullPathName().toRawUTF8(), O_WRONLY | O_CREAT | O_EXCL | O_NOFOLLOW, 0600);
    demand (fd >= 0, ErrorCode::writeFailed, "Cannot create report file");
    try
    {
        size_t left = text.getNumBytesAsUTF8(); const auto* data = text.toRawUTF8();
        while (left)
        {
            pcm24::checkpoint (hooks, "report-write", control);
            const auto n = ::write (fd, data, std::min<size_t> (left, 65536));
            if (n < 0 && errno == EINTR) continue;
            demand (n > 0, ErrorCode::writeFailed, "Report write failed"); data += n; left -= static_cast<size_t> (n);
        }
        pcm24::checkpoint (hooks, "report-flush", control);
        demand (::fsync (fd) == 0, ErrorCode::writeFailed, "Report flush failed");
        pcm24::checkpoint(hooks,"report-close",control);
    }
    catch (...) { ::close (fd); ::unlink (f.getFullPathName().toRawUTF8()); throw; }
    demand (::close (fd) == 0, ErrorCode::writeFailed, "Report close failed");
}
using FileHashes = std::vector<std::pair<juce::String, juce::String>>;
FileHashes writeReports (const juce::File& dir, const juce::var& model, bool audio, JobControl& control, const pcm24::Hooks& hooks)
{
    const auto json = reportJson (model), html = reportHtml (model);
    demand (json.getNumBytesAsUTF8() + html.getNumBytesAsUTF8() <= 16u * 1024u * 1024u,
            ErrorCode::invalidRequest, "Report exceeds 16 MiB budget");
    try
    {
        writeText (dir.getChildFile ("report.json"), json, control, hooks, "report-json");
        writeText (dir.getChildFile ("report.html"), html, control, hooks, "report-html");
        FileHashes expected {{"report.json", digest(json)}, {"report.html", digest(html)}};
        if(audio) expected.emplace_back("programme.wav", model["artifact"]["sha256"].toString());
        auto* files = new juce::DynamicObject();
        for (const auto& item : expected) files->setProperty(item.first, item.second);
        auto* manifest = new juce::DynamicObject(); manifest->setProperty ("schema", "levelscope-export-manifest-v1");
        manifest->setProperty ("kind", audio ? "verified-pcm-export" : model.hasProperty("artifact") ? "artifact-report-no-exported-audio" : "measurement-only-no-exported-audio");
        manifest->setProperty ("files_sha256", juce::var (files));
        const auto manifestText=reportJson(juce::var(manifest));
        writeText(dir.getChildFile("manifest.json"),manifestText,control,hooks,"report-manifest");
        expected.emplace_back("manifest.json",digest(manifestText));
        return expected;
    }
    catch (...)
    {
        for (const auto& n : {"report.json", "report.html", "manifest.json"}) ::unlink (dir.getChildFile (n).getFullPathName().toRawUTF8());
        throw;
    }
}
void commitDirectory (const juce::File& stage, const juce::File& dest, const FileHashes& expected, JobControl& control, const pcm24::Hooks& hooks)
{
    const int directory = ::open (stage.getFullPathName().toRawUTF8(), O_RDONLY | O_NOFOLLOW);
    demand (directory >= 0, ErrorCode::writeFailed, "Cannot flush staging directory");
    const int flushed = ::fsync (directory); ::close (directory);
    demand (flushed == 0, ErrorCode::writeFailed, "Staging directory flush failed");
    pcm24::checkpoint (hooks, "publish-rename", control);
    for(const auto& item:expected)
    {
        const auto file=stage.getChildFile(item.first);struct stat st{};
        demand(::lstat(file.getFullPathName().toRawUTF8(),&st)==0&&S_ISREG(st.st_mode)
               &&pcm24::hashFile(file,control)==item.second,ErrorCode::sourceChanged,
               "Final publication bytes differ from measured audio or generated report: "+item.first);
    }
    demand(!control.isCancelled(),ErrorCode::cancelled,"Cancelled before publication");
    demand (::renamex_np (stage.getFullPathName().toRawUTF8(), dest.getFullPathName().toRawUTF8(), RENAME_EXCL) == 0,
            ErrorCode::writeFailed, "Exclusive publication failed; destination may already exist");
    // Atomic rename is the commit point. No partial audio/report set is exposed.
}
template<typename T, typename F> ServiceResult<T> run (JobControl& control, F&& fn)
{
    try { return {fn(), ErrorCode::none, "Completed"}; }
    catch (const Failure& e) { return {{}, control.isCancelled() ? ErrorCode::cancelled : e.code, e.what()}; }
    catch (const std::exception& e) { return {{}, control.isCancelled() ? ErrorCode::cancelled : ErrorCode::writeFailed, e.what()}; }
}
}
Request capture (const StandaloneSessionModel& s)
{
    Request r; r.source = s.measurementContext (MeasurementDomain::source); r.current = s.measurementContext (MeasurementDomain::currentState);
    if (s.source.hasMeasurement && !s.sourceContentInvalid) r.sourceAnalysis = s.source.analysis;
    if (s.currentState.hasMeasurement && !s.sourceContentInvalid && s.currentStateDocument.isInitialized && !s.currentStateNeedsRemeasurement()) r.currentAnalysis = s.currentState.analysis;
    r.applicability = s.programmeApplicability; r.documentId = s.documentId; r.programmeId = s.programmeId;
    r.documentRevision = s.documentRevision; r.unsaved = s.documentDirty(); r.sessionFile = s.documentFile;
    if(s.sourceContentInvalid)r.mediaIssue="Programme media is unresolved or invalidated; recorded hashes are historical identities.";
    else if(r.current.media.sha256.isNotEmpty()&&!r.current.sourceFilesMatch())
        r.mediaIssue="A programme member is missing, changed or duplicated; reverify and remeasure.";
    r.sections = s.sectionDocument.map(); r.video = s.referenceVideo; return r;
}
juce::String specification (const Request& r, Specification& s)
{
    s = {}; const auto& c = r.current;
    auto error = OfflineAnalysisService::validateRequest (c, MeasurementDomain::currentState);
    if (error.isNotEmpty()) return error;
    if(c.media.bytes<=0)return "Verified source byte count required";
    const auto mapped = OfflineAnalysisService::measurementRoles (c.fileLayout);
    using R = measurement::Role;
    const bool mono=c.media.channels==1 && mapped.size()==1
        && (mapped[0]==R::left || mapped[0]==R::right || mapped[0]==R::centre);
    if(c.media.channels==1 && !mono) return "Mono export requires confirmed Left, Right or Centre";
    const uint32_t mask=mono ? (mapped[0]==R::left?1u:mapped[0]==R::right?2u:4u) : 0;
    if (!pcm24::describe (static_cast<int> (c.media.sampleRate), c.media.channels, c.media.frames, s.format, error,mask)) return error;
    const std::vector<R> roles = mono ? mapped : c.media.channels == 2 ? std::vector<R> {R::left, R::right}
        : std::vector<R> {R::left, R::right, R::centre, R::lfe, R::surroundLeft, R::surroundRight};
    for (auto role : roles)
    {
        const auto it = std::find (mapped.begin(), mapped.end(), role);
        if (it == mapped.end() || std::count (mapped.begin(), mapped.end(), role) != 1) return "Export requires explicit mono, stereo or conventional 5.1 roles";
        s.outputToInput.push_back (static_cast<int> (std::distance (mapped.begin(), it)));
    }
    s.layout.count = c.media.channels; s.layout.confirmed = true;
    const auto canonical = c.media.channels == 2 ? juce::AudioChannelSet::stereo() : juce::AudioChannelSet::create5point1();
    for (int n = 0; n < c.media.channels; ++n) s.layout.roles[n] = mono ? c.fileLayout.roles[0] : canonical.getTypeOfChannel (n);
    juce::String key = juce::String (pcm24::version) + ";" + c.media.sha256 + ";" + exact (c.processing.staticGainDb)
        + ";" + c.processing.engineVersion + ";" + juce::String (c.media.frames) + ";" + exact (c.media.sampleRate);
    if(c.processing.sectionGains) key += ";" + c.processing.sectionGains->audioKey();
    for (int n = 0; n < c.fileLayout.count; ++n) key += ";" + juce::String (static_cast<int> (c.fileLayout.roles[n]));
    for (auto n : s.outputToInput) key += ";map=" + juce::String (n);
    s.identity = digest (key); s.seedHex = s.identity.substring (0,16); s.seed = static_cast<uint64_t> (s.seedHex.getHexValue64()); return {};
}
juce::String audioIdentity (const Request& r)
{
    Specification s; if (specification (r, s).isNotEmpty()) return {};
    return digest (juce::String (version) + ";" + s.identity + ";"
        + OfflineAnalysisService::identityFor (r.current, MeasurementDomain::currentState, r.current.media.sha256));
}
juce::String requestIdentity (const Request& r)
{
    return digest (audioIdentity (r) + ";" + r.current.sessionId + ";" + r.documentId + ";" + r.programmeId
        + ";" + juce::String (r.documentRevision) + ";" + juce::String (r.current.ticket.sessionRevision)
        + ";" + r.current.profile.digest() + ";" + juce::String (static_cast<int> (r.applicability))
        + ";" + (r.sections.valid() ? r.sections.digest() : juce::String("no-bound-sections"))
        + ";" + r.video.sha256 + ";" + juce::String(r.video.startOnAudio.text()));
}
juce::String readiness (const Request& r)
{
    if(r.mediaIssue.isNotEmpty())return r.mediaIssue;
    Specification s; auto error = specification (r, s); if (error.isNotEmpty()) return error;
    const auto pair = trimRequest (r.sourceAnalysis.get(), r.currentAnalysis.get(), r.source, r.current, r.current.profile, r.applicability);
    if (!pair.source.currentIdentityMatches || !pair.current.currentIdentityMatches
        || pair.sourceInputIdentity.isEmpty() || pair.sourceInputIdentity != pair.currentInputIdentity)
        return "Valid matching Source/Current measurements are required; re-measure unavailable or stale results.";
    if (proposal::evaluateMeasuredProgramme (r.current.profile, pair.current).verdict != proposal::Verdict::pass)
        return "Current State must pass the selected authoritative programme rules. Measurement-only reports remain available.";
    return {};
}
juce::String preview (const Request& r, bool includeReadiness)
{
    Specification s; const auto error = specification (r,s);
    juce::String text = "Current State PCM export - render, verify actual bytes, then explicitly publish.\n";
    if (error.isNotEmpty()) return text + error;
    text += juce::String (s.format.rf64 ? "RF64" : "RIFF/WAVE") + " PCM24 | " + juce::String(s.format.rate) + " Hz | "
        + juce::String(s.format.frames) + " frames | " + juce::String(s.format.fileBytes) + " bytes\n";
    text += "Output mapping: ";
    for (int n=0;n<s.layout.count;++n) text += (n ? "; " : "") + juce::AudioChannelSet::getAbbreviatedChannelTypeName(s.layout.roles[n])
        + " <- Channel " + juce::String(s.outputToInput[static_cast<size_t>(n)]+1);
    text += "\nTPDF dither once at final quantization. Zero-gain output is not a bit-identical copy.\nNo monitor gain/routing, video audio, clipping, limiting or automatic correction.\n";
    if(!includeReadiness)return text;
    const auto why = readiness(r); return text + (why.isEmpty() ? "Ready to render and verify." : why);
}
Candidate::~Candidate() { if (!published) removeOwnedFolder(directory,directoryDevice,directoryInode); }
bool Candidate::matchesAudio (const Request& r) const
{ return measured && !published && runtimeId == r.current.sessionId && originAudioIdentity.isNotEmpty() && originAudioIdentity == audioIdentity(r); }
proposal::Evaluation artifactEvaluation (const Candidate& c, const Request& r)
{
    if (!c.matchesAudio(r)) return evaluateSnapshot(nullptr,r.current,MeasurementDomain::renderedOutput,r.current.profile,r.applicability);
    auto measured=c.measuredInput();measured.applicability=r.applicability;
    return proposal::evaluateMeasuredProgramme(r.current.profile,measured);
}
juce::String approvalIdentity (const Candidate& c, const Request& r)
{
    if (readiness(r).isNotEmpty() || artifactEvaluation(c,r).verdict != proposal::Verdict::pass) return {};
    return digest(requestIdentity(r)+";"+c.hash()+";"+artifactEvaluation(c,r).identity+";"+c.spec().identity);
}
ServiceResult<std::shared_ptr<Candidate>> render (const Request& r, const juce::File& parent, JobControl& control,
                                                const pcm24::Hooks& hooks, ProcessedAudioObserver* observer, std::shared_ptr<const VerifiedSource> reusableSource)
{
    return run<std::shared_ptr<Candidate>>(control,[&]
    {
        demand(readiness(r).isEmpty(),ErrorCode::invalidRequest,readiness(r));
        Specification spec; demand(specification(r,spec).isEmpty(),ErrorCode::invalidRequest,"Invalid export specification");
        SnapshotStorage storage; storage.availableBytes=hooks.availableBytes;
        struct stat destinationStat{},temporaryStat{};
        demand(::stat(parent.getFullPathName().toRawUTF8(),&destinationStat)==0
               &&::stat(storage.directory.getFullPathName().toRawUTF8(),&temporaryStat)==0,
               ErrorCode::invalidRequest,"Destination and temporary folders must be available");
        const uint64_t copyBytes=reusableSource?0:static_cast<uint64_t>(r.current.media.bytes);
        const bool sharedVolume=destinationStat.st_dev==temporaryStat.st_dev;
        space(parent,spec.format.fileBytes+(sharedVolume?copyBytes:0),hooks);
        if(copyBytes&&!sharedVolume)space(storage.directory,copyBytes,hooks);
        pcm24::checkpoint(hooks,"prepare",control);
        JobControl preparation;preparation.parent=&control;preparation.parentProgress=&control.progress;preparation.progressScale=.2f;
        VerifiedSourceResult prepared;
        if(reusableSource)
        {
            demand(reusableSource->matches(r.current.media),ErrorCode::sourceChanged,"Verified copy is incompatible with this programme");
            const auto copies=reusableSource->files();const auto members=originals(r);
            demand(copies.size()==members.size(),ErrorCode::sourceChanged,"Incomplete reusable programme copy");
            for(size_t i=0;i<copies.size();++i)demand(pcm24::hashFile(copies[i],control)==members[i].second,ErrorCode::sourceChanged,"Reusable verified copy changed");
            prepared.value=reusableSource;preparation.setProgress(1);
        }
        else prepared=prepareVerifiedSource(r.current,preparation,storage);
        demand(prepared.value.has_value(),prepared.error,prepared.message);
        Folder stage(parent); auto c=std::shared_ptr<Candidate>(new Candidate());
        c->directory=stage.file;c->directoryDevice=stage.device;c->directoryInode=stage.inode;stage.released=true;
        c->sourceCopy=*prepared.value;
        c->output=spec;c->runtimeId=r.current.sessionId;c->originAudioIdentity=audioIdentity(r);
        auto reader=(*prepared.value)->createReader(); demand(reader!=nullptr,ErrorCode::decodeFailed,"Cannot decode complete verified programme");
        pcm24::Writer writer(c->file(),spec.format,hooks,control);
        juce::AudioBuffer<float> audio(r.current.media.channels,4096); ProgrammeGainProcessor processing(r.current.processing);
        for(int64_t first=0;first<spec.format.frames;)
        {
            pcm24::checkpoint(hooks,"render-block",control); const int count=static_cast<int>(std::min<int64_t>(4096,spec.format.frames-first));
            demand(reader->read(audio.getArrayOfWritePointers(),audio.getNumChannels(),first,count),ErrorCode::decodeFailed,"Incomplete programme read; no silence substitution");
            demand(processing.process(audio.getArrayOfWritePointers(),audio.getNumChannels(),count,first),ErrorCode::nonFiniteAudio,"Non-finite processed audio");
            if(observer)observer->observe(audio.getArrayOfReadPointers(),audio.getNumChannels(),count,first);
            writer.write(audio.getArrayOfReadPointers(),spec.outputToInput,count,spec.seed,control); first+=count;
            control.setProgress(.2f+.35f*static_cast<float>(double(first)/spec.format.frames));
        }
        writer.finish(control); pcm24::checkpoint(hooks,"verify",control);
        JobControl validationControl;validationControl.parent=&control;validationControl.parentProgress=&control.progress;validationControl.progressBase=.55f;validationControl.progressScale=.1f;
        const auto validation=validateMedia(c->file(),validationControl); demand(validation.succeeded(),validation.error,validation.message);
        ServiceContext artifact=r.current;artifact.assembly.reset();artifact.composition.reset();artifact.sourceFile=c->file();artifact.media=validation.media;
        artifact.fileLayout=spec.layout;artifact.processing={};artifact.range={0,spec.format.frames};
        JobControl measurementControl;measurementControl.parent=&control;measurementControl.parentProgress=&control.progress;measurementControl.progressBase=.65f;measurementControl.progressScale=.35f;
        auto measured=OfflineAnalysisService().verifyEncodedArtifact(artifact,measurementControl);
        demand(measured.value.has_value(),measured.error,measured.message);
        c->sha256=validation.media.sha256;c->measured=std::make_shared<const AnalysisSnapshot>(std::move(*measured.value));
        c->artifactInput=measuredProgrammeForSnapshot(c->measured.get(),c->measured->context,MeasurementDomain::renderedOutput,r.applicability);
        pcm24::checkpoint(hooks,"verified",control);control.setProgress(1);return c;
    });
}
ServiceResult<juce::File> publish (Candidate& c,const Request& fresh,const juce::File& dest,const juce::String& approval,
                                 JobControl& control,const pcm24::Hooks& hooks)
{
    return run<juce::File>(control,[&]
    {
        demand(approval.isNotEmpty()&&approval==approvalIdentity(c,fresh),ErrorCode::staleRevision,"Export approval is missing or stale");
        destinationCheck(fresh,dest);
        demand(dest.getParentDirectory()==c.directory.getParentDirectory()&&sameFolder(c.directory,c.directoryDevice,c.directoryInode),
               ErrorCode::invalidRequest,"Publish to a new folder in the chosen staging filesystem");
        pcm24::checkpoint(hooks,"publish-recheck",control);
        for(const auto& item:originals(fresh))demand(pcm24::hashFile(item.first,control)==item.second,ErrorCode::sourceChanged,"Programme member changed before publication");
        demand(pcm24::hashFile(c.file(),control)==c.hash(),ErrorCode::sourceChanged,"Artifact changed after verification");
        space(dest.getParentDirectory(),0,hooks);
        const auto expected=writeReports(c.directory,report(fresh,&c,true),true,control,hooks);
        try
        {
            commitDirectory(c.directory,dest,expected,control,hooks);
        }
        catch(...){for(const auto& n:{"report.json","report.html","manifest.json"})::unlink(c.directory.getChildFile(n).getFullPathName().toRawUTF8());throw;}
        c.directory=dest;c.published=true;control.setProgress(1);return dest;
    });
}
ServiceResult<juce::File> measurementReport (const Request& r,const juce::File& dest,JobControl& control,const pcm24::Hooks& hooks)
{
    return run<juce::File>(control,[&]
    {
        destinationCheck(r,dest);space(dest.getParentDirectory(),0,hooks);Folder stage(dest.getParentDirectory());
        pcm24::checkpoint(hooks,"measurement-report",control);
        const auto validated=checkedReportMedia(r,control);
        const auto expected=writeReports(stage.file,report(validated),false,control,hooks);
        commitDirectory(stage.file,dest,expected,control,hooks);stage.released=true;control.setProgress(1);return dest;
    });
}
ServiceResult<juce::File> artifactReport(const Candidate& candidate,const Request& request,const juce::File& dest,JobControl& control,const pcm24::Hooks& hooks)
{
    return run<juce::File>(control,[&]
    {
        demand(candidate.matchesAudio(request),ErrorCode::staleRevision,"Staged artifact no longer matches this programme/runtime");
        demand(pcm24::hashFile(candidate.file(),control)==candidate.hash(),ErrorCode::sourceChanged,"Staged artifact bytes changed after measurement");
        destinationCheck(request,dest);space(dest.getParentDirectory(),0,hooks);Folder stage(dest.getParentDirectory());
        const auto expected=writeReports(stage.file,report(checkedReportMedia(request,control),&candidate),false,control,hooks);
        commitDirectory(stage.file,dest,expected,control,hooks);stage.released=true;control.setProgress(1);return dest;
    });
}

}
