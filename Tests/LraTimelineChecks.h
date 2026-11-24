// Included by DiagnosticsTests.cpp inside its fixture namespace.
namespace lt = proposal::lra_timeline;

void finishLraFixture(measurement::Result& r)
{
    const double abs = energy(-70.);
    double sum = 0;
    size_t count = 0;
    for (const auto& w : r.lraWindows)
        if (w.energy >= abs) { sum += w.energy; ++count; }
    const double relative = sum / static_cast<double>(count) / 100.;
    std::vector<double> sorted;
    for (const auto& w : r.lraWindows)
        if (w.energy >= abs && w.energy >= relative)
            sorted.push_back(measurement::OfflineMeasurement::loudness(w.energy));
    std::sort(sorted.begin(),sorted.end());
    r.lraRelativeGate = measurement::OfflineMeasurement::loudness(relative);
    r.lraIncluded = sorted.size();
    const auto order = [&](double p) { return sorted.at(static_cast<size_t>(std::floor((sorted.size()-1)*p+.5))); };
    r.p10 = order(.10); r.p95 = order(.95);
    r.lra = {*r.p95-*r.p10, measurement::Status::available};
    r.lraStable = r.frames >= int64_t(r.sampleRate)*60;
}
proposal::TargetProfile targetProfile(double ceiling)
{
    auto p = ebuR128ProgrammeQc();
    p.limits.hasMaxLoudnessRange = true; p.limits.maxLoudnessRangeLu = ceiling; p.lraProvenance = "Fixture job ceiling";
    return p;
}
void measuredLaneDefinitions()
{
    // The curves and peaks are from actual PCM; replace only the LRA sequence in this
    // explicit estimator fixture to exercise exact gates, ties and drain membership.
    auto r = signal(44100,1,61.037,0);
    size_t n=0;
    for (auto& w : r.lraWindows)
    {
        const double levels[]{-14,-14,-14,-14,-14,-22,-22,-22,-22,-22,-22,-22,-22,-22,-22,
                              -30,-30,-30,-30,-30,-34,-52,-80};
        w.energy = n % 57 == 0 ? energy(-9.) : n % 59 == 0 ? 0 : energy(levels[n % 23]);
        ++n;
    }
    finishLraFixture(r);
    const auto result = analyzeProgramme(r,"tie-and-gate-fixture");
    const auto index = lt::Index::build(result);
    const auto app = Applicability::preProducedLinearProgramme;
    auto target = index->target(targetProfile(10.),app);
    expect(index && target.state==lt::TargetState::exceeds && target.upperLevelLufs==*r.p10+10,
           "target points to P10 plus the selected ceiling only when stable measured LRA exceeds that ceiling");
    expect(index->upperRankCount()==r.lraIncluded-1-static_cast<size_t>(std::floor((r.lraIncluded-1)*.95+.5)),
           "upper rank count excludes the P95 rank and uses the qualified nearest-rank convention");
    unsigned observed=0; size_t upper=0, strict=0, ties=0;
    std::vector<std::pair<double,size_t>> order;
    for (size_t i=0;i<n;++i) if (result->lraWindows[i].included) order.emplace_back(r.lraWindows[i].energy,i);
    std::sort(order.begin(),order.end());
    std::vector<bool> rankTail(n,false);
    for (size_t i=order.size()-index->upperRankCount();i<order.size();++i) rankTail[order[i].second]=true;
    bool exact=true, firstTieSeen=false, lastTieSeen=false;
    for(size_t i=0;i<n;++i)
    {
        const auto& w=result->lraWindows[i];
        const auto windows=index->inspect(w.begin,target);
        const auto match=std::find_if(windows.begin(),windows.end(),[&](const auto& v){return v.index==i;});
        if(match==windows.end()){exact=false;continue;}
        const auto& v=*match;
        exact &= v.lufs==w.lufs && v.begin==w.begin && v.end==std::min(r.frames,w.end) && v.measuredEnd==w.end;
        exact &= bool(v.flags & (lt::aboveP95|lt::upperRankTie))==rankTail[i];
        if(w.included)
        {
            exact &= bool(v.flags & lt::aboveP95)==(*w.lufs>*r.p95);
            exact &= bool(v.flags & lt::belowP10)==(*w.lufs<*r.p10);
            exact &= bool(v.flags & lt::targetPassage)==(*w.lufs>*target.upperLevelLufs);
            if(w.lufs==r.p95){firstTieSeen |= !rankTail[i]; lastTieSeen |= rankTail[i];}
        }
        else
        {
            exact &= (v.flags & (lt::eligible|lt::aboveP95|lt::upperRankTie|lt::belowP10|lt::targetPassage))==0;
            exact &= bool(v.flags & lt::absoluteRejected)==(w.energy<energy(-70));
            exact &= bool(v.flags & lt::relativeRejected)==(w.energy>=energy(-70));
        }
        observed|=v.flags;
        upper+=bool(v.flags&(lt::aboveP95|lt::upperRankTie)); strict+=bool(v.flags & lt::aboveP95);ties+=bool(v.flags & lt::upperRankTie);
    }
    expect(exact && upper==index->upperRankCount(),"every full-precision level, support interval, percentile rank and gate matches the independent sorted fixture");
    expect(strict>0 && ties>0 && firstTieSeen && lastTieSeen,"equal P95 observations across the cut are identified; some upper ranks are equally loud");
    performance->setProperty("eligible_windows",static_cast<int>(r.lraIncluded));
    performance->setProperty("upper_ranks",static_cast<int>(upper));
    performance->setProperty("strictly_above_p95",static_cast<int>(strict));
    performance->setProperty("upper_rank_ties_at_p95",static_cast<int>(ties));
    performance->setProperty("p10_lufs",*r.p10);
    performance->setProperty("p95_lufs",*r.p95);
    performance->setProperty("index_bytes_excluding_existing_measurement",static_cast<int>(index->indexBytes()));
    constexpr unsigned all = lt::absoluteRejected|lt::relativeRejected|lt::zeroEnergy|lt::belowP10|lt::eligible|lt::atP10|lt::atP95|lt::padded|lt::targetPassage;
    expect((observed&all)==all,"fixture distinguishes silence, both gates, strict low, boundary ties, padded windows and target-relative passages");
    const auto profile=targetProfile(*r.lra.value);
    const auto within=index->target(profile,app);
    expect(within.state==lt::TargetState::within && !(index->columns(0,r.frames,1,within)[0].flags & lt::targetPassage),
           "a programme exactly at the ceiling emits no excessive-LRA passage");
    auto none=index->target(ebuR128ProgrammeQc(),app);
    expect(none.state==lt::TargetState::noCeiling && !none.upperLevelLufs,"no selected LRA ceiling has no target-relative assertion");
    expect(index->columns(0,r.frames,1,none)[0].flags==(observed&~lt::targetPassage),"target changes preserve all gate and percentile evidence");
    auto invalid=targetProfile(std::numeric_limits<double>::infinity());
    expect(index->target(invalid,app).state==lt::TargetState::inapplicable &&
           index->target(profile,Applicability::unconfirmed).state==lt::TargetState::inapplicable,"invalid and inapplicable requirements cannot publish a target-relative lane");
    auto foreign=target;foreign.snapshotIdentity="old measurement";
    expect(!(index->columns(0,r.frames,1,foreign)[0].flags & lt::targetPassage),"a foreign target interpretation cannot mark this measurement");
    expect(!lt::Index::build(result,[]{return true;}),"cancelled construction cannot publish an index");
    const auto hop=r.sampleRate/10;
    expect(index->inspect(0,target).size()==1 && index->inspect(29*hop,target).size()==30 &&
           index->inspect(30*hop,target).front().index==1,"half-open support enters and leaves on exact source frames across 30 overlapping windows");
    expect(index->inspect(-1,target).empty() && index->inspect(r.frames,target).empty(),"inspection never extends into synthetic drain time");
    const auto end=index->inspect(r.frames-1,target);
    expect(!end.empty() && end.back().end==r.frames && end.back().measuredEnd>r.frames &&
           (end.back().flags & lt::padded),"non-hop-aligned EOF clips display spans and preserves measured padded end");
    bool columnsExact=true;
    for (int pixels : {1,11,147,1024})
    {
        const int64_t begin=1251, limit=r.frames-135;
        const auto cols=index->columns(begin,limit,pixels,target);
        for(int x=0;x<pixels;++x)
        {
            const auto left=begin+(limit-begin)*x/pixels, right=begin+((limit-begin)*(x+1)+pixels-1)/pixels;
            unsigned bits=0;
            std::optional<double> low,high;
            for(const auto& v : result->lraWindows)
                if(v.begin<right && v.end>left && v.lufs)
                {if(!low||*v.lufs<*low)low=v.lufs;if(!high||*v.lufs>*high)high=v.lufs;}
            // A window that starts inside this pixel or overlaps its first frame contributes its flags.
            for(auto f=left;f<right;f=(f/hop+1)*hop)
                for(const auto& v:index->inspect(f,target))bits|=v.flags;
            columnsExact &= cols[x].minLufs==low && cols[x].maxLufs==high && cols[x].flags==bits;
        }
    }
    expect(columnsExact,"bounded viewport index equals a full-resolution overlap union, including mixed categories and extrema");
    expect(index->columns(-1,r.frames+1,5000,target).size()==4096 && index->columns(0,0,8,target).empty(),"viewport limits, clipping and empty spans are explicit");
    const auto atFrame=index->inspect(2345,target);
    const auto subframe=index->columns(2345,2346,32,target);
    expect(subframe.size()==32 && subframe.front().flags==atFrame.front().flags && subframe.back().flags==subframe.front().flags,"subframe-scale pixels preserve the single source frame's evidence");
    std::fill(r.lraWindows.begin(),r.lraWindows.end(),measurement::Window{}); // restore an exact regular clock below
    for(size_t i=0;i<n;++i)r.lraWindows[i]={3*r.sampleRate+int64_t(i)*hop,0};
    r.lraWindows[0].energy=1.;r.lraWindows[1].energy=std::nextafter(1.,0.);
    r.lraWindows[2].energy=1.;r.lraWindows[3].energy=397.;
    finishLraFixture(r);
    auto edge=lt::Index::build(analyzeProgramme(r,"gate-edge"));
    const auto gate=edge->inspect(hop,none);
    expect((gate[0].flags & lt::eligible) && (gate[1].flags & lt::relativeRejected),"exact relative gate included, one ULP below rejected without LUFS rounding");
    r.lraWindows[0].energy=measurement::absoluteGateEnergy();
    r.lraWindows[1].energy=std::nextafter(measurement::absoluteGateEnergy(),0.);
    r.lraWindows[2].energy=r.lraWindows[3].energy=0.;
    finishLraFixture(r);
    edge=lt::Index::build(analyzeProgramme(r,"absolute-edge"));
    const auto abs=edge->inspect(hop,none);
    expect((abs[0].flags & lt::eligible) && (abs[1].flags & lt::absoluteRejected),"exact absolute gate included and one ULP below rejected");
    for(auto& w:r.lraWindows)w.energy=energy(-23);
    finishLraFixture(r);
    auto flat=lt::Index::build(analyzeProgramme(r,"all-equal"));
    auto flags=flat->columns(0,r.frames,1,flat->target(targetProfile(0),app))[0].flags;
    expect(flat->upperRankCount()>0 && (flags & lt::upperRankTie) && !(flags&(lt::aboveP95|lt::belowP10|lt::targetPassage)),
           "constant levels have upper ranks but no distinctly louder, quieter or excessive passage");
    r.lraWindows[0].energy=std::numeric_limits<double>::quiet_NaN();
    rejects([&]{analyzeProgramme(r,"non-finite window");},"non-finite measured input is rejected before lane construction");
    r.lraWindows[0].energy=energy(-23);
    r.lraStable=false;
    rejects([&]{lt::Index::build(analyzeProgramme(r,"bad-stability"));},"inconsistent saved stability cannot enable or suppress a target verdict");
}
void measuredLaneAvailability()
{
    const auto profile=targetProfile(2);
    for(int channels : {1,2,6})
    {
        const auto r=signal(48000,channels,4.147,0);
        const auto a=analyzeProgramme(r,"pcm-"+juce::String(channels));
        const auto index=lt::Index::build(a);
        const auto t=index->target(profile,Applicability::preProducedLinearProgramme);
        expect(a->p10==r.p10 && a->p95==r.p95 && index->measured().lra.value==r.lra.value &&
               t.state==lt::TargetState::unstable && !t.upperLevelLufs,"measured mono/stereo/5.1 retains its exact short-programme percentiles without a target assertion");
        expect(index->columns(0,r.frames,1,t)[0].flags & lt::eligible,"completed unstable measurement still supplies descriptive gate/percentile information");
    }
    for(int kind : {1,2})
    {
        auto r=signal(48000,1,4.01,kind);
        const auto index=lt::Index::build(analyzeProgramme(r,"quiet-pcm"));
        const auto t=index->target(profile,Applicability::preProducedLinearProgramme);
        auto c=index->columns(0,r.frames,1,t).front();
        expect(t.state==lt::TargetState::unavailable && !r.p10 && !r.p95 && (c.flags & lt::absoluteRejected) && !(c.flags & lt::eligible),"silence/below-gate programmes have no percentile population");
        expect(bool(c.flags & lt::zeroEnergy)==(kind==1) && c.minLufs.has_value()==(kind!=1),"zero energy has no fabricated finite level; finite below-gate loudness is retained");
    }
    for(double duration : {.5,2.9})
    {
        auto r=signal(48000,1,duration,0);
        auto index=lt::Index::build(analyzeProgramme(r,"short-pcm"));
        auto cols=index->columns(0,r.frames,1,{});
        expect(!r.lra.value && (duration<1 ? !cols.front().flags : bool(cols.front().flags & lt::insufficientDuration)) &&
               !(cols.front().flags&(lt::eligible|lt::relativeRejected|lt::absoluteRejected)),
               "short material and drain-only windows never acquire gate or percentile verdicts through padding");
    }
}
void measuredLaneLifecycle()
{
    using State=TimelineMeasurementState;
    Fixture f; StandaloneSessionModel s; JobControl token;
    const auto domain=MeasurementDomain::source, current=MeasurementDomain::currentState;
    expect(timelineMeasurement(s,domain).state==State::unavailable,"empty session has no measured lanes");
    const auto media=validateMedia(f.file,token);
    expect(media.succeeded()&&s.setValidatedSource(f.file,media.media),"import a short real stereo fixture for lane identity checks");
    s.setChannelRole(0,juce::AudioChannelSet::left);s.setChannelRole(1,juce::AudioChannelSet::right);
    expect(s.confirmLayout(),"confirm source roles for diagnostic input");
    expect(timelineMeasurement(s,domain).state==State::unmeasured && timelineMeasurement(s,domain,true).state==State::measuring,
           "provisional analysis is distinguishable and publishes no completed lane");
    measured(s,domain); measured(s,current);
    const auto original=s.source.analysis, applied=s.currentState.analysis;
    DiagnosticsController controller(true);
    auto select=[&](MeasurementDomain d){controller.select(timelineMeasurement(s,d).snapshot,s.measurementContext(d),d);};
    select(current);awaitDiagnostics(controller);
    const auto first=controller.result();
    expect(first && first->lraTimeline && first->snapshot==applied && timelineMeasurement(s,current).state==State::measured,"accepted Current owns its validated, measurement-specific index");
    const auto copies=s.verifiedAssets->counts();
    const auto peak=s.currentState.analysis->detail.truePeakEighthFrames;
    s.setSelectedTargetProfileIndex(1);
    auto p=s.getSelectedTargetProfile();const auto target=first->lraTimeline->target(p,s.programmeApplicability);
    auto settings=s.jobRequirements;settings.maximumLraTenths=31;juce::String error;
    expect(s.setJobRequirements(settings,error),"change only the selected LRA requirement");
    p=s.getSelectedTargetProfile();auto updated=first->lraTimeline->target(p,s.programmeApplicability);
    select(current);controller.poll();
    expect(!target.ceilingLu && updated.ceilingLu==3.1 && controller.result()==first && controller.buildsStarted==1,
           "changing the ceiling reinterprets the same index without a new analysis");
    expect(s.source.analysis==original && s.currentState.analysis==applied && s.currentState.analysis->detail.truePeakEighthFrames==peak &&
           s.verifiedAssets->counts().copies==copies.copies,"lane queries preserve both measurements, peak positions and verified source copy count");
    expect(s.setCurrentGainDb(-1),"deliberately edit Current processing");
    select(current);
    expect(!controller.result() && !timelineMeasurement(s,current).snapshot && timelineMeasurement(s,current).state==State::unverified &&
           timelineMeasurement(s,domain).snapshot==original,"an edit immediately invalidates the Current lane and keeps the Original");
    controller.select(applied,s.measurementContext(current),current);awaitDiagnostics(controller);
    expect(!controller.result(),"late old Current cannot publish after processing changes");
    select(domain);awaitDiagnostics(controller);
    expect(controller.result() && controller.result()->snapshot==original && controller.result()->analysis->snapshotIdentity!=first->analysis->snapshotIdentity,
           "domain switching uses independently identified measurements");
    auto other=s.measurementContext(domain);other.sessionId="another runtime";controller.select(original,other,domain);
    expect(!controller.result(),"the same immutable pointer cannot retain acceptance in a different runtime");
    awaitDiagnostics(controller);expect(!controller.result(),"foreign runtime completion cannot publish");
    s.sourceContentInvalid=true;
    expect(!timelineMeasurement(s,domain).snapshot && timelineMeasurement(s,domain).state==State::unavailable,"a changed source immediately removes the previous measured layer");
    expect(!timelineMeasurement(s,MeasurementDomain::privateCandidate).snapshot,"a privately measured candidate cannot be a timeline domain");
}
