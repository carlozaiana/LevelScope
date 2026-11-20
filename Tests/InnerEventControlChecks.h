// Deterministic feature-cell controls isolate grouping from the PCM front end.
enum class InnerControl { recurring, equalHits, click, weak, missingPredecessor, missingCurrent, conflictingSource, maskedRole, clipped, differentTail };
EventEvidence controlEvidence(InnerControl mode,int hop,int shift)
{
    const int rate=48000,step=rate*hop/1000;
    const auto p=unity(7,rate);auto a=analysis(p);const int entry=2*rate+shift*rate/1000;
    measurement::Result src;src.sampleRate=rate;src.frames=p.frames;
    for(int t=rate*4/10;t<=p.frames;t+=rate/10)
    {
        const bool high=t>=entry+rate*4/10&&t<=entry+rate*2;
        DiagnosticWindow w;w.begin=t-rate*4/10;w.end=t;w.lufs=high&&mode!=InnerControl::click?-15.:-25.;
        if(mode!=InnerControl::missingCurrent||t!=((entry+rate-1)/(rate/10))*(rate/10))a->momentary.push_back(w);
        src.momentary.push_back({t,high&&mode!=InnerControl::conflictingSource?.01:.001});
    }
    auto c=area(p,{rate,6*rate});c.settings=RegionPlanningSettings::innerDefaults(RegionMode::automatic,p.inputIdentity);c.settings.eventFeatures.hopMs=hop;
    EventEvidence evidence;evidence.roles={measurement::Role::left,measurement::Role::right};evidence.frames=p.frames;evidence.incompleteRoots={false};evidence.retainedContexts={{0,p.frames}};
    EventGroup parent;parent.root=0;parent.id="continuing";evidence.groups={parent};
    InnerEventTracker tracker(c,*a,&src,5,evidence);
    const int first=entry-rate*8/10,second=entry-rate*4/10,last=entry+rate*8/10,tail=last+rate*3/10;
    for(int t=step;t<=p.frames;t+=step)
    {
        const auto pulse=[&](int edge){return t>=edge&&t<=edge+rate/10;};
        const bool earlier=pulse(first)||pulse(second),inner=pulse(entry)||pulse(last),trailing=mode==InnerControl::differentTail&&pulse(tail);
        const bool attack=t==first||(t==second&&mode!=InnerControl::missingPredecessor)||t==entry||t==last||(trailing&&t==tail);
        double energy=earlier?.001:inner?(mode==InnerControl::equalHits?.001:mode==InnerControl::weak?.002:.01):trailing?.0005:.0001;
        InnerEventTracker::Bands fast{},background{};fast[0]=fast[1]=energy;background.fill(.0001);
        if(mode==InnerControl::maskedRole){fast[1]=inner?.1:.0001;fast[0]=.25;energy=(fast[0]+fast[1])/2;}
        EventCell cell;cell.begin=t-step;cell.end=t;cell.broadband=cell.low=energy;
        tracker.begin(cell,t-1,1,fast,fast,background,background);
        if(mode!=InnerControl::clipped||t<entry+rate/20)tracker.observe(0,0,attack,attack,attack,t-step);
        tracker.advance();
    }
    tracker.finish(nullptr,[]{return false;});return evidence;
}
void innerControls()
{
    for(int hop:{5,10,20})for(int shift:{0,400})
    {
        const auto e=controlEvidence(InnerControl::recurring,hop,shift);
        const int frame=int(48000*(2000+shift)/1000);
        expect(e.groups.size()==2&&e.groups.back().entryDecision==frame&&e.groups.back().lastAttack.end==frame+38400&&e.groups.back().endingSupported,"shifted entrance, internal trough and final comparable hit across hop settings");
        expect(!e.groups.back().innerNominee,"missing bound Current sidecar excludes nomination");
    }
    for(auto control:{InnerControl::equalHits,InnerControl::click,InnerControl::weak,InnerControl::missingPredecessor,InnerControl::clipped})
        expect(controlEvidence(control,10,0).groups.size()==1,"equal recurring hits, isolated click, weak pickup, absent predecessor or clipped detail cannot open a new passage");
    const auto missing=controlEvidence(InnerControl::missingCurrent,10,0);
    expect(missing.groups.size()==1,"incomplete contextual Current population cannot corroborate a step");
    const auto conflict=controlEvidence(InnerControl::conflictingSource,10,0);
    expect(conflict.groups.size()==2&&conflict.groups.back().reason.contains("Source momentary"),"Source and Current disagreement is explicit");
    const auto role=controlEvidence(InnerControl::maskedRole,10,0);
    expect(role.groups.size()==2&&role.groups.back().reason.contains("role-specific")&&role.groups.back().incomplete,"masked non-LFE witness stays tentative");
    const auto tail=controlEvidence(InnerControl::differentTail,10,0);
    expect(tail.groups.size()==2&&tail.groups.back().lastAttack.end==134400&&tail.groups.back().endingSupported,"weaker trailing sound does not replace the final comparable hit");
#if defined(__APPLE__)
    event_allocation::Scope meter;
    auto ptr=std::malloc(100003);if(!ptr)throw std::bad_alloc();const size_t one=malloc_size(ptr);
    auto second=std::malloc(12347);if(!second)throw std::bad_alloc();const size_t two=malloc_size(second);
    std::free(ptr);std::free(second);
    const auto measured=meter.stop();expect(meter.valid()&&measured==one+two,"allocator probe counts simultaneous blocks at their allocator sizes");
#endif
}
