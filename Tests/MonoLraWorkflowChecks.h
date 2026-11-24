// Native Propose on a single physical channel, from explicit role confirmation
// through private assessment, applied Current and the finalized PCM24 file.
void monoSelectedRanges(StandaloneSessionModel s)
{
    const auto rate = s.source.analysis->detail.sampleRate;
    const RegionSpan selection {30LL * rate, 36LL * rate};
    auto requirements = s.jobRequirements;
    requirements.maximumLraTenths = 49;
    juce::String error;
    expect(s.setJobRequirements(requirements,error), "set mono selected-range requirement");
    auto& sections = s.sectionDocument;
    expect(sections.split(sections.map().sections.front().id, selection.begin, error)
           && sections.split(sections.map().sections.back().id, selection.end, error),
           "create a six-second editorial section around the known loud passage");
    for (const auto defaults : {RegionPlanningSettings::defaults, RegionPlanningSettings::eventDefaults, RegionPlanningSettings::innerDefaults})
    {
        for (auto mode : {RegionMode::selectedSection, RegionMode::selection})
        {
            auto settings = defaults(mode, sectionInputIdentity(s.measurementContext(MeasurementDomain::source)));
            if (mode == RegionMode::selectedSection) settings.sectionId = sections.map().sections[1].id;
            else settings.selection = selection;
            settings.leftMarginMs = 1000;
            settings.rightMarginMs = 1000;
            settings.allowGlobal = true;
            expect(s.setRegionPlanning(settings,error), "accept mono selected scope and explicit margins");
            JobControl job;
            const auto result = generateSectionProposal(s.sectionProposalRequest(),job);
            expect(result.value && (*result.value)->methodVersion() == settings.methodVersion() &&
                   (*result.value)->passes() > 0 && (*result.value)->passes() <= 3,
                   "selected section and temporary selection both reach complete private mono assessment");
            const auto review = (*result.value)->regionReview();
            expect(review && review->constraints.roots.size() == 1 && review->constraints.roots[0].hint == selection,
                   "selected mono method retains the exact requested frames");
            const auto allowed = review->constraints.roots.front().allowed;
            expect(allowed.begin == selection.begin-rate && allowed.end == selection.end+rate, "mono extensions respect each one-second allowance");
            for (const auto& attempt : review->attempts)
                for (const auto& span : attempt.plan.activity.support)
                    expect(span.begin >= allowed.begin && span.end <= allowed.end, "local correction stays within the permitted mono selection");
        }
    }
}

void monoLraWorkflowChecks(const juce::File& root)
{
    using C = juce::AudioChannelSet;
    using Phase = ProposalOperation::Phase;
    StandaloneMainComponent ui(emptyOutputCatalog(), noOutputSettings());
    Access::stopTimer(ui);
    ui.setSize(1320, 1000);
    auto& s = Access::model(ui);
    auto& op = Access::operation(ui);
    const auto original = programme(root.getChildFile("source"), 1);
    const auto source = original.source.analysis;
    JobControl job;
    juce::String error;
    const std::array<ChannelRole, 4> roles { C::centre, C::left, C::right, C::centre };
    const std::array<juce::String, 4> methods { lraTrimMethodVersion, regionSearchVersion, eventSearchVersion, innerSearchVersion };
    for (int schema = 0; schema < 4; ++schema)
    {
        auto dir = root.getChildFile("method-" + juce::String(schema));
        expect(dir.createDirectory().wasOk(), "new mono method results directory");
        installWorkflowFixture(ui, original);
        s.setChannelRole(0, roles[schema]);
        expect(s.confirmLayout(), "explicit single-channel role confirmed");
        measuredComposition(s, MeasurementDomain::source);
        if (!s.currentState.hasMeasurement)
        {
            expect(s.initializeCurrentStateFromSource(), "initialize mono Current after confirming its role");
            if (!s.currentState.hasMeasurement) measuredComposition(s, MeasurementDomain::currentState);
        }
        s.synchronizeSections();
        const auto& baseline = s.source.analysis->detail;
        expect(baseline.integrated.value == source->detail.integrated.value && baseline.lra.value == source->detail.lra.value
               && baseline.channelTruePeaks == source->detail.channelTruePeaks && baseline.truePeakEighthFrames == source->detail.truePeakEighthFrames,
               "Left, Right and Centre have identical single-channel loudness and peak interpretation");
        expect(baseline.channelSamplePeaks.size() == 1 && baseline.channelTruePeaks.size() == 1
               && baseline.truePeakEighthFrames.size() == 1 && std::abs(baseline.channelSamplePeaks.front() - 0.06) < 0.000001
               && baseline.integrated.value && *baseline.integrated.value < -31 && *baseline.integrated.value > -34,
               "calibrated mono input has one uncompensated peak and single-channel integrated loudness");
        auto req = s.jobRequirements;
        req.maximumLraTenths = 49;
        expect(s.setJobRequirements(req, error), "enable measured mono LRA ceiling");
        if (schema != 0)
        {
            const auto defaults = schema == 1 ? RegionPlanningSettings::defaults
                                : schema == 2 ? RegionPlanningSettings::eventDefaults : RegionPlanningSettings::innerDefaults;
            auto settings = defaults(RegionMode::automatic, sectionInputIdentity(s.measurementContext(MeasurementDomain::source)));
            settings.compareAlternatives = true;
            expect(s.setRegionPlanning(settings, error), "choose existing boundary method on mono");
        }
        const auto before = s.manualState();
        Access::refresh(ui);
        click(Access::proposal(ui), "Propose");
        completeOperation(ui);
        if (op.result) dir.getChildFile("search.txt").replaceWithText(op.result->text());
        expect(op.phase == Phase::verified && s.lastAppliedExplanation && s.lastAppliedExplanation->method == methods[schema],
               "selected mono method automatically applies and independently verifies");
        const auto& result = s.currentState.analysis->detail;
        expect(s.currentState.analysis->domain == MeasurementDomain::currentState && s.evaluation(MeasurementDomain::currentState).verdict == Verdict::pass
               && result.lra.value && *result.lra.value < *baseline.lra.value && *result.lra.value <= 4.9
               && s.currentStateDocument.sectionGains && result.channelTruePeaks.size() == 1
               && s.lastAppliedExplanation->candidate && s.lastAppliedExplanation->candidate->identity != s.currentState.analysis->identity,
               "local mono correction reduces LRA and independent Current owns PASS");
        if (schema == 1)
        {
            expect(op.canChoose(s) && op.choices.size() > 1, "distinct measured mono alternatives available");
            const size_t other = (op.selectedChoice + 1) % op.choices.size();
            choose(Access::proposal(ui), "Measured passing alternatives", int(other + 1));
            completeOperation(ui);
            expect(op.phase == Phase::verified && op.selectedChoice == other && s.lastAppliedExplanation->beforeIdentity == proposalPlanIdentity(op.baseline.currentContext.processing)
                   && s.currentStateDocument.sectionGains->audioKey() == op.selected->plan().sectionGains->audioKey()
                   && s.currentStateDocument.staticGainDb == op.selected->plan().staticGainDb,
                   "alternative replaces the absolute mono plan from its common original baseline");
        }

        const auto saved = dir.getChildFile("passing-mono.levelscope");
        const auto history = encodeAppliedExplanation(*s.lastAppliedExplanation);
        expect(sf::save(sf::capture(s), saved, job).succeeded(), "save passing mono session and full explanation");
        auto document = sf::read(saved, job);
        expect(bool(document.document), "read the saved mono session");
        auto reopened = sf::verify(*document.document, job);
        expect(reopened.prepared && reopened.prepared->model.currentState.hasMeasurement && !reopened.prepared->model.sectionProposal
               && reopened.prepared->model.sourceDocument.layout.roles[0] == roles[schema]
               && reopened.prepared->model.regionPlanning.has_value() == bool(schema)
               && (!schema || reopened.prepared->model.regionPlanning->digest() == s.regionPlanning->digest())
               && reopened.prepared->model.lastAppliedExplanation && encodeAppliedExplanation(*reopened.prepared->model.lastAppliedExplanation) == history,
               "reopen exact mono role, settings, measured Current and explanation without live Apply authority");

        const auto exportRequest = ex::capture(s);
        auto rendered = ex::render(exportRequest, dir, job);
        expect(rendered.value && ex::artifactEvaluation(**rendered.value, exportRequest).verdict == Verdict::pass,
               "encoded mono artifact independently measured and passing");
        auto artifact = *rendered.value;
        const uint32_t mask = roles[schema] == C::left ? 1u : roles[schema] == C::right ? 2u : 4u;
        expect(artifact->spec().format.channels == 1 && artifact->spec().format.mask == mask
               && artifact->spec().format.frames == baseline.frames && artifact->spec().outputToInput == std::vector<int>{0}
               && artifact->analysis().domain == MeasurementDomain::renderedOutput && artifact->analysis().detail.truePeakEighthFrames.size() == 1,
               "finalized audio has one channel, the original confirmed role and one peak position");
        const auto destination = dir.getChildFile("exported");
        auto delivery = ex::publish(*artifact, exportRequest, destination, ex::approvalIdentity(*artifact, exportRequest), job);
        expect(delivery.value && juce::SHA256(destination.getChildFile("programme.wav")).toHexString() == artifact->hash(),
               "published WAV bytes match the measured mono export");

        if (schema == 0)
        {
            const auto applied = s.manualState();
            expect(s.undoManual() && s.manualState() == before && !s.currentState.hasMeasurement,
                   "Undo restores original mono baseline without restoring PASS");
            expect(s.undoManual(true) && s.manualState() == applied && !s.currentState.hasMeasurement,
                   "Redo restores the exact plan and requires measurement");
        }
        if (schema == 1)
            expect(op.recover(s) && s.manualState() == before && !s.currentState.hasMeasurement,
                   "restore before Propose recovers the original mono baseline");
        if (schema == 3)
        {
            const auto cut = s.currentStateDocument.sectionGains->intervals[1].begin;
            expect(s.beginAppliedBoundary(1) && s.previewAppliedBoundary(s.appliedBoundaryDrag->token, cut + 20)
                   && s.commitAppliedBoundary(s.appliedBoundaryDrag->token), "manually refine an applied mono boundary");
            op.observe(s);
            expect(op.phase == Phase::stale && !op.canChoose(s) && !s.currentState.hasMeasurement
                   && encodeAppliedExplanation(*s.lastAppliedExplanation) == history, "editing invalidates alternatives and retains the historical explanation");
        }
    }
    monoSelectedRanges(original);
    installWorkflowFixture(ui, original);
    auto strict = s.jobRequirements;
    strict.maximumLraTenths = 0;
    expect(s.setJobRequirements(strict, error), "deliberately strict mono ceiling");
    const auto before = s.manualState();
    Access::workflow(ui);
    Access::cancel(ui);
    completeOperation(ui);
    expect(op.phase == Phase::cancelled && !op.applied && s.manualState() == before, "cancelled mono Propose preserves the baseline");
    Access::workflow(ui);
    completeOperation(ui);
    expect(op.phase == Phase::noPass && !op.applied && s.manualState() == before && !s.lastAppliedExplanation,
           "unsuccessful bounded mono search leaves Current and history untouched");
}
