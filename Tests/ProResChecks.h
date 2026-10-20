// Included inside VideoTests.cpp's test namespace; never in the application.
void proResChecks (const juce::File& file, int64_t count, int64_t scale, int64_t step, bool generated, bool vfr)
{
    beginVideoReaderAudit();
    auto decoder = makeDecoder();
    const auto asset = decoder->load (file, 0, [] { return false; });
    expect (asset->frames.size() == static_cast<size_t> (count), "ProRes exact frame count");
    expect (asset->streams.front().codec.startsWith ("ProRes"), "ProRes recognized by exact profile");
    expect (asset->streams.front().width == 1920 && asset->streams.front().height == 1080, "required 1080 raster");
    expect (!asset->sha256.isEmpty(), "source hash recorded");
    if (!generated)
        expect (asset->timecodeStatus == "Origin 00:00:00:00 (non-drop)", "actual movie timecode origin preserved");
    auto verifyBurn = [&] (const Picture& picture, size_t index)
    {
        if (!generated) return;
        int burn = 0;
        for (int bit = 0; bit < 8; ++bit)
            if (picture.image.getPixelAt ((bit * 2 + 1) * picture.image.getWidth() / 16,
                                         picture.image.getHeight() * 7 / 8).getBrightness() > 0.5f) burn |= 1 << bit;
        expect (burn == static_cast<int> (index & 255), "ProRes independent burn-in matches indexed frame");
    };
    std::vector<size_t> samples { 0, 1, static_cast<size_t> (count / 2), static_cast<size_t> (count - 1), 20, 2 };
    std::vector<juce::Image> reference;
    decoder->seek (0);
    for (int64_t i = 0; i < count; ++i)
    {
        const auto ticks = i + (vfr ? i / 3 : 0);
        expect (asset->frames[static_cast<size_t> (i)].pts == Time (ticks * step, scale), "ProRes rational PTS matches independent fixture timeline");
        const auto endTicks = vfr && i == count - 1 ? ticks + 1 : i + 1 + (vfr ? (i + 1) / 3 : 0);
        expect (asset->frames[static_cast<size_t> (i)].end == Time (endTicks * step, scale), "ProRes exact frame coverage and VFR gaps");
        if (count > 3000) continue; // Long-file test indexes all samples; decode representative samples below.
        auto picture = decoder->next();
        expect (picture && picture->index == static_cast<size_t> (i), "ProRes sequential decode index");
        verifyBurn (*picture, static_cast<size_t> (i));
        if (std::find (samples.begin(), samples.end(), static_cast<size_t> (i)) != samples.end())
            reference.push_back (picture->image);
        if (generated && i == 0)
        {
            const int y = picture->image.getHeight() * 5 / 12, w = picture->image.getWidth();
            const auto red = picture->image.getPixelAt (w / 12, y);
            const auto green = picture->image.getPixelAt (w * 3 / 12, y);
            const auto blue = picture->image.getPixelAt (w * 5 / 12, y);
            expect (red.getRed() > 210 && red.getGreen() < 40 && red.getBlue() < 40, "Rec709 red conversion tolerance");
            expect (green.getGreen() > 210 && green.getRed() < 40 && green.getBlue() < 40, "Rec709 green conversion tolerance");
            expect (blue.getBlue() > 210 && blue.getRed() < 40 && blue.getGreen() < 40, "Rec709 blue conversion tolerance");
            const auto grey = picture->image.getPixelAt (w * 11 / 12, y);
            expect (std::abs (int (grey.getRed()) - 128) <= 12 && std::abs (int (grey.getGreen()) - 128) <= 12
                    && std::abs (int (grey.getBlue()) - 128) <= 12, "Rec709 grey conversion tolerance");
        }
    }
    for (const auto index : samples)
    {
        decoder->seek (index); const auto picture = decoder->next();
        expect (picture && picture->index == index, "ProRes random seek selects exact PTS");
        verifyBurn (*picture, index);
        // Real user material has no generated binary burn-in. Compare pixels of
        // independent sequential and random access instead of inventing its content.
        if (!generated)
        {
            auto ordered = samples; std::sort (ordered.begin(), ordered.end());
            const auto position = static_cast<size_t> (std::find (ordered.begin(), ordered.end(), index) - ordered.begin());
            const auto& expected = reference.at (position);
            for (int y = 10; y < expected.getHeight(); y += 97)
                for (int x = 10; x < expected.getWidth(); x += 101)
                    expect (expected.getPixelAt (x, y) == picture->image.getPixelAt (x, y), "actual movie random/sequential pixel agreement");
        }
    }
    Player player; player.attach (file);
    uint64_t generation = 0;
    for (const auto index : samples)
    {
        const auto waitSeconds = static_cast<int> (30 + asset->bytes / (64 * 1024 * 1024) + count / 10000);
        auto view = awaitPicture (player, asset->frames[index].pts, ++generation, waitSeconds);
        expect (view.picture->index == index && view.picture->generation == generation, "ProRes worker seek generation");
        expect (view.queued <= 8, "ProRes bounded queue");
    }
    player.clear(); expect (!player.update ({}, ++generation).picture, "ProRes clear removes picture");
    expect (videoReaderAuditCount (0) > 0 && videoReaderAuditCount (1) == 0, "observed AVAssetReader outputs are video only");
    std::cout << "PRORES " << file.getFileName() << " sha256=" << asset->sha256 << " frames=" << count
              << " hash_ms=" << asset->hashMilliseconds << " index_ms=" << asset->indexMilliseconds
              << " observed_video_outputs=" << videoReaderAuditCount (0) << " observed_audio_outputs=" << videoReaderAuditCount (1)
              << " timecode=" << asset->timecodeStatus << '\n';
    // Audit stays installed for remaining reader destruction; counters reset per test.
}

void isolationChecks (const juce::File& movie)
{
    auto selectedFile = juce::File::getSpecialLocation (juce::File::tempDirectory).getNonexistentChildFile ("prores-selected-nonzero", ".wav", false);
    struct Remove { juce::File f; ~Remove() { f.deleteFile(); } } remove { selectedFile };
    {
        juce::WavAudioFormat wav; auto stream = selectedFile.createOutputStream();
        std::unique_ptr<juce::AudioFormatWriter> writer (wav.createWriterFor (stream.get(), 48000, 2, 24, {}, 0));
        expect (writer != nullptr, "isolation source writer"); stream.release();
        juce::AudioBuffer<float> signal (2, 48000);
        for (int n=0;n<48000;++n) {
            signal.setSample (0,n,0.25f * std::sin (float (n) * 2.0f * juce::MathConstants<float>::pi * 337.0f / 48000.0f));
            signal.setSample (1,n,0.125f * std::sin (float (n) * 2.0f * juce::MathConstants<float>::pi * 613.0f / 48000.0f));
        }
        expect (writer->writeFromAudioSampleBuffer (signal,0,48000), "nonzero selected fixture");
    }
    JobControl job; auto validation=validateMedia (selectedFile,job); expect (validation.succeeded(),"isolation source valid");
    auto media=validation.media; auto layout=media.suggestedLayout;layout.roles[0]=juce::AudioChannelSet::left;layout.roles[1]=juce::AudioChannelSet::right;layout.confirmed=true;
    OutputRouting routing;routing.fileToDevice[0]=0;routing.fileToDevice[1]=1;
    auto device=std::make_unique<FakeDevice>();auto* fake=device.get();
    StandaloneAudition audio (std::move(device));audio.sourceChanged(media);
    expect (audio.connect(48000).isEmpty(),"controlled sink connected");
    auto capture=[&] {
        expect (audio.play(selectedFile,media,layout,routing,true).isEmpty(),"selected audio restart");
        for(int n=0;n<2000 && !fake->isPlaying();++n){audio.poll();std::this_thread::sleep_for(std::chrono::milliseconds(1));}
        expect(fake->isPlaying(),"controlled sink started");
        juce::AudioBuffer<float> block(2,256);fake->tick(block);audio.stop();
        expect(block.getMagnitude(0,256)>1e-6f,"captured signal is independently nonzero before physical output");
        return block;
    };
    auto baseline=capture();beginVideoReaderAudit();
    const auto originalRouting=routing;const auto originalHash=media.sha256;
    Player player;player.attach(movie);auto view=awaitPicture(player,Time(0,1),1);
    expect(view.asset->embeddedAudioTracks>0,"isolation movie actually contains embedded audio");
    for(int pass=0;pass<4;++pass) {
        if(pass==1)(void)awaitPicture(player,Time(1,1),2);
        if(pass==2){
            const auto stream = view.asset->streams.back().id;
            player.attach(movie, stream);
            const auto switched = awaitPicture(player,Time(0,1),3);
            expect(switched.asset->selectedStream==stream, "explicit ProRes stream selection preserved");
            expect(switched.asset->sha256==view.asset->sha256, "stream switch preserves movie identity");
        }
        if(pass==3)player.clear();
        auto withVideo=capture();
        for(int c=0;c<2;++c)for(int n=0;n<256;++n)
            expect(withVideo.getSample(c,n)==baseline.getSample(c,n),"embedded signal absent: captured selected samples bit-identical");
        expect(media.sha256==originalHash && routing.fileToDevice==originalRouting.fileToDevice,"source identity and routing unchanged");
    }
    expect(videoReaderAuditCount(0)>0 && videoReaderAuditCount(1)==0,"runtime audit: no embedded audio reader output");
    independence();
    std::cout<<"ISOLATION PASS nonzero_selected_peak="<<baseline.getMagnitude(0,256)<<" exact_compared_samples=2048 physical_outputs=0 observed_audio_reader_outputs="<<videoReaderAuditCount(1)<<'\n';
}

void loadingChecks (const juce::File& file)
{
    const auto bytes=file.getSize();
    expect(bytes>8LL*1024*1024*1024,"realistically large ProRes file exceeds 8 GiB");
    auto now=[] {return std::chrono::steady_clock::now();};
    auto seconds=[&](auto start){return std::chrono::duration<double>(now()-start).count();};
    const auto start=now();auto decoder=makeDecoder();
    auto asset=decoder->load(file,0,[]{return false;});
    const auto loaded=seconds(start);
    decoder->seek(0);const auto picture=decoder->next();
    expect(picture && picture->index==0,"large file first decoded picture");
    const auto firstPicture=seconds(start);
    const double budget=30.0+double(file.getSize())/(64*1024*1024)+10.0*asset->frames.size()/100000.0;
    expect(firstPicture<=budget,"large attachment within predeclared size-scaled budget");
    const auto hashMs=asset->hashMilliseconds,indexMs=asset->indexMilliseconds;
    decoder.reset();asset.reset();
    // Asynchronous user cancellation while hashing a large file.
    std::atomic<bool> cancel{false};bool rejected=false;
    std::thread worker([&]{
        try{auto d=makeDecoder();d->load(file,0,[&]{return cancel.load();});}
        catch(const std::exception&){rejected=true;}
    });
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    auto cancelStart=now();cancel.store(true);worker.join();const auto hashCancel=seconds(cancelStart);
    expect(rejected && hashCancel<=2.0,"hash cancellation completes within two seconds");
    // SHA256 consumes 64-byte blocks. Request cancellation only after it has
    // completed, in the index loop; this test deliberately does a second hash.
    int64_t polls=0;bool requested=false;auto indexRequested=now();rejected=false;
    try{
        auto d=makeDecoder();
        d->load(file,0,[&]{
            if(++polls > bytes/64+30){
                if(!requested){requested=true;indexRequested=now();}
                return true;
            }
            return false;
        });
    }catch(const std::exception&){rejected=true;}
    const auto indexCancel=seconds(indexRequested);
    expect(requested && rejected && indexCancel<=2.0,"index cancellation completes within two seconds");
    const auto clearStart=now();
    {
        Player player;player.attach(file);
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        player.clear();expect(player.update({},1).state==Player::State::empty,"clear removes loading reference immediately");
    }
    const auto cleared=seconds(clearStart)-0.1;
    expect(cleared<=2.0,"clear and destruction stop large-file worker promptly");
    std::cout<<"LOADING bytes="<<file.getSize()<<" hash_ms="<<hashMs<<" index_ms="<<indexMs
             <<" loaded_s="<<loaded<<" first_picture_s="<<firstPicture<<" budget_s="<<budget
             <<" hash_cancel_s="<<hashCancel<<" index_cancel_s="<<indexCancel<<" clear_join_s="<<cleared
             <<" cache_state=uncontrolled_OS_cache\n";
}
