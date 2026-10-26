#pragma once
// Independent test decoding: never calls the production reader/quantizer for its oracle.
namespace
{
uint64_t little(const uint8_t* p,int bytes){uint64_t n=0;for(int i=0;i<bytes;++i)n|=uint64_t(p[i])<<(i*8);return n;}
int32_t pcmCode(const uint8_t* p){const auto n=static_cast<uint32_t>(little(p,3));return static_cast<int32_t>((n&0x800000)?(n|0xff000000):n);}
void rawEncoderMatrix(const juce::File& root)
{
    juce::Array<juce::var> manifest;
    for(int channels:{2,6})for(int rate:{44100,48000,96000,192000})for(double gain:{0.,6.020599913279624,-6.020599913279624,-3.305074118944,-60.,24.})
    {
        const int frames=1025;juce::AudioBuffer<float> audio(channels,frames);
        for(int ch=0;ch<channels;++ch)for(int n=0;n<frames;++n)
            audio.setSample(ch,n,n==0?0.f:float(((n*37+ch*101)%1001-500)*(ch+1))/8388608.f);
        const auto original=audio;StaticGainProcessor processor(gain);
        expect(processor.process(audio.getArrayOfWritePointers(),channels,frames),"raw matrix shared processing");
        bool exact=true;for(int ch=0;ch<channels;++ch)for(int n=0;n<frames;++n)
            exact&=audio.getSample(ch,n)==float(double(original.getSample(ch,n))*std::pow(10.,gain/20.));
        expect(exact,"raw matrix full-precision shared gain oracle");
        pcm::Format format;juce::String error;expect(pcm::describe(rate,channels,frames,format,error),"raw PCM supported rate/layout");
        std::vector<int> map;for(int ch=0;ch<channels;++ch)map.push_back(channels-1-ch);
        const uint64_t seed=UINT64_C(0x1a2b3c4d5e6f7788);const auto base=juce::String(manifest.size());
        JobControl control;
        for(int block:{1,997,4096})
        {
            const auto file=root.getChildFile(base+"-"+juce::String(block)+".wav");pcm::Writer writer(file,format,{},control);
            for(int first=0;first<frames;){const int n=std::min(block,frames-first);std::vector<const float*> ptr;for(int ch=0;ch<channels;++ch)ptr.push_back(audio.getReadPointer(ch,first));writer.write(ptr.data(),map,n,seed,control);first+=n;}
            writer.finish(control);
        }
        const auto file=root.getChildFile(base+"-1.wav");
        expect(pcm::hashFile(file,control)==pcm::hashFile(root.getChildFile(base+"-997.wav"),control)
            &&pcm::hashFile(file,control)==pcm::hashFile(root.getChildFile(base+"-4096.wav"),control),"PCM output exact across block sizes and repeat writes");
        auto out=root.getChildFile(base+".f32").createOutputStream();expect(out!=nullptr,"oracle capture open");
        for(int ch=0;ch<channels;++ch)expect(out->write(audio.getReadPointer(ch),frames*sizeof(float)),"oracle exact pre-quantization floats");out->flush();out.reset();
        auto row=juce::var(new juce::DynamicObject());auto* o=row.getDynamicObject();o->setProperty("file",file.getFileName());o->setProperty("capture",base+".f32");o->setProperty("rate",rate);o->setProperty("channels",channels);o->setProperty("frames",frames);o->setProperty("gain_db",gain);o->setProperty("seed_hex","1a2b3c4d5e6f7788");juce::Array<juce::var> mapping;for(auto ch:map)mapping.add(ch);o->setProperty("map",mapping);manifest.add(row);
        root.getChildFile(base+"-997.wav").deleteFile();root.getChildFile(base+"-4096.wav").deleteFile();
    }
    expect(root.getChildFile("raw-oracles.json").replaceWithText(ex::reportJson(manifest)),"persist independent quantizer oracle matrix");
}
void containerFailures(const juce::File& root)
{
    pcm::Format format;juce::String error;JobControl control;
    const int64_t limit=(UINT64_C(0xfffffffe)-60)/6;
    expect(pcm::describe(192000,2,limit,format,error)&&!format.rf64,"total RIFF size below limit");
    expect(pcm::describe(192000,2,limit+1,format,error)&&format.rf64,"total RIFF size above limit uses RF64 before data overflow");
    auto bad=format;++bad.dataBytes;bool rejected=false;try{pcm::header(bad);}catch(...){rejected=true;}expect(rejected,"inconsistent encoder size rejected");
    // A sparse test file exercises >4 GiB header arithmetic; real full writes are a separate long gate.
    auto file=root.getChildFile("sparse-rf64.wav");auto output=file.createOutputStream();const auto header=pcm::header(format);
    expect(output&&output->write(header.data(),header.size())&&output->setPosition(static_cast<int64_t>(format.fileBytes-1))&&output->writeByte(0),"virtual large container fixture");output->flush();output.reset();
    pcm::Format parsed;expect(pcm::inspect(file,parsed,error,control)&&parsed.frames==format.frames,"RF64 64-bit length parser");
    for(auto offset:{uint64_t(20),uint64_t(28),uint64_t(36),uint64_t(44),uint64_t(76),uint64_t(80),uint64_t(100)})
    {
        auto stream=file.createOutputStream();expect(stream&&stream->setPosition(offset)&&stream->writeByte(0x7f),"malformed RF64 byte injection");stream.reset();
        expect(!pcm::inspect(file,parsed,error,control),"malformed RF64 size/count/table/mask/GUID/data rejected");
        stream=file.createOutputStream();expect(stream&&stream->setPosition(offset)&&stream->writeByte(static_cast<char>(header[offset])),"restore RF64 test byte");stream.reset();
    }
    file.deleteFile();
}
}
