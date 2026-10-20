#include "StandaloneVideo.h"
#include <juce_cryptography/juce_cryptography.h>
#import <AVFoundation/AVFoundation.h>
#include <chrono>

namespace levelscope::standalone::video
{
namespace
{
Time fromCM (CMTime time)
{
    if (!CMTIME_IS_NUMERIC (time) || time.epoch != 0 || time.timescale <= 0)
        throw std::runtime_error ("Unsupported or invalid media timestamp");
    return Time (time.value, time.timescale);
}
CMTime toCM (Time time)
{
    if (time.denominator > INT32_MAX) throw std::runtime_error ("Media timebase exceeds Core Media limits");
    return CMTimeMake (time.numerator, static_cast<int32_t> (time.denominator));
}
[[noreturn]] void fail (NSError* error, const char* fallback)
{
    throw std::runtime_error (error ? error.localizedDescription.UTF8String : fallback);
}
struct Sample
{
    CMSampleBufferRef value = nullptr;
    explicit Sample (CMSampleBufferRef s) : value (s) {}
    ~Sample() { if (value) CFRelease (value); }
};
class HashInput final : public juce::InputStream
{
public:
    HashInput (juce::InputStream& in, const Cancelled& fn) : input (in), cancelled (fn) {}
    juce::int64 getTotalLength() override { return input.getTotalLength(); }
    juce::int64 getPosition() override { return input.getPosition(); }
    bool setPosition (juce::int64 n) override { return input.setPosition (n); }
    bool isExhausted() override { return cancelled() || input.isExhausted(); }
    int read (void* data, int n) override { return cancelled() ? 0 : input.read (data, std::min (n, 65536)); }
private:
    juce::InputStream& input;
    const Cancelled& cancelled;
};


bool isProRes (CMVideoCodecType codec)
{
    return codec == kCMVideoCodecType_AppleProRes422Proxy
        || codec == kCMVideoCodecType_AppleProRes422LT
        || codec == kCMVideoCodecType_AppleProRes422;
}
bool is709 (CFTypeRef value)
{
    return value && CFGetTypeID (value) == CFStringGetTypeID()
        && CFEqual (value, kCVImageBufferColorPrimaries_ITU_R_709_2);
}
// ProRes uses an explicit, conservative interpretation: all three effective
// format-description colour tags must identify Rec.709. Missing tags are not
// inferred from the codec, raster, filename, or bit depth.
juce::String proResRejection (CMFormatDescriptionRef description)
{
    // AVFoundation can supply default effective tags when the sample entry has
    // none. Require the actual nclc/nclx declaration as well, so missing container
    // tags cannot silently become an assumed Rec.709 interpretation.
    const auto raw = CMFormatDescriptionGetExtension (description, kCMFormatDescriptionExtension_VerbatimSampleDescription);
    bool declared709 = false;
    if (raw && CFGetTypeID (raw) == CFDataGetTypeID())
    {
        const auto data = (CFDataRef) raw;
        const auto* bytes = CFDataGetBytePtr (data);
        const auto length = static_cast<size_t> (CFDataGetLength (data));
        auto u32 = [&] (size_t p) { return (uint32_t (bytes[p]) << 24) | (uint32_t (bytes[p+1]) << 16) | (uint32_t (bytes[p+2]) << 8) | bytes[p+3]; };
        int colourAtoms = 0;
        for (size_t p = 86; p + 8 <= length;)
        {
            const auto size = u32 (p);
            if (size < 8 || size > length - p) break;
            if (!std::memcmp (bytes + p + 4, "colr", 4))
            {
                ++colourAtoms;
                declared709 = size >= 18 && (!std::memcmp (bytes + p + 8, "nclc", 4)
                    || (size >= 19 && !std::memcmp (bytes + p + 8, "nclx", 4) && !(bytes[p+18] & 0x80)))
                    && bytes[p+12] == 0 && bytes[p+13] == 1
                    && bytes[p+14] == 0 && bytes[p+15] == 1
                    && bytes[p+16] == 0 && bytes[p+17] == 1;
            }
            p += size;
        }
        declared709 = declared709 && colourAtoms == 1;
    }
    if (!declared709)
        return "ProRes requires one explicit Rec.709 SDR colr declaration (primaries, transfer and matrix 1/1/1); missing or contradictory colour metadata is unsupported";
    for (const auto key : { kCMFormatDescriptionExtension_ColorPrimaries,
                            kCMFormatDescriptionExtension_TransferFunction,
                            kCMFormatDescriptionExtension_YCbCrMatrix })
        if (!is709 (CMFormatDescriptionGetExtension (description, key)))
            return "ProRes requires explicit, consistent Rec.709 primaries, transfer and matrix tags; missing, HDR or other colour tags are unsupported";
    const auto field = CMFormatDescriptionGetExtension (description, kCMFormatDescriptionExtension_FieldCount);
    int count = 0;
    if (!field || CFGetTypeID (field) != CFNumberGetTypeID()
        || !CFNumberGetValue ((CFNumberRef) field, kCFNumberIntType, &count) || count != 1)
        return "ProRes requires explicit progressive (one-field) metadata; interlaced or unspecified field structure is unsupported";
    const auto dimensions = CMVideoFormatDescriptionGetDimensions (description);
    const auto aperture = CMVideoFormatDescriptionGetCleanAperture (description, false);
    const auto presentation = CMVideoFormatDescriptionGetPresentationDimensions (description, true, true);
    if (std::abs (aperture.size.width - dimensions.width) > 0.01
        || std::abs (aperture.size.height - dimensions.height) > 0.01
        || std::abs (aperture.origin.x) > 0.01 || std::abs (aperture.origin.y) > 0.01
        || std::abs (presentation.width - dimensions.width) > 0.01
        || std::abs (presentation.height - dimensions.height) > 0.01)
        return "ProRes requires square pixels and full, centred aperture";
    return {};
}
double elapsedMs (std::chrono::steady_clock::time_point start)
{
    return std::chrono::duration<double, std::milli> (std::chrono::steady_clock::now() - start).count();
}
class AppleDecoder final : public Decoder
{
public:
    ~AppleDecoder() override { @autoreleasepool { [reader cancelReading]; } }
    std::shared_ptr<const Asset> load (const juce::File& file, int streamId, const Cancelled& cancelled) override
    {
        @autoreleasepool
        {
            if (!file.existsAsFile()) throw std::runtime_error ("Local video file is missing");
            auto result = std::make_shared<Asset>(); result->file = file;
            result->bytes = file.getSize(); result->modifiedMs = file.getLastModificationTime().toMilliseconds();
            if (!file.hasFileExtension ("mov;mp4;m4v")) throw std::runtime_error ("Supported containers: MOV with SDR H.264 or Rec.709 ProRes Proxy/LT/422; MP4/M4V with SDR H.264");
            const auto hashStart = std::chrono::steady_clock::now();
            juce::FileInputStream input (file);
            // SHA256 requests 64-byte blocks; amortize filesystem reads without
            // changing the complete-file hash or cancellation between blocks.
            juce::BufferedInputStream buffered (input, 65536);
            HashInput hashInput (buffered, cancelled);
            result->sha256 = juce::SHA256 (hashInput).toHexString();
            result->hashMilliseconds = elapsedMs (hashStart);
            if (cancelled()) throw std::runtime_error ("Video load cancelled");
            if (input.getStatus().failed() || input.getPosition() != result->bytes) throw std::runtime_error ("Cannot hash the complete video file");
            nativeAsset = [AVURLAsset URLAssetWithURL:[NSURL fileURLWithPath:[NSString stringWithUTF8String:file.getFullPathName().toRawUTF8()]]
                options:@{ AVURLAssetPreferPreciseDurationAndTimingKey: @YES }];
            if (nativeAsset.hasProtectedContent) throw std::runtime_error ("Protected video is unsupported");
            result->embeddedAudioTracks = static_cast<int> ([nativeAsset tracksWithMediaType:AVMediaTypeAudio].count);
            auto tracks = [nativeAsset tracksWithMediaType:AVMediaTypeVideo];
            if (tracks.count == 0) throw std::runtime_error ("No readable video stream");
            for (AVAssetTrack* candidate in tracks)
            {
                Stream info; info.id = candidate.trackID; info.timeScale = candidate.naturalTimeScale;
                info.width = static_cast<int> (candidate.naturalSize.width); info.height = static_cast<int> (candidate.naturalSize.height);
                if (CMTIME_IS_NUMERIC (candidate.minFrameDuration) && candidate.minFrameDuration.value > 0)
                    info.nominalRate = Time (candidate.minFrameDuration.timescale, candidate.minFrameDuration.value);
                bool acceptable = candidate.formatDescriptions.count == 1;
                for (id format in candidate.formatDescriptions)
                {
                    auto description = (__bridge CMFormatDescriptionRef) format;
                    const auto codec = CMFormatDescriptionGetMediaSubType (description);
                    char code[] { char (codec >> 24), char (codec >> 16), char (codec >> 8), char (codec), 0 };
                    info.codec = juce::String (code);
                    const auto presentation = CMVideoFormatDescriptionGetPresentationDimensions (description, true, true);
                    const auto dimensions = CMVideoFormatDescriptionGetDimensions (description);
                    auto transfer = CMFormatDescriptionGetExtension (description, kCMFormatDescriptionExtension_TransferFunction);
                    const bool hdr = transfer && (CFEqual (transfer, CFSTR ("SMPTE_ST_2084_PQ")) || CFEqual (transfer, CFSTR ("ITU_R_2100_HLG")));
                    const bool prores = isProRes (codec);
                    if (prores)
                    {
                        info.codec = codec == kCMVideoCodecType_AppleProRes422Proxy ? "ProRes 422 Proxy"
                                   : codec == kCMVideoCodecType_AppleProRes422LT ? "ProRes 422 LT" : "ProRes 422";
                        info.unsupportedReason = file.hasFileExtension ("mov") ? proResRejection (description)
                            : "ProRes reference picture requires a MOV container";
                    }
                    else if (codec != kCMVideoCodecType_H264)
                        info.unsupportedReason = "Supported codecs: SDR H.264 and Rec.709 ProRes 422 Proxy, LT or standard 422";
                    acceptable = acceptable && (codec == kCMVideoCodecType_H264 || prores)
                        && info.unsupportedReason.isEmpty() && !hdr
                        && std::abs (presentation.width - dimensions.width) < 0.01
                        && std::abs (presentation.height - dimensions.height) < 0.01;
                }
                info.supported = acceptable && info.width > 0 && info.height > 0 && info.width <= 3840 && info.height <= 2160;
                result->streams.push_back (info);
            }
            if (streamId == 0)
                for (const auto& info : result->streams) if (info.supported) { streamId = info.id; break; }
            const auto selected = std::find_if (result->streams.begin(), result->streams.end(), [streamId] (const Stream& s) { return s.id == streamId; });
            if (selected == result->streams.end() || !selected->supported)
            {
                const auto& rejected = selected != result->streams.end() ? *selected : result->streams.front();
                throw std::runtime_error (rejected.unsupportedReason.isNotEmpty() ? rejected.unsupportedReason.toStdString()
                    : "Unsupported stream: use square-pixel full-aperture SDR H.264 or progressive Rec.709 ProRes Proxy/LT/422 up to 3840x2160");
            }
            selectedProRes = selected->codec.startsWith ("ProRes");
            track = [nativeAsset trackWithTrackID:streamId]; result->selectedStream = streamId;
            const auto t = track.preferredTransform;
            transform = juce::AffineTransform (static_cast<float> (t.a), static_cast<float> (t.c), static_cast<float> (t.tx),
                                               static_cast<float> (t.b), static_cast<float> (t.d), static_cast<float> (t.ty));
            if (!std::isfinite (t.a + t.b + t.c + t.d + t.tx + t.ty) || std::abs (t.a * t.d - t.b * t.c) < 0.5)
                throw std::runtime_error ("Unsupported video transform");
            const auto indexStart = std::chrono::steady_clock::now();
            NSError* error = nil;
            auto indexReader = [[AVAssetReader alloc] initWithAsset:nativeAsset error:&error];
            auto indexOutput = [[AVAssetReaderTrackOutput alloc] initWithTrack:track outputSettings:nil];
            if (!indexReader || !indexOutput || ![indexReader canAddOutput:indexOutput]) fail (error, "Cannot index the video stream");
            [indexReader addOutput:indexOutput];
            if (![indexReader startReading]) fail (indexReader.error, "Cannot start the video index");
            while (!cancelled())
            {
                Sample sample ([indexOutput copyNextSampleBuffer]);
                if (!sample.value) break;
                // AVAssetReader may batch compressed frames and emit empty markers.
                // Output timing includes container edit-list shifts and trims. Raw
                // compressed PTS can differ from the decoded presentation timeline.
                const auto samples = CMSampleBufferGetNumSamples (sample.value);
                if (samples == 0) continue;
                if (samples < 0 || static_cast<uint64_t> (samples) > 2000000 - result->frames.size())
                    throw std::runtime_error ("Reference picture exceeds the two-million-frame index limit");
                CMSampleTimingInfo single {};
                CMItemCount entries = 0;
                std::vector<CMSampleTimingInfo> timings;
                auto status = CMSampleBufferGetOutputSampleTimingInfoArray (sample.value, 1, &single, &entries);
                if (status == kCMSampleBufferError_ArrayTooSmall && entries == samples)
                {
                    timings.resize (static_cast<size_t> (entries));
                    status = CMSampleBufferGetOutputSampleTimingInfoArray (sample.value, entries, timings.data(), &entries);
                }
                if (status || (entries != 1 && entries != samples))
                    throw std::runtime_error ("Cannot read per-frame video output timing");
                for (CMItemIndex i = 0; i < samples; ++i)
                {
                    const auto& timing = timings.empty() ? single : timings[static_cast<size_t> (i)];
                    auto pts = fromCM (timing.presentationTimeStamp);
                    const auto duration = timing.duration;
                    if (entries == 1 && i != 0)
                        pts = pts + fromCM (duration) * Time (i, 1);
                    const auto end = CMTIME_IS_NUMERIC (duration) && duration.value > 0 ? pts + fromCM (duration) : pts;
                    result->frames.push_back ({ pts, end });
                }
            }
            if (cancelled()) { [indexReader cancelReading]; throw std::runtime_error ("Video load cancelled"); }
            if (indexReader.status != AVAssetReaderStatusCompleted) fail (indexReader.error, "Video index failed");
            auto& frames = result->frames;
            if (frames.empty()) throw std::runtime_error ("Empty video stream");
            std::sort (frames.begin(), frames.end(), [] (const FrameTime& a, const FrameTime& b) { return a.pts < b.pts; });
            const auto trackEnd = fromCM (CMTimeRangeGetEnd (track.timeRange));
            for (size_t i = 0; i < frames.size(); ++i)
            {
                if (i && !(frames[i - 1].pts < frames[i].pts)) throw std::runtime_error ("Duplicate or non-monotonic video PTS");
                const auto next = i + 1 < frames.size() ? frames[i + 1].pts : trackEnd;
                if (!(frames[i].pts < frames[i].end) || next < frames[i].end) frames[i].end = next;
                if (!(frames[i].pts < frames[i].end)) throw std::runtime_error ("Invalid video frame duration");
                if (i && frames[i].end - frames[i].pts != frames[0].end - frames[0].pts) result->variableRate = true;
            }
            readTimecode (*result, cancelled);
            result->indexMilliseconds = elapsedMs (indexStart);
            if (!result->unchanged()) throw std::runtime_error ("Video file changed during indexing");
            asset = result;
            return result;
        }
    }
    void seek (size_t index) override
    {
        @autoreleasepool
        {
            [reader cancelReading];
            NSError* error = nil;
            reader = [[AVAssetReader alloc] initWithAsset:nativeAsset error:&error];
            output = [[AVAssetReaderTrackOutput alloc] initWithTrack:track outputSettings:
                @{ (id) kCVPixelBufferPixelFormatTypeKey: @(kCVPixelFormatType_32BGRA) }];
            output.alwaysCopiesSampleData = NO;
            if (!reader || !output || ![reader canAddOutput:output]) fail (error, "Cannot create video decoder");
            [reader addOutput:output];
            // Arbitrary times clip the first PTS. Always seek to an indexed source PTS.
            reader.timeRange = CMTimeRangeMake (toCM (asset->frames.at (index).pts), kCMTimePositiveInfinity);
            if (![reader startReading]) fail (reader.error, "Video decoder could not start");
        }
    }
    std::optional<Picture> next() override
    {
        @autoreleasepool
        {
            Sample sample ([output copyNextSampleBuffer]);
            if (!sample.value)
            {
                if (reader.status != AVAssetReaderStatusCompleted) fail (reader.error, "Video decode failed");
                return {};
            }
            const auto pts = fromCM (CMSampleBufferGetOutputPresentationTimeStamp (sample.value));
            const auto index = frameAt (asset->frames, pts);
            if (!index || asset->frames[*index].pts != pts) throw std::runtime_error ("Decoded PTS is absent from the source index");
            auto pixel = CMSampleBufferGetImageBuffer (sample.value);
            if (!pixel || CVPixelBufferGetPixelFormatType (pixel) != kCVPixelFormatType_32BGRA) throw std::runtime_error ("Unsupported decoded pixel format");
            if (selectedProRes)
                for (const auto key : { kCVImageBufferColorPrimariesKey, kCVImageBufferTransferFunctionKey, kCVImageBufferYCbCrMatrixKey })
                {
                    const auto tag = CVBufferGetAttachment (pixel, key, nullptr);
                    if (tag && !is709 (tag))
                        throw std::runtime_error ("Decoded ProRes colour tags contradict the required Rec.709 stream metadata");
                }
            const auto width = static_cast<int> (CVPixelBufferGetWidth (pixel)), height = static_cast<int> (CVPixelBufferGetHeight (pixel));
            if (width <= 0 || height <= 0 || width > 3840 || height > 2160) throw std::runtime_error ("Decoded picture exceeds size limit");
            // CoreGraphics copies a SoftwareImage into a new CFData on every
            // paint. Native images retain immutable backing storage and cache
            // their CGImage instead, avoiding repeated 1080p presentation copies.
            // Each decoded picture still owns separate pixels; never reuse or
            // mutate a published image while the UI may retain it.
            auto makePictureImage = [this] (int w, int h)
            {
                if (selectedProRes)
                    return juce::Image (juce::Image::ARGB, w, h, true, juce::NativeImageType());
                return juce::Image (juce::Image::ARGB, w, h, true, juce::SoftwareImageType());
            };
            juce::Image image = makePictureImage (width, height);
            if (CVPixelBufferLockBaseAddress (pixel, kCVPixelBufferLock_ReadOnly) != kCVReturnSuccess) throw std::runtime_error ("Cannot access decoded picture");
            {
                juce::Image::BitmapData target (image, juce::Image::BitmapData::writeOnly);
                const auto* base = static_cast<const uint8_t*> (CVPixelBufferGetBaseAddress (pixel));
                const auto stride = CVPixelBufferGetBytesPerRow (pixel);
                for (int y = 0; y < height; ++y) std::memcpy (target.getLinePointer (y), base + y * stride, static_cast<size_t> (width) * 4);
            }
            CVPixelBufferUnlockBaseAddress (pixel, kCVPixelBufferLock_ReadOnly);
            if (!transform.isIdentity())
            {
                const auto bounds = juce::Rectangle<float> (0, 0, static_cast<float> (width), static_cast<float> (height)).transformedBy (transform).getSmallestIntegerContainer();
                if (bounds.getWidth() > 3840 || bounds.getHeight() > 3840 || bounds.isEmpty()) throw std::runtime_error ("Transformed picture exceeds size limit");
                juce::Image rotated = makePictureImage (bounds.getWidth(), bounds.getHeight());
                juce::Graphics g (rotated); g.drawImageTransformed (image, transform.translated (-static_cast<float> (bounds.getX()), -static_cast<float> (bounds.getY())));
                image = rotated;
            }
            return Picture { std::move (image), *index, 0, 0 };
        }
    }
private:
    void readTimecode (Asset& result, const Cancelled& cancelled)
    {
        auto tracks = [track associatedTracksOfType:AVTrackAssociationTypeTimecode];
        if (tracks.count == 0) return;
        result.timecodeStatus = "Unsupported timecode metadata";
        if (tracks.count != 1) return;
        NSError* error = nil;
        auto tcReader = [[AVAssetReader alloc] initWithAsset:nativeAsset error:&error];
        auto tcOutput = [[AVAssetReaderTrackOutput alloc] initWithTrack:tracks[0] outputSettings:nil];
        if (!tcReader || !tcOutput || ![tcReader canAddOutput:tcOutput]) return;
        [tcReader addOutput:tcOutput]; if (![tcReader startReading]) return;
        std::optional<Timecode> parsed;
        while (!cancelled())
        {
            Sample sample ([tcOutput copyNextSampleBuffer]);
            if (!sample.value) break;
            const auto count = CMSampleBufferGetNumSamples (sample.value);
            if (count == 0) continue; // Empty reader markers have neither format nor data.
            if (count != 1 || parsed) { result.timecodeStatus = "Multiple timecode origins are unsupported"; return; }
            auto format = CMSampleBufferGetFormatDescription (sample.value);
            auto buffer = CMSampleBufferGetDataBuffer (sample.value);
            if (!format || !buffer) return;
            const auto type = CMFormatDescriptionGetMediaSubType (format);
            int64_t origin = 0;
            if (type == kCMTimeCodeFormatType_TimeCode32)
            {
                uint32_t v = 0; if (CMBlockBufferCopyDataBytes (buffer, 0, 4, &v)) return;
                origin = static_cast<int32_t> (CFSwapInt32BigToHost (v));
            }
            else if (type == kCMTimeCodeFormatType_TimeCode64)
            {
                uint64_t v = 0; if (CMBlockBufferCopyDataBytes (buffer, 0, 8, &v)) return;
                origin = static_cast<int64_t> (CFSwapInt64BigToHost (v));
            }
            else return;
            const auto duration = CMTimeCodeFormatDescriptionGetFrameDuration (format);
            if (!CMTIME_IS_NUMERIC (duration) || duration.value <= 0) return;
            const auto flags = CMTimeCodeFormatDescriptionGetTimeCodeFlags (format);
            parsed = Timecode { origin, fromCM (CMSampleBufferGetOutputPresentationTimeStamp (sample.value)),
                Time (duration.timescale, duration.value), (flags & kCMTimeCodeFlag_DropFrame) != 0,
                (flags & kCMTimeCodeFlag_24HourMax) != 0 };
        }
        if (cancelled()) { [tcReader cancelReading]; throw std::runtime_error ("Video load cancelled"); }
        if (tcReader.status != AVAssetReaderStatusCompleted || !parsed) return;
        const auto tc = *parsed;
        const auto origin = tc.startFrame;
        try { (void) timecodeLabel (origin, tc.rate, tc.dropFrame, tc.wraps24Hours); }
        catch (...) { return; }
        result.timecode = tc;
        result.timecodeStatus = "Origin " + juce::String (timecodeLabel (origin, tc.rate, tc.dropFrame, tc.wraps24Hours))
                             + (tc.dropFrame ? " (drop-frame)" : " (non-drop)");
    }
    AVURLAsset* nativeAsset = nil;
    AVAssetTrack* track = nil;
    AVAssetReader* reader = nil;
    AVAssetReaderTrackOutput* output = nil;
    std::shared_ptr<const Asset> asset;
    juce::AffineTransform transform;
    bool selectedProRes = false;
};
}
std::unique_ptr<Decoder> makeDecoder() { return std::make_unique<AppleDecoder>(); }
}
