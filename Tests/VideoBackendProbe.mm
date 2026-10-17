// Standalone compatibility spike. No audio device, AVPlayer or audio output exists.
#import <AVFoundation/AVFoundation.h>
#import <AppKit/AppKit.h>
#include <cstdio>
#include <thread>
#include <chrono>

static bool makeMovie (NSString* path, bool vfr, int frameCount = 120, int timescale = 30000, int step = 1001, bool timecode = false)
{
    NSError* error = nil;
    auto writer = [[AVAssetWriter alloc] initWithURL:[NSURL fileURLWithPath:path]
                                          fileType:AVFileTypeQuickTimeMovie error:&error];
    writer.movieTimeScale = timescale;
    auto input = [AVAssetWriterInput assetWriterInputWithMediaType:AVMediaTypeVideo outputSettings:
        @{ AVVideoCodecKey: AVVideoCodecTypeH264, AVVideoWidthKey: @320, AVVideoHeightKey: @180,
           AVVideoCompressionPropertiesKey: @{ AVVideoAverageBitRateKey: @1000000,
               AVVideoMaxKeyFrameIntervalKey: @30, AVVideoAllowFrameReorderingKey: @YES } }];
    input.mediaTimeScale = timescale;
    AVAssetWriterInput* tcInput = nil;
    if (timecode)
    {
        tcInput = [AVAssetWriterInput assetWriterInputWithMediaType:AVMediaTypeTimecode outputSettings:nil];
        [input addTrackAssociationWithTrackOfInput:tcInput type:AVTrackAssociationTypeTimecode];
        [writer addInput:tcInput];
    }
    auto adaptor = [AVAssetWriterInputPixelBufferAdaptor assetWriterInputPixelBufferAdaptorWithAssetWriterInput:input
        sourcePixelBufferAttributes:@{ (id) kCVPixelBufferPixelFormatTypeKey: @(kCVPixelFormatType_32BGRA),
            (id) kCVPixelBufferWidthKey: @320, (id) kCVPixelBufferHeightKey: @180,
            (id) kCVPixelBufferCGImageCompatibilityKey: @YES, (id) kCVPixelBufferCGBitmapContextCompatibilityKey: @YES }];
    [writer addInput:input];
    if (![writer startWriting]) { fprintf (stderr, "%s\n", writer.error.description.UTF8String); return false; }
    [writer startSessionAtSourceTime:kCMTimeZero];
    if (timecode)
    {
        CMTimeCodeFormatDescriptionRef description = nullptr;
        if (CMTimeCodeFormatDescriptionCreate (nullptr, kCMTimeCodeFormatType_TimeCode32, CMTimeMake (step, timescale),
            30, kCMTimeCodeFlag_DropFrame | kCMTimeCodeFlag_24HourMax, nullptr, &description)) return false;
        CMBlockBufferRef block = nullptr;
        CMBlockBufferCreateWithMemoryBlock (nullptr, nullptr, 4, nullptr, nullptr, 0, 4, 0, &block);
        uint32_t origin = CFSwapInt32HostToBig (107892); // 01:00:00;00 at 30000/1001
        CMBlockBufferReplaceDataBytes (&origin, block, 0, 4);
        CMSampleTimingInfo timing { CMTimeMake (static_cast<int64_t> (frameCount) * step, timescale), kCMTimeZero, kCMTimeInvalid };
        size_t size = 4; CMSampleBufferRef sample = nullptr;
        auto error = CMSampleBufferCreateReady (nullptr, block, description, 1, 1, &timing, 1, &size, &sample);
        if (error || ![tcInput appendSampleBuffer:sample]) return false;
        CFRelease (sample); CFRelease (block); CFRelease (description); [tcInput markAsFinished];
    }
    int64_t pts = 0;
    for (int frame = 0; frame < frameCount; ++frame)
    {
        while (!input.readyForMoreMediaData && writer.status == AVAssetWriterStatusWriting)
            std::this_thread::sleep_for (std::chrono::milliseconds (1));
        CVPixelBufferRef pixel = nullptr;
        if (CVPixelBufferPoolCreatePixelBuffer (nullptr, adaptor.pixelBufferPool, &pixel) != kCVReturnSuccess) return false;
        CVPixelBufferLockBaseAddress (pixel, 0);
        auto colour = CGColorSpaceCreateDeviceRGB();
        auto context = CGBitmapContextCreate (CVPixelBufferGetBaseAddress (pixel), 320, 180, 8,
            CVPixelBufferGetBytesPerRow (pixel), colour, kCGBitmapByteOrder32Little | kCGImageAlphaPremultipliedFirst);
        CGContextSetRGBFillColor (context, frame % 30 == 0 ? 1 : 0.05, 0.1, (frame % 120) / 120.0, 1);
        CGContextFillRect (context, CGRectMake (0, 0, 320, 180));
        CGContextSetRGBFillColor (context, 1, 1, 1, 1);
        CGContextSelectFont (context, "Helvetica", 24, kCGEncodingMacRoman);
        char text[80];
        snprintf (text, sizeof text, "F%05d %lld/%d", frame, (long long) pts, timescale);
        CGContextShowTextAtPoint (context, 12, 85, text, strlen (text));
        // Encode frame index as 8 black/white squares for machine inspection.
        for (int bit = 0; bit < 8; ++bit)
        {
            const auto value = (frame >> bit) & 1;
            CGContextSetRGBFillColor (context, value, value, value, 1);
            CGContextFillRect (context, CGRectMake (10 + 30 * bit, 10, 20, 20));
        }
        CGContextRelease (context); CGColorSpaceRelease (colour);
        CVPixelBufferUnlockBaseAddress (pixel, 0);
        const bool ok = [adaptor appendPixelBuffer:pixel withPresentationTime:CMTimeMake (pts, timescale)];
        CVPixelBufferRelease (pixel);
        if (!ok) { fprintf (stderr, "%s\n", writer.error.description.UTF8String); return false; }
        pts += vfr && frame % 3 == 2 ? 2 * step : step;
    }
    [input markAsFinished];
    [writer endSessionAtSourceTime:CMTimeMake (pts, timescale)];
    dispatch_semaphore_t done = dispatch_semaphore_create (0);
    [writer finishWritingWithCompletionHandler:^{ dispatch_semaphore_signal (done); }];
    dispatch_semaphore_wait (done, DISPATCH_TIME_FOREVER);
    return writer.status == AVAssetWriterStatusCompleted;
}

static bool probe (NSString* path, bool present)
{
    auto asset = [AVURLAsset URLAssetWithURL:[NSURL fileURLWithPath:path]
                                   options:@{ AVURLAssetPreferPreciseDurationAndTimingKey: @YES }];
    auto tracks = [asset tracksWithMediaType:AVMediaTypeVideo];
    if (tracks.count == 0) return false;
    auto track = tracks[0];
    printf ("video_tracks=%lu embedded_audio_tracks=%lu track_id=%d dimensions=%.0fx%.0f rate=%.9g min_duration=%lld/%d\n",
        (unsigned long) tracks.count, (unsigned long) [asset tracksWithMediaType:AVMediaTypeAudio].count,
        track.trackID, track.naturalSize.width, track.naturalSize.height, track.nominalFrameRate,
        (long long) track.minFrameDuration.value, track.minFrameDuration.timescale);
    NSError* error = nil;
    NSWindow* window = nil;
    if (present)
    {
        [NSApplication sharedApplication];
        [NSApp setActivationPolicy:NSApplicationActivationPolicyRegular];
        window = [[NSWindow alloc] initWithContentRect:NSMakeRect (150, 150, 640, 360)
            styleMask:NSWindowStyleMaskTitled | NSWindowStyleMaskClosable backing:NSBackingStoreBuffered defer:NO];
        window.title = @"LevelScope video compatibility probe — SILENT";
        window.contentView.wantsLayer = YES;
        [window makeKeyAndOrderFront:nil]; [NSApp activateIgnoringOtherApps:YES];
    }
    int count = 0;
    for (int pass = 0; pass < 2; ++pass)
    {
        auto reader = [[AVAssetReader alloc] initWithAsset:asset error:&error];
        auto output = [[AVAssetReaderTrackOutput alloc] initWithTrack:track outputSettings:
            @{ (id) kCVPixelBufferPixelFormatTypeKey: @(kCVPixelFormatType_32BGRA) }];
        output.alwaysCopiesSampleData = NO;
        [reader addOutput:output];
        if (pass) reader.timeRange = CMTimeRangeMake (CMTimeMake (1001 * 45, 30000), kCMTimePositiveInfinity);
        if (![reader startReading]) return false;
        CMTime previous = kCMTimeInvalid;
        int passCount = 0;
        while (auto sample = [output copyNextSampleBuffer])
        {
            auto pts = CMSampleBufferGetPresentationTimeStamp (sample);
            auto duration = CMSampleBufferGetDuration (sample);
            auto pixel = CMSampleBufferGetImageBuffer (sample);
            if (!pixel || (CMTIME_IS_VALID (previous) && CMTimeCompare (previous, pts) >= 0)) return false;
            if (passCount < 3) printf ("pass=%d pts=%lld/%d duration=%lld/%d\n", pass,
                (long long) pts.value, pts.timescale, (long long) duration.value, duration.timescale);
            previous = pts;
            if (present && pass == 0)
            {
                CVPixelBufferLockBaseAddress (pixel, kCVPixelBufferLock_ReadOnly);
                auto colour = CGColorSpaceCreateDeviceRGB();
                auto context = CGBitmapContextCreate (CVPixelBufferGetBaseAddress (pixel),
                    CVPixelBufferGetWidth (pixel), CVPixelBufferGetHeight (pixel), 8,
                    CVPixelBufferGetBytesPerRow (pixel), colour, kCGBitmapByteOrder32Little | kCGImageAlphaPremultipliedFirst);
                auto picture = CGBitmapContextCreateImage (context);
                window.contentView.layer.contents = (__bridge id) picture;
                [window.contentView displayIfNeeded];
                CGImageRelease (picture); CGContextRelease (context); CGColorSpaceRelease (colour);
                CVPixelBufferUnlockBaseAddress (pixel, kCVPixelBufferLock_ReadOnly);
                [[NSRunLoop currentRunLoop] runUntilDate:[NSDate dateWithTimeIntervalSinceNow:0.01]];
            }
            CFRelease (sample); ++passCount; ++count;
        }
        if (reader.status != AVAssetReaderStatusCompleted || passCount == 0)
        { fprintf (stderr, "%s\n", reader.error.description.UTF8String); return false; }
        printf ("pass=%d frames=%d reader_outputs=%lu audio_outputs=0\n", pass, passCount, (unsigned long) reader.outputs.count);
    }
    if (window) { [window orderOut:nil]; printf ("native_layer_presentations=120 (submission only; no physical sync claim)\n"); }
    return count > 120;
}

static bool mux (NSString* movie, NSString* audio, NSString* outputPath)
{
    auto source = [AVURLAsset URLAssetWithURL:[NSURL fileURLWithPath:movie] options:nil];
    auto soundtrack = [AVURLAsset URLAssetWithURL:[NSURL fileURLWithPath:audio] options:nil];
    auto composition = [AVMutableComposition composition];
    NSError* error = nil;
    for (int i = 0; i < 2; ++i)
    {
        auto v = [composition addMutableTrackWithMediaType:AVMediaTypeVideo preferredTrackID:kCMPersistentTrackID_Invalid];
        if (![v insertTimeRange:CMTimeRangeMake (kCMTimeZero, source.duration) ofTrack:[source tracksWithMediaType:AVMediaTypeVideo][0]
                         atTime:kCMTimeZero error:&error]) return false;
    }
    auto a = [composition addMutableTrackWithMediaType:AVMediaTypeAudio preferredTrackID:kCMPersistentTrackID_Invalid];
    if (![a insertTimeRange:CMTimeRangeMake (kCMTimeZero, CMTimeMinimum (source.duration, soundtrack.duration))
                   ofTrack:[soundtrack tracksWithMediaType:AVMediaTypeAudio][0] atTime:kCMTimeZero error:&error]) return false;
    auto exporter = [[AVAssetExportSession alloc] initWithAsset:composition presetName:AVAssetExportPresetPassthrough];
    exporter.outputURL = [NSURL fileURLWithPath:outputPath]; exporter.outputFileType = AVFileTypeQuickTimeMovie;
    dispatch_semaphore_t done = dispatch_semaphore_create (0);
    [exporter exportAsynchronouslyWithCompletionHandler:^{ dispatch_semaphore_signal (done); }];
    dispatch_semaphore_wait (done, DISPATCH_TIME_FOREVER);
    if (exporter.status != AVAssetExportSessionStatusCompleted) { fprintf (stderr, "%s\n", exporter.error.description.UTF8String); return false; }
    return true;
}

static bool remuxMP4 (NSString* inputPath, NSString* outputPath)
{
    auto asset = [AVURLAsset URLAssetWithURL:[NSURL fileURLWithPath:inputPath] options:nil];
    auto exporter = [[AVAssetExportSession alloc] initWithAsset:asset presetName:AVAssetExportPresetPassthrough];
    exporter.outputURL = [NSURL fileURLWithPath:outputPath]; exporter.outputFileType = AVFileTypeMPEG4;
    dispatch_semaphore_t done = dispatch_semaphore_create (0);
    [exporter exportAsynchronouslyWithCompletionHandler:^{ dispatch_semaphore_signal (done); }];
    dispatch_semaphore_wait (done, DISPATCH_TIME_FOREVER);
    if (exporter.status != AVAssetExportSessionStatusCompleted) { fprintf (stderr, "%s\n", exporter.error.description.UTF8String); return false; }
    return true;
}

int main (int argc, char** argv)
{
    @autoreleasepool
    {
        if (argc < 3) { fprintf (stderr, "Usage: probe generate|generate-vfr|decode|present file.mov\n"); return 2; }
        const auto mode = [NSString stringWithUTF8String:argv[1]];
        auto path = [NSString stringWithUTF8String:argv[2]];
        const bool ok = [mode isEqualToString:@"remux-mp4"] && argc == 4
            ? remuxMP4 (path, [NSString stringWithUTF8String:argv[3]])
            : [mode isEqualToString:@"mux"] && argc == 5
            ? mux (path, [NSString stringWithUTF8String:argv[3]], [NSString stringWithUTF8String:argv[4]])
            : [mode hasPrefix:@"generate"] ? makeMovie (path, [mode isEqualToString:@"generate-vfr"],
                argc > 3 ? atoi (argv[3]) : 120, argc > 4 ? atoi (argv[4]) : 30000,
                argc > 5 ? atoi (argv[5]) : 1001, [mode isEqualToString:@"generate-tc"])
                                                    : probe (path, [mode isEqualToString:@"present"]);
        printf ("%s\n", ok ? "PASS" : "FAIL");
        return ok ? 0 : 1;
    }
}
