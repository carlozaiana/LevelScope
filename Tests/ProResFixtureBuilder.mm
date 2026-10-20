// Offline test media only. No audio device or player is constructed.
#import <AVFoundation/AVFoundation.h>
#include <cstdio>
#include <thread>
#include <chrono>
#include <cstring>

static bool finish(AVAssetWriter* writer, AVAssetWriterInput* input, CMTime end) {
    [input markAsFinished]; [writer endSessionAtSourceTime:end];
    dispatch_semaphore_t done=dispatch_semaphore_create(0);
    [writer finishWritingWithCompletionHandler:^{dispatch_semaphore_signal(done);}];
    dispatch_semaphore_wait(done,DISPATCH_TIME_FOREVER);
    if(writer.status!=AVAssetWriterStatusCompleted) fprintf(stderr,"%s\n",writer.error.description.UTF8String);
    return writer.status==AVAssetWriterStatusCompleted;
}
// 256 unique, detailed moving frames with 8-bit frame identity, colour patches
// and a ramp. Long fixtures repeat compressed frames but preserve unique PTS.
static bool generate(NSString* path, NSString* codec, int width,int height,int count,int scale,int step,NSString* variant) {
    NSError* error=nil;
    auto writer=[[AVAssetWriter alloc]initWithURL:[NSURL fileURLWithPath:path] fileType:AVFileTypeQuickTimeMovie error:&error];
    writer.movieTimeScale=scale;
    auto colour=[@{AVVideoColorPrimariesKey:AVVideoColorPrimaries_ITU_R_709_2,
                   AVVideoTransferFunctionKey:AVVideoTransferFunction_ITU_R_709_2,
                   AVVideoYCbCrMatrixKey:AVVideoYCbCrMatrix_ITU_R_709_2} mutableCopy];
    if([variant isEqual:@"conflict"]) colour[AVVideoYCbCrMatrixKey]=AVVideoYCbCrMatrix_ITU_R_601_4;
    if([variant isEqual:@"hdr"]) colour[AVVideoTransferFunctionKey]=@"SMPTE_ST_2084_PQ";
    auto settings=[@{AVVideoCodecKey:codec,AVVideoWidthKey:@(width),AVVideoHeightKey:@(height),
       AVVideoColorPropertiesKey:colour,
       AVVideoPixelAspectRatioKey:@{AVVideoPixelAspectRatioHorizontalSpacingKey:@1,AVVideoPixelAspectRatioVerticalSpacingKey:@1}} mutableCopy];
    auto input=[AVAssetWriterInput assetWriterInputWithMediaType:AVMediaTypeVideo outputSettings:settings];
    input.mediaTimeScale=scale;
    auto adaptor=[AVAssetWriterInputPixelBufferAdaptor assetWriterInputPixelBufferAdaptorWithAssetWriterInput:input
        sourcePixelBufferAttributes:@{(id)kCVPixelBufferPixelFormatTypeKey:@(kCVPixelFormatType_32BGRA),
        (id)kCVPixelBufferWidthKey:@(width),(id)kCVPixelBufferHeightKey:@(height)}];
    if(![writer canAddInput:input]) return false;
    [writer addInput:input]; if(![writer startWriting]) return false;
    [writer startSessionAtSourceTime:kCMTimeZero];
    int64_t ticks=0;
    for(int f=0;f<count;++f) { @autoreleasepool {
        while(!input.readyForMoreMediaData && writer.status==AVAssetWriterStatusWriting) std::this_thread::sleep_for(std::chrono::milliseconds(1));
        CVPixelBufferRef pixel=nullptr;
        if(CVPixelBufferPoolCreatePixelBuffer(nullptr,adaptor.pixelBufferPool,&pixel)) return false;
        CVPixelBufferLockBaseAddress(pixel,0);
        auto* base=(unsigned char*)CVPixelBufferGetBaseAddress(pixel);auto stride=CVPixelBufferGetBytesPerRow(pixel);
        for(int y=0;y<height;++y) for(int x=0;x<width;++x) {
            auto* p=base+y*stride+x*4;
            unsigned n=(unsigned)(x+f*7)*1664525u ^ (unsigned)(y+f*3)*1013904223u;
            p[0]=32+((n>>9)%192);p[1]=32+((n>>17)%192);p[2]=32+((n>>23)%192);p[3]=255;
            // Top quarter: neutral ramp; lower quarter: binary frame bits.
            if(y<height/4) p[0]=p[1]=p[2]=(unsigned char)(16+219*x/width);
            if(y>=height*3/4) p[0]=p[1]=p[2]=((f>>(x*8/width))&1)?235:16;
            // Middle band: RGB/white/black/grey patches.
            if(y>=height/3 && y<height/2) {
                int patch=x*6/width; p[0]=p[1]=p[2]=16;
                if(patch<3)p[2-patch]=235;
                if(patch==3)p[0]=p[1]=p[2]=235;
                if(patch==5)p[0]=p[1]=p[2]=128;
            }
        }
        CVBufferSetAttachment(pixel,kCVImageBufferColorPrimariesKey,kCVImageBufferColorPrimaries_ITU_R_709_2,kCVAttachmentMode_ShouldPropagate);
        CVBufferSetAttachment(pixel,kCVImageBufferTransferFunctionKey,kCVImageBufferTransferFunction_ITU_R_709_2,kCVAttachmentMode_ShouldPropagate);
        CVBufferSetAttachment(pixel,kCVImageBufferYCbCrMatrixKey,kCVImageBufferYCbCrMatrix_ITU_R_709_2,kCVAttachmentMode_ShouldPropagate);
        CVPixelBufferUnlockBaseAddress(pixel,0);
        bool ok=[adaptor appendPixelBuffer:pixel withPresentationTime:CMTimeMake(ticks,scale)];CVPixelBufferRelease(pixel);
        if(!ok){fprintf(stderr,"%s\n",writer.error.description.UTF8String);return false;}
        ticks+=([variant isEqual:@"vfr"] && f%3==2)?2*step:step;
    }}
    return finish(writer,input,CMTimeMake(ticks,scale));
}
static bool repeat(NSString* seed,NSString* path,int count,int scale,int step) {
    auto asset=[AVURLAsset URLAssetWithURL:[NSURL fileURLWithPath:seed] options:nil];
    auto track=[asset tracksWithMediaType:AVMediaTypeVideo].firstObject; if(!track)return false;
    NSError* error=nil;
    auto writer=[[AVAssetWriter alloc]initWithURL:[NSURL fileURLWithPath:path] fileType:AVFileTypeQuickTimeMovie error:&error];
    writer.movieTimeScale=scale;
    auto input=[AVAssetWriterInput assetWriterInputWithMediaType:AVMediaTypeVideo outputSettings:nil sourceFormatHint:(__bridge CMFormatDescriptionRef)track.formatDescriptions.firstObject];
    input.mediaTimeScale=scale; [writer addInput:input];
    if(![writer startWriting])return false;[writer startSessionAtSourceTime:kCMTimeZero];
    int index=0;
    while(index<count) { @autoreleasepool {
        auto reader=[[AVAssetReader alloc]initWithAsset:asset error:&error];
        auto output=[[AVAssetReaderTrackOutput alloc]initWithTrack:track outputSettings:nil];[reader addOutput:output];
        if(![reader startReading])return false;
        while(index<count) { @autoreleasepool {
            auto sample=[output copyNextSampleBuffer];if(!sample)break;
            if(CMSampleBufferGetNumSamples(sample)==0){CFRelease(sample);continue;}
            if(CMSampleBufferGetNumSamples(sample)!=1){CFRelease(sample);return false;}
            CMSampleTimingInfo timing={CMTimeMake(step,scale),CMTimeMake((int64_t)index*step,scale),kCMTimeInvalid};
            CMSampleBufferRef retimed=nullptr;
            auto code=CMSampleBufferCreateCopyWithNewTiming(nullptr,sample,1,&timing,&retimed);CFRelease(sample);
            if(code)return false;
            while(!input.readyForMoreMediaData && writer.status==AVAssetWriterStatusWriting) std::this_thread::sleep_for(std::chrono::milliseconds(1));
            bool ok=[input appendSampleBuffer:retimed];CFRelease(retimed);if(!ok)return false;++index;
        }}
        [reader cancelReading];
    }}
    return finish(writer,input,CMTimeMake((int64_t)count*step,scale));
}
int main(int argc,char** argv) {
 @autoreleasepool {
    bool ok=false;
    if(argc==10 && !strcmp(argv[1],"generate"))
        ok=generate(@(argv[2]),@(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),@(argv[9]));
    else if(argc==7 && !strcmp(argv[1],"repeat"))
        ok=repeat(@(argv[2]),@(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]));
    else {fprintf(stderr,"generate output.mov apco|apcs|apcn width height frames timescale step normal|vfr|conflict|hdr\nrepeat seed.mov output.mov frames timescale step\n");return 2;}
    puts(ok?"PASS":"FAIL");return ok?0:1;
 }
}
