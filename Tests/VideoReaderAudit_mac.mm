// Test-only observation of every AVAssetReaderTrackOutput construction, including
// background threads. The application never links this Objective-C interception.
#import <AVFoundation/AVFoundation.h>
#import <objc/runtime.h>
#include <atomic>
static std::atomic<int> counts[2];
static IMP original = nullptr;
static id audited(id self,SEL selector,AVAssetTrack* track,NSDictionary* settings)
{
    if([track.mediaType isEqualToString:AVMediaTypeVideo]) ++counts[0];
    if([track.mediaType isEqualToString:AVMediaTypeAudio]) ++counts[1];
    return ((id(*)(id,SEL,AVAssetTrack*,NSDictionary*))original)(self,selector,track,settings);
}
void beginVideoReaderAudit()
{
    counts[0]=0;counts[1]=0;
    if(!original) original=method_setImplementation(class_getInstanceMethod([AVAssetReaderTrackOutput class],@selector(initWithTrack:outputSettings:)),(IMP)audited);
}
int videoReaderAuditCount(int kind){return counts[kind].load();}
