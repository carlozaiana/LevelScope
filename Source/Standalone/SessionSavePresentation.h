#pragma once
#include "SessionDocument.h"
namespace levelscope::standalone
{
// UI presentation derives compatibility from the prospective serialized snapshot,
// including history-only method requirements. It does not choose a writer version.
struct SessionSavePresentation
{
    int writer=0;
    bool upgradesExisting=false;
    juce::String message,title;
    juce::File suggestedCopy;
    static SessionSavePresentation forDocument(const sessionfile::Document& document,int openedSchema,const juce::File& existing)
    {
        SessionSavePresentation p;p.writer=sessionfile::writerSchema(document);
        p.upgradesExisting=existing!=juce::File()&&openedSchema<p.writer;
        p.message="Saving writes session format "+juce::String(p.writer)+". Readers supporting format "+juce::String(p.writer-1)+" or earlier cannot open it. "
            +juce::String(p.writer==sessionfile::innerSchemaVersion?"Competing-boundary settings or their applied explanation require this version. ":"")
            +"Choose Save As to keep the original older-readable file. Applied audio settings and historical explanations are preserved.";
        p.title="Save LevelScope session (requires format "+juce::String(p.writer)+" reader)";
        p.suggestedCopy=existing==juce::File()?juce::File::getSpecialLocation(juce::File::userDocumentsDirectory).getChildFile("Untitled.levelscope"):
            existing.getSiblingFile(existing.getFileNameWithoutExtension()+(p.upgradesExisting?"-format-"+juce::String(p.writer):"-copy")+".levelscope").getNonexistentSibling();
        return p;
    }
};
}
