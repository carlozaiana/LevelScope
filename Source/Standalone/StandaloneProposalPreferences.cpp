#include "StandaloneProposalPreferences.h"
namespace levelscope::standalone
{
ProposalInspectionPanel::ProposalInspectionPanel(juce::String text)
{
    // At most 10,000 correction rows; reports enforce their separate 16 MiB bound.
    const auto bytes=text.toStdString();
    for(size_t begin=0;begin<bytes.size();)
    {
        size_t end=std::min(begin+12000,bytes.size());
        if(end<bytes.size())
        {
            const auto line=bytes.rfind('\n',end-1);
            if(line!=std::string::npos&&line>=begin)end=line+1;
            else while(end>begin&&(static_cast<unsigned char>(bytes[end])&0xc0)==0x80)--end;
        }
        pages.push_back(juce::String::fromUTF8(bytes.data()+begin,int(end-begin)));begin=end;
    }
    if(pages.empty())pages.push_back("No proposal available.");
    view.setMultiLine(true);view.setReadOnly(true);view.setScrollbarsShown(true);view.setName("Complete proposal - paginated");
    for(auto* c:std::initializer_list<juce::Component*>{&view,&location,&previous,&next})addAndMakeVisible(c);
    previous.onClick=[this]{if(page>0)--page;refresh();};next.onClick=[this]{nextPage();};setSize(850,650);refresh();
}
juce::String ProposalInspectionPanel::allText()const {juce::String s;for(const auto& p:pages)s+=p;return s;}
void ProposalInspectionPanel::nextPage(){if(page+1<int(pages.size()))++page;refresh();}
void ProposalInspectionPanel::refresh()
{view.setText(pages[size_t(page)],false);view.moveCaretToTop(false);location.setText("Complete inspection - page "+juce::String(page+1)+" / "+juce::String(pages.size())+". No corrections omitted.",juce::dontSendNotification);previous.setEnabled(page>0);next.setEnabled(page+1<int(pages.size()));}
void ProposalInspectionPanel::resized()
{auto r=getLocalBounds().reduced(12);auto footer=r.removeFromBottom(32);previous.setBounds(footer.removeFromLeft(140));next.setBounds(footer.removeFromRight(140));location.setBounds(footer);view.setBounds(r.reduced(0,5));}
std::optional<int> ProposalPreferencesPanel::exactUnits(juce::String text,bool tenths)
{
    text=text.trim();if(text.isEmpty()||text.length()>8)return {};
    const auto dot=text.indexOfChar('.');
    if(!tenths&&dot>=0)return {};
    if(dot>=0)
    {
        if(dot==0||dot!=text.length()-2||!text.substring(0,dot).containsOnly("0123456789")||!text.substring(dot+1).containsOnly("0123456789"))return {};
        return text.substring(0,dot).getIntValue()*10+text.substring(dot+1).getIntValue();
    }
    if(!text.containsOnly("0123456789"))return {};return text.getIntValue()*(tenths?10:1);
}
ProposalPreferencesPanel::ProposalPreferencesPanel(StandaloneSessionModel& s,std::function<void()> c):session(s),changed(std::move(c)),runtime(s.sessionId),baseDigest(s.proposalPolicy.digest())
{
    const char* names[]={"Optional maximum changed processing regions (1-10,000)","Newly affected duration % (0-100; preservation preference)","Additional attenuation / request dB (0-6)","Resulting cumulative local attenuation dB (0-12)","Requested reserve dB (0-1, inside additional/cumulative caps)"};
    addAndMakeVisible(unrestricted);unrestricted.onClick=[this]{fields[0].setEnabled(!unrestricted.getToggleState());};
    for(size_t n=0;n<fields.size();++n){labels[n].setText(names[n],juce::dontSendNotification);fields[n].setName(names[n]);addAndMakeVisible(labels[n]);addAndMakeVisible(fields[n]);}
    status.setMultiLine(true);status.setReadOnly(true);status.setScrollbarsShown(true);
    for(auto* component:std::initializer_list<juce::Component*>{&use,&cancel,&defaults,&status})addAndMakeVisible(component);
    load(s.proposalPolicy);
    cancel.onClick=[this]{load(session.proposalPolicy);baseDigest=session.proposalPolicy.digest();status.setText("Edits cancelled. Active preferences unchanged.",false);};
    defaults.onClick=[this]{load({});status.setText("Draft defaults only. Use these limits to confirm; audio remains unchanged.",false);};
    use.onClick=[this]
    {
        if(runtime!=session.sessionId||baseDigest!=session.proposalPolicy.digest()){status.setText("Session/preferences changed. Close and reopen this editor.",false);return;}
        proposal::ProposalPolicy p;
        for(size_t n=0;n<fields.size();++n)
        {
            if(n==0&&unrestricted.getToggleState())continue;
            const auto value=exactUnits(fields[n].getText(),n!=0);
            if(!value){status.setText("Use exact non-negative decimal values with at most one fractional digit (integer for count). No rounding or clamping.",false);return;}
            switch(n){case 0:p.regionLimit=*value;break;case 1:p.durationThousandths=*value;break;case 2:p.additionalTenths=*value;break;case 3:p.cumulativeTenths=*value;break;case 4:p.reserveTenths=*value;break;}
        }
        juce::String error;if(!session.setProposalPolicy(p,error)){status.setText(error,false);return;}
        baseDigest=p.digest();status.setText("Active preferences confirmed. Measurements/audio unchanged; old proposals are stale.\n"+p.text(),false);if(changed)changed();
    };
    setSize(800,570);
}
void ProposalPreferencesPanel::load(const proposal::ProposalPolicy& p)
{
    unrestricted.setToggleState(!p.regionLimit,juce::dontSendNotification);fields[0].setText(p.regionLimit?juce::String(*p.regionLimit):juce::String(),false);fields[0].setEnabled(bool(p.regionLimit));
    const int values[]={p.durationThousandths,p.additionalTenths,p.cumulativeTenths,p.reserveTenths};
    for(size_t n=1;n<5;++n)fields[n].setText(juce::String(values[n-1]/10.,1),false);
    status.setText("Draft / active snapshot: explicit Use these limits required. No audio changes or automatic proposal.\n"+p.text(),false);
}
void ProposalPreferencesPanel::resized()
{
    auto r=getLocalBounds().reduced(16);unrestricted.setBounds(r.removeFromTop(32));
    for(size_t n=0;n<5;++n){auto row=r.removeFromTop(40);fields[n].setBounds(row.removeFromRight(110).reduced(0,4));labels[n].setBounds(row);}
    auto buttons=r.removeFromTop(36);use.setBounds(buttons.removeFromLeft(200));cancel.setBounds(buttons.removeFromLeft(180));defaults.setBounds(buttons.removeFromLeft(180));status.setBounds(r.reduced(0,10));
}
}
