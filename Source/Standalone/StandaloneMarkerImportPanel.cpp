#include "StandaloneMarkerImportPanel.h"
#include <charconv>
namespace levelscope::standalone
{
namespace mk=proposal::markers;
class MarkerImportWindow::Panel final : public juce::Component,private juce::Timer,private juce::ListBoxModel
{
public:
    Panel(StandaloneSessionModel& s,bool pasted,std::function<bool()> allowed,std::function<void()> changed)
        :session(s),canCommit(std::move(allowed)),onChanged(std::move(changed))
    {
        format.addItem("Nuendo CSV point markers",1);format.addItem("Pro Tools Session Info text",2);format.addItem("Pasted positions",3);format.setSelectedId(pasted?3:1,juce::dontSendNotification);
        units.addItem("Audio sample frames",1);units.addItem("Exact seconds",2);units.addItem("Video frame counts",3);units.addItem("Non-drop HH:MM:SS:FF",4);units.setSelectedId(3,juce::dontSendNotification);
        rate.addItem("24/1",1);rate.addItem("25/1",2);rate.addItem("30/1",3);rate.addItem("24000/1001",4);rate.addItem("30000/1001",5);rate.setSelectedId(2,juce::dontSendNotification);
        delimiter.addItem("Comma",1);delimiter.addItem("Semicolon",2);delimiter.addItem("Tab",3);delimiter.setSelectedId(1,juce::dontSendNotification);
        view.addItem("Original marker rows and exact mapping",1);view.addItem("Complete proposed section map",2);view.addItem("Complete changes and disclosures",3);view.setSelectedId(1,juce::dontSendNotification);
        origin.setText("0");positionColumn.setText("1");nameColumn.setText("2");typeColumn.setText("0");initialName.setText("Programme start");header.setToggleState(true,juce::dontSendNotification);
        paste.setMultiLine(true);paste.setReturnKeyStartsNewLine(true);paste.setInputRestrictions(int(mk::maximumBytes));paste.setTextToShowWhenEmpty("Paste one position per line, or a comma-separated integer frame list",juce::Colours::grey);
        details.setMultiLine(true);details.setReadOnly(true);details.setFont(13.f);status.setJustificationType(juce::Justification::topLeft);
        for(auto* c:std::initializer_list<juce::Component*>{&format,&units,&rate,&delimiter,&view,&origin,&positionColumn,&nameColumn,&typeColumn,&initialName,&paste,&details,&header,&pointsConfirmed,&unitsConfirmed,&choose,&previewButton,&cancelButton,&confirmButton,&fileLabel,&mappingLabel,&columnsLabel,&originLabel,&initialLabel,&status,&list})addAndMakeVisible(c);
        format.setName("Marker export format");units.setName("Marker position units");rate.setName("Exact marker frame rate");delimiter.setName("Marker CSV separator");origin.setName("Programme origin in marker units");paste.setName("Pasted marker positions");positionColumn.setName("Position column one based");nameColumn.setName("Name column zero means absent");typeColumn.setName("Type column zero means absent");initialName.setName("Initial section name");view.setName("Marker preview view");details.setName("Selected marker or section detail");status.setName("Marker import status");
        mappingLabel.setText("Position units / exact rate",juce::dontSendNotification);originLabel.setText("Origin (audio sample zero)",juce::dontSendNotification);columnsLabel.setText("CSV columns: position / name / type (1-based; 0 = absent except position)",juce::dontSendNotification);initialLabel.setText("Initial section name",juce::dontSendNotification);
        format.onChange=[this]{invalidate();resized();};units.onChange=[this]{invalidate();};rate.onChange=[this]{invalidate();};delimiter.onChange=[this]{invalidate();};header.onClick=[this]{invalidate();};pointsConfirmed.onClick=[this]{invalidate();};unitsConfirmed.onClick=[this]{invalidate();};
        for(auto* editor:{&origin,&positionColumn,&nameColumn,&typeColumn,&initialName,&paste})editor->onTextChange=[this]{invalidate();};
        view.onChange=[this]{list.deselectAllRows();list.updateContent();list.repaint();details.clear();};
        choose.onClick=[this]{chooser=std::make_unique<juce::FileChooser>("Choose point-marker export",juce::File(),"*.csv;*.tsv;*.txt");auto safe=juce::Component::SafePointer<Panel>(this);chooser->launchAsync(juce::FileBrowserComponent::openMode|juce::FileBrowserComponent::canSelectFiles,[safe](const juce::FileChooser& c){if(safe&&c.getResult()!=juce::File()){safe->file=c.getResult();safe->fileLabel.setText(safe->file.getFileName(),juce::dontSendNotification);safe->invalidate();}});};
        previewButton.onClick=[this]{begin(false);};confirmButton.onClick=[this]{begin(true);};cancelButton.onClick=[this]{cancel();};
        list.setRowHeight(25);setSize(1120,790);setStatus("Confirm the exported marker units and point subset. Real DAW compatibility requires tested exports; no audio or video is opened.");startTimerHz(20);resized();
    }
    ~Panel()override{stopTimer();if(control)control->store(true);if(worker.valid())worker.wait();}
    void cancel(){invalidate();setStatus("Marker import cancelled; annotations, audio and measurements unchanged.");}
    void paint(juce::Graphics& g)override{g.fillAll(juce::Colour(0xff151b24));}
    void resized()override
    {
        auto b=getLocalBounds().reduced(12),r=b.removeFromTop(30);format.setBounds(r.removeFromLeft(280));choose.setBounds(r.removeFromLeft(170));fileLabel.setBounds(r);b.removeFromTop(5);
        const bool pasted=format.getSelectedId()==3;paste.setVisible(pasted);choose.setEnabled(!pasted);fileLabel.setVisible(!pasted);if(pasted){paste.setBounds(b.removeFromTop(64));b.removeFromTop(5);}
        r=b.removeFromTop(30);mappingLabel.setBounds(r.removeFromLeft(170));units.setBounds(r.removeFromLeft(215));rate.setBounds(r.removeFromLeft(130));originLabel.setBounds(r.removeFromLeft(180));origin.setBounds(r);
        r=b.removeFromTop(30);columnsLabel.setBounds(r.removeFromLeft(540));positionColumn.setBounds(r.removeFromLeft(60));nameColumn.setBounds(r.removeFromLeft(60));typeColumn.setBounds(r.removeFromLeft(60));header.setBounds(r.removeFromLeft(90));delimiter.setBounds(r);
        const bool csv=format.getSelectedId()==1;
        columnsLabel.setText(format.getSelectedId()==2?"Pro Tools: audio frames = TIME REFERENCE / Samples; other units = LOCATION":"CSV columns: position / name / type (1-based; 0 = absent except position)",juce::dontSendNotification);
        for(auto* c:std::initializer_list<juce::Component*>{&positionColumn,&nameColumn,&typeColumn,&header,&delimiter})c->setEnabled(csv);
        r=b.removeFromTop(30);initialLabel.setBounds(r.removeFromLeft(170));initialName.setBounds(r.removeFromLeft(280));pointsConfirmed.setBounds(r);
        unitsConfirmed.setBounds(b.removeFromTop(30));r=b.removeFromTop(32);previewButton.setBounds(r.removeFromLeft(175));cancelButton.setBounds(r.removeFromLeft(175));confirmButton.setBounds(r.removeFromLeft(250));view.setBounds(r);
        status.setBounds(b.removeFromTop(62));details.setBounds(b.removeFromBottom(95));b.removeFromBottom(5);list.setBounds(b);
    }
private:
    StandaloneSessionModel& session;std::function<bool()> canCommit;std::function<void()> onChanged;
    juce::ComboBox format,units,rate,delimiter,view;
    juce::TextEditor origin,positionColumn,nameColumn,typeColumn,initialName,paste,details;
    juce::ToggleButton header{"Header"},pointsConfirmed{"I confirm this export contains point markers only"},unitsConfirmed{"I confirm the actual position-field units, rate and origin; no inference from filenames or general DAW settings"};
    juce::TextButton choose{"Choose marker file..."},previewButton{"Preview mapping"},cancelButton{"Cancel import"},confirmButton{"Confirm: replace annotations"};
    juce::Label fileLabel,mappingLabel,columnsLabel,originLabel,initialLabel,status;
    juce::ListBox list{"Marker import preview",this};juce::File file;std::unique_ptr<juce::FileChooser> chooser;
    std::shared_ptr<std::atomic<bool>> control;
    struct Result{std::shared_ptr<const mk::Preview> preview;juce::String hash,error;uint64_t generation=0;bool confirm=false;};
    std::future<Result> worker;std::shared_ptr<const mk::Preview> accepted;uint64_t generation=0;
    void setStatus(const juce::String& text){status.setText(text,juce::dontSendNotification);status.setTooltip(text);}
    void invalidate(){++generation;if(control)control->store(true);accepted.reset();list.updateContent();list.repaint();details.clear();confirmButton.setEnabled(false);}
    int column(const juce::TextEditor& e){const auto s=e.getText().trim().toStdString();int n=0;auto r=std::from_chars(s.data(),s.data()+s.size(),n);if(r.ec!=std::errc()||r.ptr!=s.data()+s.size()||n<0||n>128)throw std::invalid_argument("Columns must be integers from 0 to 128");return n-1;}
    mk::Options options()
    {
        mk::Options o;o.format=mk::Format(format.getSelectedId()-1);o.unit=mk::Unit(units.getSelectedId()-1);o.separator=delimiter.getSelectedId()==1?',':delimiter.getSelectedId()==2?';':'\t';o.header=header.getToggleState();o.positionColumn=column(positionColumn);o.nameColumn=column(nameColumn);o.typeColumn=column(typeColumn);
        const int rates[]={24,25,30,24000,30000};o.rateNumerator=rates[rate.getSelectedId()-1];o.rateDenominator=rate.getSelectedId()>3?1001:1;o.origin=origin.getText();o.initialName=initialName.getText();o.confirmedPointMarkers=pointsConfirmed.getToggleState();o.confirmedUnits=unitsConfirmed.getToggleState();return o;
    }
    static std::string read(const juce::File& file,const std::string& pasted,bool usePaste,const std::shared_ptr<std::atomic<bool>>& cancel)
    {
        if(usePaste)return pasted;auto stream=file.createInputStream();if(!stream||stream->getTotalLength()<=0||stream->getTotalLength()>int64_t(mk::maximumBytes))throw std::invalid_argument("Marker file is missing, empty or exceeds 8 MiB");
        std::string bytes;std::array<char,16384> block{};for(;;){if(cancel->load())throw std::invalid_argument("Marker import cancelled");int n=stream->read(block.data(),int(block.size()));if(n<=0)break;bytes.append(block.data(),size_t(n));if(bytes.size()>mk::maximumBytes)throw std::invalid_argument("Marker file exceeds 8 MiB");}if(stream->getStatus().failed())throw std::invalid_argument("Marker file read failed");return bytes;
    }
    void begin(bool confirm)
    {
        if(worker.valid()){setStatus("Wait for the cancelled operation to finish, then preview again.");return;}
        if(session.audioEditsBlocked()||session.job.isActive()||(canCommit&&!canCommit())){setStatus("Stop playback and finish the active file, measurement, proposal or export operation first.");return;}
        if(confirm&&!accepted)return;
        try
        {
            auto o=options();const auto text=o.format==mk::Format::pastedPositions?paste.getText().toStdString():std::string();
            auto snapshot=session;auto old=accepted;const auto token=++generation;control=std::make_shared<std::atomic<bool>>(false);auto flag=control;const auto inputFile=file;
            setStatus(confirm?"Checking marker bytes and current session before atomic confirmation...":"Reading and mapping point markers...");
            worker=std::async(std::launch::async,[o,text,snapshot,old,token,flag,inputFile,confirm]()mutable
            {
                Result result;result.generation=token;result.confirm=confirm;
                try
                {
                    auto bytes=read(inputFile,text,o.format==mk::Format::pastedPositions,flag);result.hash=mk::contentHash(bytes);
                    if(confirm)result.preview=old;
                    else result.preview=std::make_shared<const mk::Preview>(mk::prepare(bytes,o.format==mk::Format::pastedPositions?"Pasted positions":inputFile.getFileName(),o,snapshot.sectionDocument.map(),snapshot.sessionId,[flag]{return flag->load();}));
                    if(!preflightMarkerImport(snapshot,*result.preview,result.error))result.preview.reset();
                }catch(const std::exception& e){result.error=e.what();}return result;
            });
        }catch(const std::exception& e){setStatus(e.what());}
    }
    void timerCallback()override
    {
        if(worker.valid()&&worker.wait_for(std::chrono::seconds(0))==std::future_status::ready)
        {
            auto result=worker.get();if(result.generation==generation&&control&&!control->load())
            {
                if(!result.preview){accepted.reset();setStatus(result.error);}
                else if(result.confirm)
                {
                    juce::String error;if((!canCommit||canCommit())&&commitMarkerImport(session,*result.preview,result.hash,error)){accepted.reset();setStatus(error);if(onChanged)onChanged();}
                    else setStatus(error.isNotEmpty()?error:"Active operation prevents confirmation; preview retained.");
                }
                else if(session.sessionId!=result.preview->runtime||!session.sectionDocument.map().valid()||session.sectionDocument.map().digest()!=result.preview->previousMap)
                {accepted.reset();setStatus("Stale preview: session or annotations changed; preview again.");}
                else{accepted=result.preview;setStatus("Preview only: "+juce::String(int(accepted->receipt.rows.size()))+" input rows; "+juce::String(int(accepted->receipt.imported.sections.size()))+" proposed sections. Review all three views; confirm replaces annotations only. Quantization is the exact amount moved earlier, in samples.");}
                list.deselectAllRows();list.updateContent();list.repaint();details.clear();
            }
            control.reset();
        }
        const bool allowed=!worker.valid()&&!session.audioEditsBlocked()&&!session.job.isActive()&&(!canCommit||canCommit());
        previewButton.setEnabled(allowed);confirmButton.setEnabled(allowed&&accepted!=nullptr);
        if(accepted&&(session.sessionId!=accepted->runtime||!session.sectionDocument.map().valid()||session.sectionDocument.map().digest()!=accepted->previousMap))
        {invalidate();setStatus("Preview became stale; annotations unchanged. Preview again.");}
    }
    int getNumRows()override
    {if(!accepted)return 0;if(view.getSelectedId()==1)return int(accepted->receipt.rows.size());if(view.getSelectedId()==2)return int(accepted->receipt.imported.sections.size());return accepted->diff.size()+accepted->receipt.notices.size();}
    juce::String rowText(int n)const
    {
        if(!accepted||n<0)return {};
        if(view.getSelectedId()==1)
        {if(size_t(n)>=accepted->receipt.rows.size())return {};const auto& r=accepted->receipt.rows[size_t(n)];return "Row "+juce::String(r.line)+" | "+r.position+" | "+r.name+" -> sample "+juce::String(r.frame)+" | seconds "+r.seconds+" | earlier by "+r.quantization+" samples"+(r.coalesced?" | COALESCED":"")+(r.frame==0?" | START":r.frame==accepted->receipt.imported.frames?" | EOF (no following section)":"");}
        if(view.getSelectedId()==2){if(size_t(n)>=accepted->receipt.imported.sections.size())return {};const auto& s=accepted->receipt.imported.sections[size_t(n)];return "ID "+juce::String(s.id)+" ["+juce::String(s.begin)+", "+juce::String(s.end)+") | "+s.name;}
        return n<accepted->diff.size()?accepted->diff[n]:accepted->receipt.notices[n-accepted->diff.size()];
    }
    void paintListBoxItem(int n,juce::Graphics& g,int w,int h,bool selected)override{if(selected)g.fillAll(juce::Colour(0xff28405a));g.setColour(juce::Colours::white);g.setFont(13.f);g.drawText(rowText(n),5,0,w-10,h,juce::Justification::centredLeft,true);}
    void selectedRowsChanged(int row)override{details.setText(rowText(row),false);}
};
MarkerImportWindow::MarkerImportWindow(StandaloneSessionModel& s,bool paste,std::function<bool()> canCommit,std::function<void()> changed)
    :juce::DocumentWindow("Import point markers - annotations only",juce::Colour(0xff151b24),juce::DocumentWindow::closeButton)
{setUsingNativeTitleBar(true);panel=new Panel(s,paste,std::move(canCommit),std::move(changed));setContentOwned(panel,true);setResizable(true,false);centreWithSize(1120,790);setVisible(true);}
void MarkerImportWindow::closeButtonPressed(){panel->cancel();setVisible(false);}
}
