#include "MarkerImport.h"
#include <charconv>
#include <map>
#include <set>
#include <numeric>
#include <limits>
#include <stdexcept>

namespace levelscope::proposal::markers
{
namespace
{
void need(bool b,const juce::String& message){if(!b)throw std::invalid_argument(message.toStdString());}
void poll(const Cancel& c){need(!c||!c(),"Marker import cancelled; annotations unchanged");}
int64_t integer(const juce::String& text)
{
    const auto s=text.trim().toStdString();int64_t n=0;const auto r=std::from_chars(s.data(),s.data()+s.size(),n);
    need(!s.empty()&&r.ec==std::errc()&&r.ptr==s.data()+s.size()&&n>=0,"Expected a nonnegative exact integer");return n;
}
int64_t bounded(__int128 n){need(n>=0&&n<=INT64_MAX,"Position arithmetic overflow");return int64_t(n);}
struct Rational
{
    int64_t n=0,d=1;
    Rational(int64_t a=0,int64_t b=1):n(a),d(b){need(n>=0&&d>0,"Invalid rational position");auto g=std::gcd(n,d);n/=g;d/=g;}
    juce::String text()const{return juce::String(n)+"/"+juce::String(d);}
};
Rational subtract(Rational a,Rational b)
{
    auto g=std::gcd(a.d,b.d);const auto n=(__int128)a.n*(b.d/g)-(__int128)b.n*(a.d/g);
    return {bounded(n),bounded((__int128)(a.d/g)*b.d)};
}
Rational scale(Rational a,int64_t n,int64_t d)
{
    auto g=std::gcd(a.n,d);a.n/=g;d/=g;g=std::gcd(n,a.d);n/=g;a.d/=g;
    return {bounded((__int128)a.n*n),bounded((__int128)a.d*d)};
}
int nominal(const Options& o)
{
    if(o.rateDenominator==1&&(o.rateNumerator==24||o.rateNumerator==25||o.rateNumerator==30))return o.rateNumerator;
    if(o.rateDenominator==1001&&(o.rateNumerator==24000||o.rateNumerator==30000))return o.rateNumerator/1000;
    throw std::invalid_argument("Select supported exact rate: 24, 25, 30, 24000/1001 or 30000/1001; drop-frame unsupported");
}
Rational parsePosition(const juce::String& raw,const Options& o,int audioRate)
{
    const auto s=raw.trim();need(s.isNotEmpty()&&raw.length()<=128,"Position missing or longer than 128 characters");
    if(o.unit==Unit::audioFrames)return {integer(s),audioRate};
    if(o.unit==Unit::videoFrames)return scale({integer(s)},o.rateDenominator,o.rateNumerator);
    if(o.unit==Unit::nonDropTimecode)
    {
        need(!s.containsChar(';'),"Drop-frame timecode is unsupported");
        auto parts=juce::StringArray::fromTokens(s,":","");need(parts.size()==4,"Use non-drop HH:MM:SS:FF");
        for(auto& p:parts)need(p.length()==2&&p.containsOnly("0123456789"),"Use two digits in each non-drop timecode field");
        auto h=integer(parts[0]),m=integer(parts[1]),sec=integer(parts[2]),f=integer(parts[3]);const auto base=nominal(o);
        need(h<24&&m<60&&sec<60&&f<base,"Invalid timecode field or unsupported 24-hour wrap");
        return scale({((h*60+m)*60+sec)*base+f},o.rateDenominator,o.rateNumerator);
    }
    need(o.unit==Unit::seconds,"Unsupported marker unit");
    if(s.containsChar('/'))
    {auto p=juce::StringArray::fromTokens(s,"/","");need(p.size()==2,"Use one numerator/denominator");return {integer(p[0]),integer(p[1])};}
    if(s.containsChar('.'))
    {
        auto p=juce::StringArray::fromTokens(s,".","");need(p.size()==2&&p[0].isNotEmpty()&&p[1].isNotEmpty()&&p[1].length()<=9,"Use exact decimal seconds with 1-9 fractional digits, or rational seconds");
        need(p[0].containsOnly("0123456789")&&p[1].containsOnly("0123456789"),"Invalid decimal seconds");int64_t d=1;for(int i=0;i<p[1].length();++i)d*=10;
        return {bounded((__int128)integer(p[0])*d+integer(p[1])),d};
    }
    return {integer(s)};
}
void mapped(Row& row,const Options& o,const SectionMap& map)
{
    try
    {
        auto t=subtract(parsePosition(row.position,o,map.sampleRate),parsePosition(o.origin,o,map.sampleRate));
        auto samples=scale(t,map.sampleRate,1);need((__int128)samples.n<=(__int128)map.frames*samples.d,"Position outside programme duration (before quantization)");
        row.frame=samples.n/samples.d;row.seconds=t.text();
        row.quantization=Rational(samples.n%samples.d,samples.d).text(); // positive amount moved earlier, in samples
    }
    catch(const std::exception& e){throw std::invalid_argument(("Row "+juce::String(row.line)+": "+e.what()).toStdString());}
}
using Cells=std::vector<juce::String>;
struct Record{int line;Cells cells;};
std::vector<Record> csv(const std::string& bytes,char delimiter,const Cancel& cancel)
{
    need(delimiter==','||delimiter==';'||delimiter=='\t',"Choose comma, semicolon or tab separator");
    std::vector<Record> rows;Record row{1,{}};std::string cell;bool quoted=false,closed=false;int line=1;
    auto field=[&]{need(cell.size()<=1024*1024,"CSV field exceeds 1 MiB");row.cells.push_back(juce::String::fromUTF8(cell.data(),int(cell.size())));cell.clear();closed=false;need(row.cells.size()<=128,"More than 128 CSV columns");};
    for(size_t i=0;i<bytes.size();++i)
    {
        if((i&4095)==0)poll(cancel);char c=bytes[i];
        if(quoted){if(c=='"'){if(i+1<bytes.size()&&bytes[i+1]=='"'){cell+='"';++i;}else{quoted=false;closed=true;}}else{cell+=c;if(c=='\n')++line;}continue;}
        if(c=='"'){need(cell.empty()&&!closed,"Quote inside unquoted CSV field");quoted=true;continue;}
        if(c==delimiter){field();continue;}
        if(c=='\r'||c=='\n')
        {
            if(c=='\r'&&i+1<bytes.size()&&bytes[i+1]=='\n')++i;
            field();rows.push_back(std::move(row));row={++line,{}};need(rows.size()<=maximumRows+1,"More than 10,000 marker rows");continue;
        }
        need(!closed,"Unexpected characters after closing CSV quote");cell+=c;
    }
    need(!quoted,"Unclosed CSV quote");if(!cell.empty()||closed||!row.cells.empty()){field();rows.push_back(std::move(row));}
    return rows;
}
juce::String canonical(juce::String s){return s.toLowerCase().removeCharacters(" _-\t");}
bool point(const juce::String& s){return s.equalsIgnoreCase("point")||s.equalsIgnoreCase("marker")||s.equalsIgnoreCase("position");}
void label(const juce::String& s){need((s.isEmpty()||s.trim().isNotEmpty())&&s.length()<=128&&!s.containsChar('\n')&&!s.containsChar('\r')&&!s.containsChar('\t'),"Marker name must be single-line text of at most 128 characters");}
std::vector<Row> readNuendo(const std::string& bytes,const Options& o,juce::StringArray& notices,const Cancel& cancel)
{
    auto records=csv(bytes,o.separator,cancel);need(!records.empty(),"No marker rows");const auto columns=records.front().cells.size();
    need(o.positionColumn>=0&&size_t(o.positionColumn)<columns&&o.nameColumn>=-1&&(o.nameColumn<0||size_t(o.nameColumn)<columns)
         &&o.typeColumn>=-1&&(o.typeColumn<0||size_t(o.typeColumn)<columns),"Selected CSV column does not exist");
    need(o.positionColumn!=o.nameColumn&&o.positionColumn!=o.typeColumn&&(o.typeColumn<0||o.typeColumn!=o.nameColumn),"Select distinct position, name and type columns");
    std::vector<int> endColumns,durationColumns,trackColumns;std::map<int,juce::String> tracks;
    for(size_t c=0;c<columns;++c)
    {
        auto name=o.header?records.front().cells[c]:"Column "+juce::String(int(c+1));const auto key=canonical(name);
        if(o.header&&(key=="type"||key=="markertype"))need(int(c)==o.typeColumn,"Select the exported marker-type column explicitly");
        if(o.header&&(key=="track"||key=="markertrack"||key=="trackname"))trackColumns.push_back(int(c));
        if(o.header&&(key=="end"||key=="endposition"||key=="timecodeout"))endColumns.push_back(int(c));
        if(o.header&&(key=="duration"||key=="length"))durationColumns.push_back(int(c));
        if(int(c)!=o.positionColumn&&int(c)!=o.nameColumn&&int(c)!=o.typeColumn)notices.add("Unused column: "+name);
    }
    std::vector<Row> result;
    for(size_t n=o.header?1:0;n<records.size();++n)
    {
        poll(cancel);const auto& r=records[n];need(r.cells.size()==columns,"Row "+juce::String(r.line)+": column count differs from header/first row");
        for(auto c:trackColumns){auto track=r.cells[size_t(c)].trim();auto known=tracks.find(c);need(known==tracks.end()||known->second==track,"Multiple marker tracks are unsupported; export one point-marker track");tracks.emplace(c,track);}
        if(o.typeColumn>=0)need(point(r.cells[size_t(o.typeColumn)].trim()),"Row "+juce::String(r.line)+": unsupported/ambiguous marker type; point markers only");
        for(auto c:endColumns)need(r.cells[size_t(c)].trim().isEmpty()||r.cells[size_t(c)].trim()==r.cells[size_t(o.positionColumn)].trim(),"Range/cycle end differs from position; point-only export required");
        for(auto c:durationColumns){auto d=r.cells[size_t(c)].trim();need(d.isEmpty()||d=="0"||d=="0.0"||d=="00:00:00:00","Nonzero or ambiguous marker duration; point-only export required");}
        Row row;row.line=r.line;row.position=r.cells[size_t(o.positionColumn)];if(o.nameColumn>=0)row.name=r.cells[size_t(o.nameColumn)];label(row.name);result.push_back(std::move(row));
    }
    return result;
}
std::vector<Row> readProTools(const std::string& bytes,const Options& o,int audioRate,juce::StringArray& notices,const Cancel& cancel)
{
    // Export Session Info as Text is a sectioned, tab-delimited document. Never
    // interpret track event tables as markers or infer units from a session heading.
    auto lines=juce::StringArray::fromLines(juce::String::fromUTF8(bytes.data(),int(bytes.size())));
    int section=-1;for(int n=0;n<lines.size();++n)if(canonical(lines[n].trim())=="markers"||canonical(lines[n].trim())=="markerslisting")
    {need(section<0,"Multiple MARKERS sections are ambiguous");section=n;}
    need(section>=0,"Pro Tools text has no MARKERS section; export Session Info as Text with Include Markers");
    int header=section+1;while(header<lines.size()&&lines[header].trim().isEmpty())++header;
    need(header<lines.size(),"Missing Pro Tools marker header");auto names=juce::StringArray::fromTokens(lines[header],"\t","");
    int location=-1,name=-1,id=-1,reference=-1,units=-1;std::vector<int> types;
    for(int c=0;c<names.size();++c)
    {
        const auto key=canonical(names[c].trim());
        if(key=="location"||key=="position"){need(location<0,"Ambiguous marker position columns");location=c;}
        else if(key=="name"){need(name<0,"Duplicate marker name column");name=c;}
        else if(key=="#"||key=="number"){need(id<0,"Duplicate marker number column");id=c;}
        else if(key=="timereference"){need(reference<0,"Duplicate TIME REFERENCE column");reference=c;}
        else if(key=="units"){need(units<0,"Duplicate UNITS column");units=c;}
        else if(key=="type"||key=="markertype")types.push_back(c);
        else if(key=="end"||key=="length"||key=="duration")throw std::invalid_argument("Pro Tools range-marker table is unsupported; point markers only");
        else notices.add("Unused Pro Tools marker column: "+names[c]);
    }
    need(location>=0&&name>=0&&id>=0,"Unsupported Pro Tools marker table: require tab-separated #, LOCATION, NAME columns");
    const bool samples=o.unit==Unit::audioFrames;
    if(samples)
    {
        need(reference>=0&&units>=0,"Audio sample frames require TIME REFERENCE and UNITS columns; no fallback to LOCATION");
        int rateHeaders=0;auto numeric=o;numeric.unit=Unit::seconds;
        for(int n=0;n<section;++n)
        {
            poll(cancel);auto colon=lines[n].indexOfChar(':');
            if(colon>=0&&canonical(lines[n].substring(0,colon).trim())=="samplerate")
            {++rateHeaders;auto rate=parsePosition(lines[n].substring(colon+1),numeric,audioRate);need(rate.n==audioRate&&rate.d==1,"Pro Tools sample rate differs from programme; sample-reference resampling is unsupported");}
        }
        need(rateHeaders==1,"Sample-reference import requires one unambiguous SAMPLE RATE header");
        notices.add("Selected TIME REFERENCE / Samples at the matching session sample rate; LOCATION is display timecode and is not used. Confirm the sample-coordinate programme origin.");
    }
    else
    {
        if(reference>=0)notices.add("Unused Pro Tools marker column: TIME REFERENCE");
        if(units>=0)notices.add("Unused Pro Tools marker column: UNITS (describes TIME REFERENCE, not LOCATION)");
        notices.add("Selected LOCATION only; confirm its actual units/rate/origin independently of session headers. Display timecode can omit subframe precision; choose audio sample frames to use an available TIME REFERENCE / Samples field.");
    }
    notices.add("Only the marker table is imported. Unrelated session/track text is excluded.");
    std::vector<Row> result;std::set<int64_t> ids;
    for(int n=header+1;n<lines.size();++n)
    {
        poll(cancel);if(lines[n].trim().isEmpty())continue;
        // Known new session sections terminate the marker table. An unknown line
        // is an error instead of silently dropping possibly meaningful markers.
        const auto key=canonical(lines[n].trim());if(key=="tracklisting"||key=="clipslist"||key=="fileslist"||key=="onlinereferences"||key=="offlinereferences")break;
        auto cells=juce::StringArray::fromTokens(lines[n],"\t","");
        need(cells.size()==names.size(),"Pro Tools row "+juce::String(n+1)+": unsupported/ambiguous marker table structure");
        need(ids.insert(integer(cells[id])).second,"Duplicate Pro Tools marker number");for(auto t:types)need(point(cells[t].trim()),"Unsupported Pro Tools marker type");
        if(samples)need(cells[units].trim().equalsIgnoreCase("Samples"),"TIME REFERENCE requires Samples on every row; mixed or musical reference units are unsupported");
        Row row;row.line=n+1;row.position=cells[samples?reference:location];row.name=cells[name];label(row.name);result.push_back(std::move(row));
        need(result.size()<=maximumRows,"More than 10,000 marker rows");
    }
    return result;
}
}
juce::String contentHash(const std::string& s){return juce::SHA256(s.data(),s.size()).toHexString();}
Preview prepare(const std::string& input,const juce::String& basename,const Options& o,const SectionMap& before,const juce::String& runtime,const Cancel& cancel)
{
    need(before.valid()&&runtime.isNotEmpty(),"Confirm a valid programme before importing markers");poll(cancel);
    need(input.size()<=maximumBytes,"Marker text exceeds 8 MiB");need(!input.empty()&&input.find('\0')==std::string::npos,"Empty or unsupported marker encoding");
    need(juce::CharPointer_UTF8::isValidString(input.data(),int(input.size())),"Marker input must be UTF-8 (optional BOM); convert the export explicitly");
    need(o.confirmedPointMarkers&&o.confirmedUnits,"Confirm point-marker subset and actual exported units/origin");nominal(o);label(o.initialName);need(o.initialName.trim().isNotEmpty(),"Initial section name required");
    auto bytes=input;if(bytes.size()>=3&&bytes.substr(0,3)=="\xef\xbb\xbf")bytes.erase(0,3);
    Preview p;p.previousMap=before.digest();p.runtime=runtime;p.receipt.options=o;p.receipt.basename=basename;p.receipt.inputHash=contentHash(input);
    need(basename.length()<=256&&!basename.containsAnyOf("/\\\r\n"),"Use marker basename, not a path");
    auto& rows=p.receipt.rows;
    if(o.format==Format::nuendoCsv)rows=readNuendo(bytes,o,p.receipt.notices,cancel);
    else if(o.format==Format::proToolsText)rows=readProTools(bytes,o,before.sampleRate,p.receipt.notices,cancel);
    else
    {
        need(o.format==Format::pastedPositions,"Unsupported marker format");
        if(bytes.find(',')!=std::string::npos){need(o.unit==Unit::audioFrames||o.unit==Unit::videoFrames,"Comma paste lists require integer frame units");for(auto& c:bytes)if(c==',')c='\n';}
        auto lines=juce::StringArray::fromLines(juce::String::fromUTF8(bytes.data(),int(bytes.size())));
        for(int i=0;i<lines.size();++i){if(i+1==lines.size()&&lines[i].isEmpty())continue;Row r;r.line=i+1;r.position=lines[i];rows.push_back(std::move(r));}
    }
    need(!rows.empty()&&rows.size()<=maximumRows,"Require 1-10,000 point-marker rows");
    std::map<int64_t,size_t> points;int64_t previous=-1;
    for(size_t n=0;n<rows.size();++n)
    {
        poll(cancel);auto& r=rows[n];mapped(r,o,before);if(r.frame<previous)p.receipt.notices.addIfNotAlreadyThere("Rows reordered by audio-source position");previous=r.frame;
        auto found=points.find(r.frame);if(found!=points.end())
        {need(rows[found->second].name==r.name,"Conflicting marker names map to audio frame "+juce::String(r.frame));r.coalesced=true;p.receipt.notices.addIfNotAlreadyThere("Duplicate positions coalesced; all original rows retained");}
        else points.emplace(r.frame,n);
    }
    SectionMap after=before;after.sections.clear();std::map<std::pair<int64_t,int64_t>,ProgrammeSection> old;
    for(const auto& s:before.sections)old.emplace(std::make_pair(s.begin,s.end),s);
    std::vector<int64_t> boundaries{0};for(const auto& point:points)if(point.first>0&&point.first<before.frames)boundaries.push_back(point.first);boundaries.push_back(before.frames);
    need(boundaries.size()-1<=maximumRows,"More than 10,000 resulting sections");
    for(size_t n=0;n+1<boundaries.size();++n)
    {
        poll(cancel);const auto begin=boundaries[n],end=boundaries[n+1];ProgrammeSection s;s.begin=begin;s.end=end;
        auto prior=old.find({begin,end});if(prior!=old.end())s.id=prior->second.id;else{need(after.nextId<UINT64_MAX,"Section IDs exhausted");s.id=after.nextId++;}
        auto point=points.find(begin);s.name=point!=points.end()?rows[point->second].name:juce::String();
        if(s.name.isEmpty())s.name=begin==0?o.initialName:"Section at frame "+juce::String(begin);
        after.sections.push_back(s);
    }
    bool same=after.sections.size()==before.sections.size();if(same)for(size_t n=0;n<after.sections.size();++n)
    {const auto& a=after.sections[n];const auto& b=before.sections[n];same&=a.id==b.id&&a.begin==b.begin&&a.end==b.end&&a.name==b.name;}
    p.mapChanged=!same;if(!same){need(after.revision<UINT64_MAX,"Annotation revision exhausted");++after.revision;}
    need(after.valid(),"Invalid imported section coverage");p.receipt.imported=after;
    std::map<std::pair<int64_t,int64_t>,ProgrammeSection> current;for(const auto& s:after.sections)current.emplace(std::make_pair(s.begin,s.end),s);
    for(const auto& s:before.sections)
    {auto f=current.find({s.begin,s.end});if(f==current.end())p.diff.add("Remove "+juce::String(s.id)+" ["+juce::String(s.begin)+", "+juce::String(s.end)+") "+s.name);
     else if(s.name!=f->second.name)p.diff.add("Rename "+juce::String(s.id)+": "+s.name+" -> "+f->second.name);}
    for(const auto& s:after.sections)if(old.find({s.begin,s.end})==old.end())p.diff.add("Add "+juce::String(s.id)+" ["+juce::String(s.begin)+", "+juce::String(s.end)+") "+s.name);
    if(p.diff.isEmpty())p.diff.add("Section map unchanged; identical provenance is a no-op");
    need(juce::JSON::toString(encode(p.receipt),true).getNumBytesAsUTF8()<=3*1024*1024,"Marker provenance exceeds bounded session-core budget");poll(cancel);return p;
}
juce::String mapDescription(const SectionMap& s){return s.valid()?s.digest():juce::String();}
namespace
{
juce::var obj(){return juce::var(new juce::DynamicObject());}
void put(juce::var& v,const char* k,const juce::var& x){v.getDynamicObject()->setProperty(k,x);}
juce::var get(const juce::var& v,const char* k){need(v.isObject()&&v.hasProperty(k),juce::String("Missing marker provenance field: ")+k);return v[k];}
juce::String text(const juce::var& v,int limit=256){need(v.isString()&&v.toString().length()<=limit,"Invalid marker provenance string");return v.toString();}
int number(const juce::var& v,int lo,int hi){need(v.isInt()&&int(v)>=lo&&int(v)<=hi,"Invalid marker provenance number");return int(v);}
bool boolean(const juce::var& v){need(v.isBool(),"Invalid marker provenance boolean");return bool(v);}
const juce::Array<juce::var>& array(const juce::var& v,size_t limit){need(v.isArray()&&size_t(v.size())<=limit,"Invalid marker provenance rows");return *v.getArray();}
}
juce::var encode(const Receipt& r)
{
    auto v=obj();put(v,"version",version);put(v,"basename",r.basename);put(v,"inputSha256",r.inputHash);const auto& o=r.options;
    auto options=obj();put(options,"format",int(o.format));put(options,"unit",int(o.unit));put(options,"separator",int(o.separator));put(options,"header",o.header);
    put(options,"pointConfirmed",o.confirmedPointMarkers);put(options,"unitsConfirmed",o.confirmedUnits);put(options,"positionColumn",o.positionColumn);put(options,"nameColumn",o.nameColumn);put(options,"typeColumn",o.typeColumn);
    put(options,"rateNumerator",o.rateNumerator);put(options,"rateDenominator",o.rateDenominator);put(options,"origin",o.origin);put(options,"initialName",o.initialName);put(options,"rounding","floor-to-audio-sample");put(v,"options",options);
    juce::Array<juce::var> rows;for(const auto& r0:r.rows){juce::Array<juce::var> a;a.add(r0.line);a.add(r0.position);a.add(r0.name);a.add(juce::String(r0.frame));a.add(r0.seconds);a.add(r0.quantization);a.add(r0.coalesced);rows.add(a);}put(v,"rows",rows);
    juce::Array<juce::var> notices;for(const auto& n:r.notices)notices.add(n);put(v,"notices",notices);
    const auto& m=r.imported;auto map=obj();put(map,"input",m.inputIdentity);put(map,"rate",m.sampleRate);put(map,"frames",juce::String(m.frames));put(map,"revision",juce::String(m.revision));put(map,"nextId",juce::String(m.nextId));
    juce::Array<juce::var> sections;for(const auto& s:m.sections){juce::Array<juce::var> a;a.add(juce::String(s.id));a.add(juce::String(s.begin));a.add(juce::String(s.end));a.add(s.name);sections.add(a);}put(map,"sections",sections);put(v,"imported",map);return v;
}
Receipt decode(const juce::var& v)
{
    need(text(get(v,"version"))==version,"Unsupported marker provenance version");Receipt r;r.basename=text(get(v,"basename"));r.inputHash=text(get(v,"inputSha256"),64);
    need(r.inputHash.length()==64&&r.inputHash.containsOnly("0123456789abcdef")&&!r.basename.containsAnyOf("/\\\r\n"),"Invalid marker provenance identity");
    auto options=get(v,"options");auto& o=r.options;o.format=Format(number(get(options,"format"),0,2));o.unit=Unit(number(get(options,"unit"),0,3));o.separator=char(number(get(options,"separator"),0,127));
    need(o.separator==','||o.separator==';'||o.separator=='\t',"Invalid saved marker delimiter");o.header=boolean(get(options,"header"));o.confirmedPointMarkers=boolean(get(options,"pointConfirmed"));o.confirmedUnits=boolean(get(options,"unitsConfirmed"));need(o.confirmedPointMarkers&&o.confirmedUnits,"Unconfirmed marker receipt");
    o.positionColumn=number(get(options,"positionColumn"),0,127);o.nameColumn=number(get(options,"nameColumn"),-1,127);o.typeColumn=number(get(options,"typeColumn"),-1,127);
    o.rateNumerator=number(get(options,"rateNumerator"),1,30000);o.rateDenominator=number(get(options,"rateDenominator"),1,1001);nominal(o);o.origin=text(get(options,"origin"),128);o.initialName=text(get(options,"initialName"),128);label(o.initialName);need(o.initialName.trim().isNotEmpty()&&text(get(options,"rounding"))=="floor-to-audio-sample","Invalid marker mapping convention");
    auto map=get(v,"imported");auto& m=r.imported;m.inputIdentity=text(get(map,"input"),64);m.sampleRate=number(get(map,"rate"),44100,192000);m.frames=integer(text(get(map,"frames")));m.revision=uint64_t(integer(text(get(map,"revision"))));m.nextId=uint64_t(integer(text(get(map,"nextId"))));
    for(const auto& row:array(get(map,"sections"),maximumRows)){const auto& a=array(row,4);need(a.size()==4,"Invalid imported section");m.sections.push_back({uint64_t(integer(text(a[0]))),integer(text(a[1])),integer(text(a[2])),text(a[3],128)});}
    need(m.valid(),"Invalid saved imported map");std::map<int64_t,juce::String> points;
    for(const auto& row:array(get(v,"rows"),maximumRows))
    {
        const auto& a=array(row,7);need(a.size()==7,"Invalid marker row");Row x;x.line=number(a[0],1,1000000);x.position=text(a[1],128);x.name=text(a[2],128);label(x.name);x.coalesced=boolean(a[6]);mapped(x,o,m);
        need(x.frame==integer(text(a[3]))&&x.seconds==text(a[4],128)&&x.quantization==text(a[5],128),"Saved marker mapping contradicts original position");auto f=points.find(x.frame);
        need(x.coalesced==(f!=points.end())&&(f==points.end()||f->second==x.name),"Invalid saved duplicate handling");points.emplace(x.frame,x.name);r.rows.push_back(x);
    }
    need(!r.rows.empty(),"Empty marker receipt");size_t internal=0;for(const auto& p:points)if(p.first>0&&p.first<m.frames)++internal;need(m.sections.size()==internal+1,"Marker rows and imported sections disagree");
    for(const auto& s:m.sections)
    {auto f=points.find(s.begin);need(s.begin==0||f!=points.end(),"Imported boundary missing marker");auto name=f!=points.end()?f->second:juce::String();if(name.isEmpty())name=s.begin==0?o.initialName:"Section at frame "+juce::String(s.begin);need(name==s.name,"Imported marker name mismatch");}
    for(const auto& n:array(get(v,"notices"),256))r.notices.add(text(n,1024));return r;
}
juce::String Receipt::identity()const{return contentHash(juce::JSON::toString(encode(*this),true).toStdString());}
}
