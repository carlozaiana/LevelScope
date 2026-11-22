"""Build a separate, instrumented test object; never change frozen product inputs."""
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys

ui = Path(__file__).resolve().parents[3]
candidate = ui/'build/m20-competing-boundary-ui-2026-11-21/candidate-01'
out = Path(sys.argv[1]).resolve()
assert not out.exists()
out.mkdir(parents=True)
compile_dir = ui/'build/m20-ui-integration/compile'
target = 'LevelScopeRegionUITests'
identity = json.loads((candidate/'identity.json').read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
for name, value in identity['source_sha256'].items():
    assert sha(candidate/'source'/name) == value, name
link = shlex.split((candidate/'configuration'/(target+'-link.txt')).read_text())
link[link.index('-o')+1] = str(out/target)
subprocess.run(link, cwd=compile_dir, check=True)
assert sha(out/target) == identity['artifacts_sha256']['bin/'+target]
source = (candidate/'source/Tests/RegionUIIntegrationTests.cpp').read_text()
source = source.replace('struct Heartbeat:private juce::Timer', '''
juce::Array<juce::var> diagnosticEvents;
void diagnosticMark(const char* name) {
    auto row=juce::var(new juce::DynamicObject());
    row.getDynamicObject()->setProperty("mark",name);
    row.getDynamicObject()->setProperty("time_ms",juce::Time::getMillisecondCounterHiRes());
    diagnosticEvents.add(row);
}
struct Heartbeat:private juce::Timer''')
source = source.replace('gaps.push_back(now-previous);previous=now;', '''
auto row=juce::var(new juce::DynamicObject());
row.getDynamicObject()->setProperty("gap_ms",now-previous);
row.getDynamicObject()->setProperty("time_ms",now);diagnosticEvents.add(row);
gaps.push_back(now-previous);previous=now;''')
source = source.replace('stopTimer();std::sort(gaps.begin(),gaps.end());', '''
stopTimer();auto diagnostic=juce::var(new juce::DynamicObject());
diagnostic.getDynamicObject()->setProperty("heartbeat_diagnostic",diagnosticEvents);
cases.add(diagnostic);std::sort(gaps.begin(),gaps.end());''')
event = (candidate/'source/Tests/EventUIIntegrationChecks.h').read_text()
replacements = [
    ('Heartbeat heartbeat;heartbeat.begin();', 'diagnosticMark("proposal");Heartbeat heartbeat;heartbeat.begin();'),
    ('"native Propose starts async event assessment after request capture");waitJob(ui);', '"native Propose starts async event assessment after request capture");waitJob(ui);diagnosticMark("review-setup");'),
    ('for(size_t n=0;n<review->attempts.size();++n)', 'diagnosticMark("browse-attempts");for(size_t n=0;n<review->attempts.size();++n)'),
    ('const auto restore=s;', 'diagnosticMark("guard-manipulations");const auto restore=s;'),
    ('const auto selected=s.sectionProposal;', 'diagnosticMark("selection-inspection");const auto selected=s.sectionProposal;'),
    ('imageOf(*panel,root.getChildFile("selected-partial-review.png"));', 'diagnosticMark("snapshot-review");imageOf(*panel,root.getChildFile("selected-partial-review.png"));diagnosticMark("snapshot-main");'),
    ('std::vector<double> switches;', 'diagnosticMark("switch-loop");std::vector<double> switches;'),
    ('heartbeat.finish();', 'diagnosticMark("finish");heartbeat.finish();'),
]
for old, new in replacements:
    assert event.count(old) == 1, old
    event = event.replace(old, new)
if '--callbacks' in sys.argv:
    marker = 'diagnosticMark("switch-loop");'
    event = event.replace(marker, '''
    diagnosticMark("isolated-read-only-callbacks");
    const auto timeCallback=[&](const char* name,auto callback) {
        std::vector<double> times;
        for(int i=0;i<10;++i){const auto begin=juce::Time::getMillisecondCounterHiRes();callback();times.push_back(juce::Time::getMillisecondCounterHiRes()-begin);}
        std::sort(times.begin(),times.end());auto row=juce::var(new juce::DynamicObject());
        row.getDynamicObject()->setProperty("callback",name);row.getDynamicObject()->setProperty("p95_ms",times.back());cases.add(row);
    };
    timeCallback("stale-check",[&]{(void)s.sectionProposalStale();});
    timeCallback("can-apply-check",[&]{(void)s.canApplySectionProposal();});
    timeCallback("report-status",[&]{view().refreshStatus();});
    timeCallback("panel-selection",[&]{panel->refreshSelection();});
    timeCallback("timeline-refresh",[&]{Access::refreshTimeline(ui);});
    diagnosticMark("switch-loop");''')
(out/'RegionUIIntegrationTests.cpp').write_text(source)
(out/'EventUIIntegrationChecks.h').write_text(event)
flags = {}
for line in (candidate/'configuration'/(target+'-flags.make')).read_text().splitlines():
    if ' = ' in line:
        key, value = line.split(' = ',1); flags[key] = shlex.split(value)
args = flags['CXX_DEFINES']+flags['CXX_INCLUDES']+flags['CXX_FLAGS']
args = [a.replace('-I'+str(ui/'Source'),'-I'+str(candidate/'source/Source')) for a in args]
obj = out/'diagnostic.o'
subprocess.run(['/usr/bin/c++',*args,'-I'+str(candidate/'source/Tests'),'-c',str(out/'RegionUIIntegrationTests.cpp'),'-o',str(obj)],check=True)
link[link.index('-o')+1] = str(out/'heartbeat-diagnostic')
index = next(i for i, a in enumerate(link) if a.endswith('/RegionUIIntegrationTests.cpp.o'))
link[index] = str(obj)
subprocess.run(link,cwd=compile_dir,check=True)
(out/'identity.json').write_text(json.dumps(dict(candidate_identity=sha(candidate/'identity.json'),
    original_relinked=sha(out/target),diagnostic=sha(out/'heartbeat-diagnostic'),
    test_source=sha(out/'RegionUIIntegrationTests.cpp'),event_header=sha(out/'EventUIIntegrationChecks.h'),
    builder=sha(Path(__file__)),production_changes=False,thresholds_changed=False),indent=2)+'\n')
