#include "TraceBuffer.h"
#include <iostream>
#include <thread>
#include <stdexcept>
using namespace levelscope::paintdiagnostic;
static void expect(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
int main()
{
    try
    {
        Buffer<3> small;
        for(uint64_t i=0;i<5;++i) { Event e;e.hostNs=i;e.picture=static_cast<int64_t>(i);small.push(e); }
        expect(small.size()==3 && small.lost()==2,"overflow must preserve first events and count loss");
        expect(small.at(0).picture==0 && small.at(2).picture==2,"overflow must not overwrite causal history");
        Buffer<1024> ui,display;
        std::thread producer([&] { for(uint64_t i=0;i<1000;++i) { Event e;e.kind=Kind::displayTick;e.hostNs=i;display.push(e); } });
        for(uint64_t i=0;i<1000;++i) { Event e;e.kind=Kind::selection;e.picture=i;e.epoch=7;ui.push(e); }
        producer.join();
        expect(ui.size()==1000 && display.size()==1000 && !ui.lost() && !display.lost(),"independent producer histories survive");
        for(size_t i=0;i<1000;++i) expect(ui.at(i).picture==static_cast<int64_t>(i) && display.at(i).hostNs==i && ui.at(i).epoch==7,"ordered lossless fields");
        std::cout<<"PASS bounded trace, overflow preservation, independent producers and ordering\n";
        return 0;
    }
    catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
