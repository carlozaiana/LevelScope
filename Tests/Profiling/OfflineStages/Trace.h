#pragma once
#include <array>
#include <chrono>

// Linked only into the private profiling executable. One worker per operation.
// No timing calls inside the sample, channel, tap or peak-observer loops.
namespace lsprofile
{
using Clock = std::chrono::steady_clock;
enum Stage { verified, copy, copyHash, decodeCheck, validationHash, read, gain, wave, meter, publish, finish, passTotal, count };
inline constexpr const char* names[] = { "verification_total", "copy_flush", "first_hash", "payload_scan", "validation_hash",
    "read_blocks", "processing", "waveform", "meter_push", "progress_publish", "meter_finish", "pass_total" };
inline constexpr const char* phases[] = { "job", "attribution", "local", "combined" };
using Timings = std::array<std::array<double,count>,4>;
inline thread_local Timings values{};
inline thread_local size_t phase = 0;
inline double ms(Clock::time_point a, Clock::time_point b) { return std::chrono::duration<double,std::milli>(b-a).count(); }
struct Lap
{
    Clock::time_point last=Clock::now();
    void mark(Stage stage) { const auto now=Clock::now(); values[phase][stage]+=ms(last,now); last=now; }
};
struct Span
{
    Stage stage;
    Clock::time_point start=Clock::now();
    ~Span() { values[phase][stage]+=ms(start,Clock::now()); }
};
struct Pass
{
    size_t previous;
    Clock::time_point start=Clock::now();
    explicit Pass(size_t which):previous(phase) { phase=which; }
    ~Pass() { values[phase][passTotal]+=ms(start,Clock::now()); phase=previous; }
};
}
