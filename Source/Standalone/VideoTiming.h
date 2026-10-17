#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace levelscope::standalone::video
{
// All timeline arithmetic is exact. Overflow is an error, never a wrapped timestamp.
// Quantization uses floor (towards -infinity), including negative offsets.
class Time
{
public:
    int64_t numerator = 0, denominator = 1;
    Time() = default;
    Time (int64_t n, int64_t d) { *this = wide (n, d); }
    static Time wide (__int128 n, __int128 d)
    {
        if (d <= 0) throw std::invalid_argument ("Invalid timebase");
        auto a = n < 0 ? -n : n, b = d;
        while (b != 0) { const auto r = a % b; a = b; b = r; }
        if (a != 0) { n /= a; d /= a; }
        if (n < INT64_MIN || n > INT64_MAX || d > INT64_MAX)
            throw std::overflow_error ("Timestamp overflow");
        Time result; result.numerator = static_cast<int64_t> (n); result.denominator = static_cast<int64_t> (d);
        return result;
    }
    Time operator+ (Time b) const { return wide ((__int128) numerator * b.denominator + (__int128) b.numerator * denominator,
                                                  (__int128) denominator * b.denominator); }
    Time operator- (Time b) const { return wide ((__int128) numerator * b.denominator - (__int128) b.numerator * denominator,
                                                  (__int128) denominator * b.denominator); }
    Time operator* (Time b) const { return wide ((__int128) numerator * b.numerator, (__int128) denominator * b.denominator); }
    bool operator< (Time b) const { return (__int128) numerator * b.denominator < (__int128) b.numerator * denominator; }
    bool operator== (Time b) const { return numerator == b.numerator && denominator == b.denominator; }
    bool operator!= (Time b) const { return ! (*this == b); }
    int64_t floorTicks (int64_t scale) const
    {
        if (scale <= 0) throw std::invalid_argument ("Invalid clock rate");
        const auto n = (__int128) numerator * scale;
        auto q = n / denominator;
        if (n % denominator < 0) --q;
        if (q < INT64_MIN || q > INT64_MAX) throw std::overflow_error ("Clock conversion overflow");
        return static_cast<int64_t> (q);
    }
    double seconds() const { return static_cast<double> (numerator) / static_cast<double> (denominator); }
    std::string text() const { return std::to_string (numerator) + "/" + std::to_string (denominator); }
    static Time parse (const std::string& text)
    {
        // Signed decimal seconds (up to nanoseconds), or a signed integer rational.
        const auto slash = text.find ('/');
        const auto integer = [] (const std::string& s)
        {
            if (s.empty()) throw std::invalid_argument ("Enter signed seconds or numerator/denominator");
            size_t end = 0; const auto n = std::stoll (s, &end);
            if (end != s.size()) throw std::invalid_argument ("Invalid timestamp");
            return n;
        };
        if (slash != std::string::npos) return Time (integer (text.substr (0, slash)), integer (text.substr (slash + 1)));
        const auto dot = text.find ('.');
        if (dot == std::string::npos) return Time (integer (text), 1);
        auto whole = text.substr (0, dot), fraction = text.substr (dot + 1);
        if (fraction.empty() || fraction.size() > 9 || fraction.find_first_not_of ("0123456789") != std::string::npos)
            throw std::invalid_argument ("Use at most 9 decimal places");
        const bool negative = !whole.empty() && whole[0] == '-';
        if (whole.empty() || whole == "+" || whole == "-") whole += "0";
        int64_t scale = 1; for (size_t i = 0; i < fraction.size(); ++i) scale *= 10;
        return wide ((__int128) integer (whole) * scale + (negative ? -1 : 1) * integer (fraction), scale);
    }
};

struct FrameTime { Time pts, end; };

// Half-open coverage; explicit PTS preserve VFR and gaps. No nominal-rate inference.
inline std::optional<size_t> frameAt (const std::vector<FrameTime>& frames, Time time)
{
    const auto next = std::upper_bound (frames.begin(), frames.end(), time,
        [] (Time t, const FrameTime& f) { return t < f.pts; });
    if (next == frames.begin()) return {};
    const auto frame = std::prev (next);
    if (!(time < frame->end)) return {};
    return static_cast<size_t> (frame - frames.begin());
}

struct Timecode
{
    int64_t startFrame = 0;
    Time origin, rate { 30, 1 };
    bool dropFrame = false, wraps24Hours = false;
};

inline std::string timecodeLabel (int64_t frame, Time rate, bool drop, bool wrap = false)
{
    const auto rounded = ((__int128) rate.numerator + rate.denominator / 2) / rate.denominator;
    if (rounded < 1 || rounded > 120) throw std::invalid_argument ("Unsupported timecode rate");
    const auto nominal = static_cast<int64_t> (rounded);
    if (drop && rate != Time (30000, 1001) && rate != Time (60000, 1001))
        throw std::invalid_argument ("Drop-frame requires 30000/1001 or 60000/1001");
    const bool negative = frame < 0;
    __int128 labels = negative ? -(__int128) frame : frame;
    if (drop)
    {
        const int64_t dropped = nominal == 60 ? 4 : 2;
        const int64_t tenMinutes = nominal * 600 - 9 * dropped, minute = nominal * 60 - dropped;
        const auto remainder = labels % tenMinutes;
        labels += (labels / tenMinutes) * 9 * dropped;
        if (remainder >= dropped) labels += ((remainder - dropped) / minute) * dropped;
    }
    const auto frames = static_cast<int> (labels % nominal); labels /= nominal;
    const auto seconds = static_cast<int> (labels % 60); labels /= 60;
    const auto minutes = static_cast<int> (labels % 60); labels /= 60;
    if (wrap) labels %= 24;
    if (labels > INT64_MAX) throw std::overflow_error ("Timecode hours overflow");
    auto pad = [] (int64_t n) { auto s = std::to_string (n); return s.size() < 2 ? "0" + s : s; };
    return std::string (negative ? "-" : "") + pad (static_cast<int64_t> (labels)) + ":" + pad (minutes)
         + ":" + pad (seconds) + (drop ? ";" : ":") + pad (frames);
}
}
