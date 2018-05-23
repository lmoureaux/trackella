#ifndef FAST_SINCOS_H
#define FAST_SINCOS_H

#include <cmath>
#include <cstdint>

#include "compact.h"

class fast_sincos
{
public:
    using sincos_type = int;
    using angle_type = std::int16_t;

private:
    sincos_type _sin, _cos;
    angle_type _last;

public:
    explicit fast_sincos(angle_type start) :
        _sin(std::sin(compact_to_radians(start)) * (1 << 8)),
        _cos(std::cos(compact_to_radians(start)) * (1 << 8)),
        _last(start)
    {}

    void sync()
    {
        float last = compact_to_radians(_last);
        _sin = std::sin(last) * (1 << 8);
        _cos = std::cos(last) * (1 << 8);
    }

    void sync(std::int16_t angle)
    {
        _last = angle;
        sync();
    }

    void step(angle_type angle)
    {
        angle_type da = (angle - _last) >> 8;
        sincos_type oldsin = _sin;
        sincos_type oldcos = _cos;
        // First order Taylor development
        auto tmp = 402 * da;
        _sin += (oldcos * tmp) >> (7 + 7);
        _cos -= (oldsin * tmp) >> (7 + 7);
        _last = angle;
    }

    angle_type angle() const
    {
        return _last;
    }

    sincos_type sin() const
    {
        return _sin;
    }

    sincos_type cos() const
    {
        return _cos;
    }

    template<class U>
    U sin_times(U value) const
    {
        return _sin * value >> 8;
    }

    template<class U>
    U cos_times (U value) const
    {
        return _cos * value >> 8;
    }
};

class fast_float_sincos
{
public:
    using sincos_type = float;
    using angle_type = float;

private:
    sincos_type _sin, _cos;
    angle_type _last;

public:
    explicit fast_float_sincos(angle_type start) :
        _sin(std::sin(start)),
        _cos(std::cos(start)),
        _last(start)
    {}

    void sync()
    {
        _sin = std::sin(_last);
        _cos = std::cos(_last);
    }

    void sync(float angle)
    {
        _last = angle;
        sync();
    }

    void step(angle_type angle)
    {
        angle_type da = angle - _last;
        sincos_type oldsin = _sin;
        sincos_type oldcos = _cos;
        // First order Taylor development
        _sin += oldcos * da;
        _cos -= oldsin * da;
        _last = angle;
    }

    angle_type angle() const
    {
        return _last;
    }

    sincos_type sin() const
    {
        return _sin;
    }

    sincos_type cos() const
    {
        return _cos;
    }

    template<class U>
    U sin_times(U value) const
    {
        return _sin * value;
    }

    template<class U>
    U cos_times (U value) const
    {
        return _cos * value;
    }
};

#endif // FAST_SINCOS_H
