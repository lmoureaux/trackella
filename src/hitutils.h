#ifndef HITUTILS_H
#define HITUTILS_H

#include "event.h"

/**
 * \brief Hit compare function
 */
inline bool hit_less_than(const hit &a, const hit &b)
{
    return a.z < b.z && a.phi < b.phi && a.r < b.r;
}

/**
 * \brief Hit compare function
 */
inline bool hit_equal(const hit &a, const hit &b)
{
    return a.z == b.z && a.phi == b.phi && a.r == b.r;
}

/**
 * \brief Returns \c true if \c h is from the pixel barrel
 */
inline bool hit_is_pixel_barrel(const hit &h)
{
    return h.r < 20.f && std::abs(h.z) < 28.f;
}

/**
 * \brief Returns the pixel layer the hit originates from (0-3)
 * \pre `hit_is_pixel_barrel(h)` is \c true.
 */
inline int hit_pixel_barrel_layer(const hit &h)
{
    return h.r * 0.2f;
}

#endif // HITUTILS_H
