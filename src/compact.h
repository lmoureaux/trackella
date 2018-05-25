#ifndef COMPACT_H
#define COMPACT_H

#include <cstdint>

#include "geometry.h"
#include "event.h"

const constexpr float pi = 3.141592653589793f;

/**
 * \brief Converts the compact representation to radians.
 */
constexpr float compact_to_radians(std::int16_t compact)
{
    return compact / float(1 << 15) * pi;
}

/**
 * \brief Converts radians to the compact representation.
 */
constexpr std::int16_t radians_to_compact(float radians)
{
    return radians * float(1 << 15) / pi;
}


/**
 * \brief Converts the compact representation to a length (in cm).
 */
template<class T> constexpr float compact_to_length(T compact)
{
    return compact / float(1 << 14);
}

/**
 * \brief Converts a length (in cm) to the compact representation.
 */
template<class T> constexpr T length_to_compact(float length)
{
    return length * float(1 << 14);
}

/**
 * \brief Holds compact information about the beam spot.
 */
#pragma pack(push,4)
struct compact_beam_spot
{
    std::int32_t r;   ///< \brief Radius (cm)
    std::int32_t z;   ///< \brief Position along the \c z axis
    std::int16_t phi; ///< \brief Azimutal angle (compact representation)
};
static_assert(sizeof(compact_beam_spot) == 12, "compact_beam_spot size");
#pragma pack(pop)

/**
 * \brief Holds compact information about a pixel hit in the barrel.
 *
 * It is assumed that the pixel layer is known from elsewhere.
 */
struct compact_hit
{
    /**
     * \brief Difference in \c r with respect to the layer's average.
     *
     * Encoding: 1 unit = 2cm / 2^15 = 6.1um.
     *
     * Range: -2cm to 2cm.
     */
    std::int16_t dr;

    /**
     * \brief Azimutal angle \c phi.
     *
     * Encoding: 1 unit = pi / 2^15 rad.
     *
     * Range: -pi to pi.
     */
    std::int16_t phi;

    /**
     * \brief Longitudinal component.
     *
     * Encoding: 1 unit = 2cm / 2^15 = 6.1um.
     *
     * Range: -1.3km to 1.3km.
     */
    std::int32_t z;

    /**
     * \brief Constructor
     */
    explicit compact_hit(const hit &h, int layer) :
        dr(length_to_compact<std::int16_t>(h.r - geom::pixel_barrel_radius[layer])),
        phi(radians_to_compact(h.phi)),
        z(length_to_compact<int>(h.z))
    {}
};
static_assert(sizeof(compact_hit) == 8, "compact_hit size");

#endif // COMPACT_H
