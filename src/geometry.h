#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <array>

namespace geom
{
    const constexpr std::array<float, 4> pixel_barrel_radius = {
        3.0f,
        6.8f,
        11.0f,
        16.0f
    };
} // namespace geom

#endif // GEOMETRY_H
