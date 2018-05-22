#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <array>

#include "compact.h"

namespace geom
{
    namespace real
    {
        std::array<float, 4> pixel_barrel_radius = { 3.0f, 6.8f, 11.0f, 16.0f };
    }
    namespace compact
    {
        std::array<int, 4> pixel_barrel_radius = {
            length_to_compact<int>(real::pixel_barrel_radius[0]),
            length_to_compact<int>(real::pixel_barrel_radius[1]),
            length_to_compact<int>(real::pixel_barrel_radius[2]),
            length_to_compact<int>(real::pixel_barrel_radius[3])
        };
    }
} // namespace geom

#endif // GEOMETRY_H
