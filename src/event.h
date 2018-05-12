#ifndef EVENT_H
#define EVENT_H

#include <functional>
#include <vector>

#include "hit.h"

struct hit
{
    float r, phi, z;
};

struct track
{
    float pt, eta, phi;
    std::vector<hit> hits, seed;
};

struct event
{
    std::vector<hit> hits;
    std::vector<track> tracks;
};


#endif // EVENT_H
