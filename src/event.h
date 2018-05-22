#ifndef EVENT_H
#define EVENT_H

#include <functional>
#include <vector>

struct hit
{
    float r, phi, z;
};

struct beam_spot
{
    float r, phi, z;
};

struct track
{
    float pt, eta, phi, b0, z0;
    std::vector<hit> hits, seed;
};

struct event
{
    beam_spot bs;
    std::vector<hit> hits;
    std::vector<track> tracks;
    int nvtx;
};


#endif // EVENT_H
