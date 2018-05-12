#include "eventreader.h"

#include <algorithm>
#include <iostream>

#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>

struct event_reader::data
{
    TFile input;
    TTreeReader reader;

    TTreeReaderArray<float> trk_pt;
    TTreeReaderArray<float> trk_eta;
    TTreeReaderArray<float> trk_phi;
    TTreeReaderArray<int>   trk_hit_n;
    TTreeReaderArray<float> trk_hit_globalPos_x;
    TTreeReaderArray<float> trk_hit_globalPos_y;
    TTreeReaderArray<float> trk_hit_globalPos_z;
    TTreeReaderArray<int>   trk_seed_n;
    TTreeReaderArray<float> trk_seed_globalPos_x;
    TTreeReaderArray<float> trk_seed_globalPos_y;
    TTreeReaderArray<float> trk_seed_globalPos_z;

    data(const std::string &filename);
};

event_reader::data::data(const std::string &filename) :
    input(filename.c_str(), "OPEN"),
    reader("TrackTree/tree", &input),
    trk_pt(reader, "trk_pt"),
    trk_eta(reader, "trk_eta"),
    trk_phi(reader, "trk_phi"),
    trk_hit_n(reader, "trk_hit_n"),
    trk_hit_globalPos_x(reader, "trk_hit_globalPos_x"),
    trk_hit_globalPos_y(reader, "trk_hit_globalPos_y"),
    trk_hit_globalPos_z(reader, "trk_hit_globalPos_z"),
    trk_seed_n(reader, "trk_seed_n"),
    trk_seed_globalPos_x(reader, "trk_seed_globalPos_x"),
    trk_seed_globalPos_y(reader, "trk_seed_globalPos_y"),
    trk_seed_globalPos_z(reader, "trk_seed_globalPos_z")
{
}

event_reader::event_reader(const std::string &filename) :
    _d(std::make_unique<data>(filename))
{
}

event_reader::~event_reader()
{
}

namespace /* anonymous */
{
bool hit_ordering(const hit &lhs, const hit &rhs)
{
    return lhs.phi < rhs.phi && lhs.z < rhs.z && lhs.phi < rhs.phi;
}

bool hit_compare(const hit &lhs, const hit &rhs)
{
    return lhs.phi < rhs.phi && lhs.z < rhs.z && lhs.phi < rhs.phi;
}
} // namespace anonymous

std::unique_ptr<event> event_reader::get()
{
    std::vector<hit> hits;
    hits.reserve(10 * _d->trk_pt.GetSize());

    std::vector<track> tracks;
    tracks.reserve(_d->trk_pt.GetSize());

    std::size_t ihit = 0, iseed = 0;
    for (std::size_t itrk = 0; itrk < _d->trk_pt.GetSize(); ++itrk) {
        // Hits
        int hit_count = _d->trk_hit_n[itrk];
        std::vector<hit> track_hits;
        track_hits.reserve(hit_count);
        for (std::size_t i = 0; i < hit_count; ++i) {
            float x = _d->trk_hit_globalPos_x[ihit + i];
            float y = _d->trk_hit_globalPos_y[ihit + i];
            float z = _d->trk_hit_globalPos_z[ihit + i];

            float r = std::sqrt(x * x + y * y);
            float phi = std::atan2(y, x);

            hits.push_back({ r, phi, z });
            track_hits.push_back({ r, phi, z });
        }
        ihit += hit_count;

        // Seeds
        int seed_count = _d->trk_seed_n[itrk];
        std::vector<hit> track_seed;
        track_seed.reserve(seed_count);
        for (std::size_t i = 0; i < seed_count; ++i) {
            float x = _d->trk_seed_globalPos_x[iseed + i];
            float y = _d->trk_seed_globalPos_y[iseed + i];
            float z = _d->trk_seed_globalPos_z[iseed + i];

            float r = std::sqrt(x * x + y * y);
            float phi = std::atan2(y, x);

            track_seed.push_back(hit{ r, phi, z });
        }
        iseed += seed_count;

        float pt = _d->trk_pt[itrk];
        float eta = _d->trk_eta[itrk];
        float phi = _d->trk_phi[itrk];

        tracks.push_back({ pt, eta, phi, track_hits, track_seed });
    }

    std::unique_ptr<event> e = std::make_unique<event>();
    e->hits = std::move(hits);
    e->tracks = std::move(tracks);
    return e;
}

bool event_reader::next()
{
    return _d->reader.Next();
}
