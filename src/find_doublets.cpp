#include <array>
#include <iostream>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include "doublet_finder.h"
#include "eventreader.h"
#include "hitutils.h"

int main(int, char **)
{
    event_reader in("~lmoureau/data/v3.root");
    long long i = 0;
    while (in.next()) {
        std::cout << "==== Next event ====" << std::endl;

        std::unique_ptr<event> e = in.get();
        std::sort(e->hits.begin(), e->hits.end(), hit_less_than);
        e->hits.erase(std::unique(e->hits.begin(), e->hits.end(), hit_equal),
                      e->hits.end());

        std::vector<hit> pb_hits;
        std::copy_if(e->hits.begin(), e->hits.end(), std::back_inserter(pb_hits), hit_is_pixel_barrel);

        std::array<std::vector<hit>, 4> pb_hits_per_layer;
        for (const hit &h : pb_hits) {
            pb_hits_per_layer[hit_pixel_barrel_layer(h)].push_back(h);
        }
        std::cout << "Hits in 1st layer: " << pb_hits_per_layer[0].size() << std::endl;
        std::cout << "Hits in 2nd layer: " << pb_hits_per_layer[1].size() << std::endl;

        std::cout << "Making doublets..." << std::endl;

        std::vector<compact_pb_hit> layer1;
        layer1.reserve(pb_hits_per_layer[0].size());
        for (const hit &h : pb_hits_per_layer[0]) {
            assert(std::abs(h.r - 3) < 2);
            layer1.push_back({
                length_to_compact<std::int16_t>(h.r - 3),
                radians_to_compact(h.phi),
                length_to_compact<std::int32_t>(h.z)
            });
        }

        std::vector<compact_pb_hit> layer2;
        layer2.reserve(pb_hits_per_layer[1].size());
        for (const hit &h : pb_hits_per_layer[1]) {
            assert(std::abs(h.r - 6.8) < 2);
            layer2.push_back({
                length_to_compact<std::int16_t>(h.r - 6.8),
                radians_to_compact(h.phi),
                length_to_compact<std::int32_t>(h.z)
            });
        }

        pb_doublet_finder finder;
        finder.set_beam_spot({
            e->bs.r,
            e->bs.phi,
            length_to_compact<std::int32_t>(e->bs.z)
        });
        finder.set_hits(layer1, layer2);

        finder.start();

        std::vector<pb_doublet_finder::doublet> doublets;
        finder.get_doublets(doublets);
        std::cout << "Doublets: "
                  << doublets.size()
                  << "; factor: "
                  << layer1.size() * layer2.size() / doublets.size() << std::endl;
    }
}
