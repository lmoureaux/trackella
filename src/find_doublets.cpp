#include <array>
#include <chrono>
#include <iostream>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include "doublet_finder.h"
#include "eventreader.h"
#include "hitutils.h"

float deltaphi(float phi1, float phi2)
{
    float delta = phi1 - phi2;
    if (delta > pi) {
        delta -= 2 * pi;
    } else if (delta <= -pi) {
        delta += 2 * pi;
    }
    return delta;
}

int main(int, char **)
{
    std::chrono::duration<double> formatting;
    long long formatted_hits = 0;

    std::chrono::duration<double> finding;
    long long doublets_found = 0;

    TFile out("doublets.root", "RECREATE");

    TH1D doublet_phi1("doublet_phi1", ";phi1;count", 50, -pi, pi);
    TH1D doublet_phi2("doublet_phi2", ";phi2;count", 50, -pi, pi);
    TH1D doublet_phi2_phi1("doublet_phi2_phi1", ";phi2 - phi1;count", 50, -0.05, 0.05);
    TH1D doublet_z1("doublet_z1", ";z1;count", 50, -30, 30);
    TH1D doublet_z2("doublet_z2", ";z2;count", 50, -30, 30);
    TH1D doublet_z0("doublet_z0", ";z0;count", 50, -15, 15);
    TH1D doublet_b0("doublet_b0", ";b0;count", 50, -0.05, 0.05);

    event_reader in("~lmoureau/data/v3.root");
    long long i = 0;
    while (in.next()) {
        std::cout << "==== Next event ====" << std::endl;

        std::unique_ptr<event> e = in.get();
        std::sort(e->hits.begin(), e->hits.end(), hit_less_than);
        e->hits.erase(std::unique(e->hits.begin(), e->hits.end(), hit_equal),
                      e->hits.end());

        std::array<std::vector<hit>, 4> pb_hits_per_layer;
        for (const hit &h : e->hits) {
            if (hit_is_pixel_barrel(h)) {
                pb_hits_per_layer[hit_pixel_barrel_layer(h)].push_back(h);
            }
        }
        std::cout << "Hits in 1st layer: " << pb_hits_per_layer[0].size() << std::endl;
        std::cout << "Hits in 2nd layer: " << pb_hits_per_layer[1].size() << std::endl;

        std::cout << "Making doublets..." << std::endl;

        auto start = std::chrono::high_resolution_clock::now();

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
        std::sort(layer1.begin(),
                  layer1.end(),
                  [](const compact_pb_hit &a, const compact_pb_hit &b) {
                      return a.phi < b.phi;
                  });

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
        std::sort(layer2.begin(),
                  layer2.end(),
                  [](const compact_pb_hit &a, const compact_pb_hit &b) {
                      return a.phi < b.phi;
                  });

        pb_doublet_finder finder;
        finder.set_beam_spot({
            e->bs.r,
            e->bs.phi,
            length_to_compact<std::int32_t>(e->bs.z)
        });
        finder.set_hits(layer1, layer2);

        formatting += std::chrono::high_resolution_clock::now() - start;
        formatted_hits += layer1.size();
        formatted_hits += layer1.size();
        start = std::chrono::high_resolution_clock::now();

        finder.start();

        std::vector<pb_doublet_finder::doublet> doublets;
        finder.get_doublets(doublets);

        finding += std::chrono::high_resolution_clock::now() - start;
        doublets_found += doublets.size();

        if (doublets.empty()) {
            std::cout << "No doublets found!" << std::endl;
            continue;
        }
        std::cout << "Doublets: "
                  << doublets.size()
                  << "; factor: "
                  << layer1.size() * layer2.size() / doublets.size() << std::endl;

        auto previous_size = doublets.size();
        std::sort(doublets.begin(), doublets.end());
        doublets.erase(std::unique(doublets.begin(), doublets.end()), doublets.end());

        if (doublets.size() != previous_size) {
            std::cout << "Erased " << (previous_size - doublets.size())
                      << " duplicates!" << std::endl;
        }

        for (const auto &doublet : doublets) {
            const auto &h1 = layer1.at(doublet.first);
            const auto &h2 = layer2.at(doublet.second);

            doublet_phi1.Fill(compact_to_radians(h1.phi));
            doublet_phi2.Fill(compact_to_radians(h2.phi));
            doublet_phi2_phi1.Fill(compact_to_radians(h2.phi - h1.phi));
            doublet_z1.Fill(compact_to_length(h1.z));
            doublet_z2.Fill(compact_to_length(h2.z));
        }
    }

    std::cout << "==== Performance info ====" << std::endl;
    std::cout << "Formatted " << formatted_hits
              << " hits in " << formatting.count()
              << " s (" << (formatting.count() / formatted_hits * 1e6)
              << " us)" << std::endl;
    std::cout << "Found " << doublets_found
              << " doublets in " << finding.count()
              << " s (" << (finding.count() / formatted_hits * 1e6)
              << " us)" << std::endl;

    out.cd();
    out.Write();
}
