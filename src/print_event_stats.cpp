#include <array>
#include <iostream>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include "compact.h"
#include "eventreader.h"
#include "hitutils.h"

float dphi3(float phi1, float phi2, float phi3)
{
    float dphi3 = phi3 + phi1 - 2 * phi2;
    if (phi3 - phi2 < 0) {
        dphi3 = -dphi3;
    }
    return dphi3;
}

float deltaphi(float phi1, float phi2)
{
    float delta = phi1 - phi2;
    if (delta > pi) {
        delta -= 2 * pi;
    }
    return delta;
}

int main(int, char **)
{
    TFile out("output.root", "RECREATE");

    TH1D dphi("dphi3", "dphi3", 100, -0.25, 0.25);
    TH1D dz3("dz3", "dz3", 100, -10, 10);
    TH1D dz4("dz4", "dz4", 100, -10, 10);
    TH2D dphi_dphi3("dphi_dphi3", "dphi vs dphi3", 50, 0, 0.5, 60, -0.1, 0.5);
    TH2D dphi3_dz3("dphi3_dz3", "dphi3 vs dz3", 50, -0.1, 0.4, 60, 0, 30);

    TH2D track_phi1_phi2("track_phi1_phi2", ";phi1 - bs phi;phi2 - phi1", 50, -pi, pi, 50, -0.5, 0.5);
    TH2D track_proj_z_dphi("track_proj_z_diphi", ";Projected z;phi2 - phi1", 50, -20, 20, 50, -0.05, 0.05);
    TH2D track_z0_2hits("track_z0_2hits", ";z0;z0(2hits)", 50, -20, 20, 50, -20, 20);
    TH1D track_z0_resolution("track_z0_resolution", ";#Delta z_{0} (2 hits)", 50, -0.2, -0.1);
    TH1D track_z0_resolution_int("track_z0_resolution_int", ";#Delta z_{0} (2 hits)", 50, -0.2, -0.1);

    TH2D track_phi2_phi3("track_phi2_phi3", ";phi2;phi3", 50, -0.5, 0.5, 50, -0.5, 0.5);
    TH2D track_phi2_phi4("track_phi2_phi4", ";phi2;phi4", 50, -0.5, 0.5, 50, -0.5, 0.5);
    TH2D track_phi3_phi4("track_phi3_phi4", ";phi3;phi4", 50, -0.5, 0.5, 50, -0.5, 0.5);

    TH2D seed_phi2_phi3("seed_phi2_phi3", ";phi2;phi3", 50, -0.5, 0.5, 50, -0.5, 0.5);
    TH2D seed_phi2_phi4("seed_phi2_phi4", ";phi2;phi4", 50, -0.5, 0.5, 50, -0.5, 0.5);
    TH2D seed_phi3_phi4("seed_phi3_phi4", ";phi3;phi4", 50, -0.5, 0.5, 50, -0.5, 0.5);

    event_reader in("~lmoureau/data/v3.root");
    long long i = 0;
    while (in.next()) {
        std::cout << "===== Event " << i++ << " =====" << std::endl;
        std::unique_ptr<event> e = in.get();
        std::cout << "Beam spot r phi z: "
                  << e->bs.r << " "
                  << e->bs.phi << " "
                  << e->bs.z << std::endl;
        std::cout << "#hits:       " << e->hits.size() << std::endl;
        std::sort(e->hits.begin(), e->hits.end(), hit_less_than);
        e->hits.erase(std::unique(e->hits.begin(), e->hits.end(), hit_equal),
                      e->hits.end());
        std::cout << "  unique:    " << e->hits.size() << std::endl;
        std::cout << "#tracks:     " << e->tracks.size() << std::endl;

        std::array<std::vector<hit>, 4> pb_seeds_per_layer;
        for (const track &trk : e->tracks) {
            for (const hit &h : trk.seed) {
                if (hit_is_pixel_barrel(h)) {
                    pb_seeds_per_layer[hit_pixel_barrel_layer(h)].push_back(h);
                }
            }
            float phi0 = 0, phi1 = 0, phi2 = 0, phi3 = 0;
            float z0 = 0, z1 = 0, z2 = 0, z3 = 0;
            float r0 = 0, r1 = 0, r2 = 0, r3 = 0;
            for (const hit &h : trk.hits) {
                if (hit_is_pixel_barrel(h)) {
                    switch (hit_pixel_barrel_layer(h)) {
                    case 0:
                        phi0 = h.phi;
                        z0 = h.z;
                        r0 = h.r;
                        break;
                    case 1:
                        phi1 = h.phi;
                        z1 = h.z;
                        r1 = h.r;
                        break;
                    case 2:
                        phi2 = h.phi;
                        z2 = h.z;
                        r2 = h.r;
                        break;
                    case 3:
                        phi3 = h.phi;
                        z3 = h.z;
                        r3 = h.r;
                        break;
                    }
                }
            }
            if (phi0 != 0 && phi1 != 0 && phi2 != 0) {
                if (phi2 - phi1 > 0) {
                    dphi.Fill(phi2 + phi0 - 2 * phi1);
                    dphi_dphi3.Fill(phi2 - phi1, phi2 + phi0 - 2 * phi1);
                } else {
                    dphi.Fill(-(phi2 + phi0 - 2 * phi1));
                    dphi_dphi3.Fill(phi1 - phi2, -(phi2 + phi0 - 2 * phi1));
                }
                if (z0 != 0 && z1 != 0 && z2 != 0) {
                    dphi3_dz3.Fill(dphi3(phi0, phi1, phi2), std::abs(z2 + z0 - 2 * z1));
                }
                if (trk.pt > 0.8) {
                    track_phi1_phi2.Fill(deltaphi(phi1, e->bs.phi), deltaphi(phi2, phi1));
                    track_proj_z_dphi.Fill(r2 - (r2 - r1) / (z2 - z1) * z2, deltaphi(phi2, phi1));

                    {
                        // This calculation is EXACT (for pt -> inf)
                        float x1 = r0;
                        float x2 = r1 * std::cos(phi1 - phi0);
                        float y1 = 0;
                        float y2 = r1 * std::sin(phi1 - phi0);

                        float xb = e->bs.r * std::cos(e->bs.phi - phi0);
                        float yb = e->bs.r * std::sin(e->bs.phi - phi0);

                        float num = (x2 - x1) * (x1 - xb) + (y2 - y1) * (y1 - yb);
                        float den = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
                        float xi = -num / den;

                        track_z0_2hits.Fill(trk.z0, z0 + (z1 - z0) * xi);
                    }
                    {
                        // This calculation is APPROXIMATE but very precise too (and much simpler)
                        float rb_proj = e->bs.r * std::cos(e->bs.phi - phi0);

                        float num = (rb_proj - r0);
                        float den = (r1 - r0);
                        float xi = num / den;

                        track_z0_resolution.Fill(trk.z0 - (z0 + (z1 - z0) * xi));
                    }
                    {
                        // This calculation is the approx one but without floats (except for the cos)
                        int ir0 = length_to_compact<int>(r0);
                        int ir1 = length_to_compact<int>(r1);
                        int dr = ir1 - ir0;

                        int rb_proj = length_to_compact<int>(e->bs.r * std::cos(e->bs.phi - phi0));
                        int num = ir0 - rb_proj;

                        int xi = -(num << 12) / dr;

                        int iz0 = length_to_compact<int>(z0);
                        int iz1 = length_to_compact<int>(z1);

                        int z0_est = iz0 + (((iz1 - iz0) * xi) >> 12);

                        track_z0_resolution_int.Fill(trk.z0 - compact_to_length(z0_est));
                    }
                }
                track_phi2_phi3.Fill(phi2 - phi3, phi1 - phi3);
                track_phi2_phi4.Fill(phi2 - phi3, phi0 - phi3);
                track_phi3_phi4.Fill(phi1 - phi3, phi0 - phi3);
            }
            if (z0 != 0 && z1 != 0 && z2 != 0) {
                dz3.Fill(z2 + z0 - 2 * z1);
                dz4.Fill(z0 - z1 - z2 + z3);
            }
        }
        for (unsigned layer = 0; layer < pb_seeds_per_layer.size(); ++layer) {
            std::cout << "  layer " << layer << ":   " << pb_seeds_per_layer[layer].size() << std::endl;
        }

        std::vector<hit> pb_hits;
        std::copy_if(e->hits.begin(), e->hits.end(), std::back_inserter(pb_hits), hit_is_pixel_barrel);
        std::cout << "#pb hits:    " << pb_hits.size() << std::endl;

        std::array<std::vector<hit>, 4> pb_hits_per_layer;
        for (const hit &h : pb_hits) {
            pb_hits_per_layer[hit_pixel_barrel_layer(h)].push_back(h);
        }
        for (unsigned layer = 0; layer < pb_hits_per_layer.size(); ++layer) {
            std::cout << "  layer " << layer << ":   " << pb_hits_per_layer[layer].size() << std::endl;
        }

        continue; // Don't run seed finding

        std::cout << "Searching for hit pairs in layers 3 and 4..." << std::flush;
        std::vector<std::vector<hit>> candidates;
        for (const hit &h3 : pb_hits_per_layer[3]) {
            for (const hit &h2 : pb_hits_per_layer[2]) {
                if (std::abs(h3.phi - h2.phi) < 0.2) {
                    candidates.push_back({ h3, h2 });
                }
            }
        }
        std::cout << " " << candidates.size() << " candidates." << std::endl;

        std::cout << "Matching to layer 2..." << std::flush;
        std::vector<std::vector<hit>> new_candidates;
        for (const auto &pair : candidates) {
            const hit &h3 = pair[0];
            const hit &h2 = pair[1];

            for (const hit &h1 : pb_hits_per_layer[1]) {
                float dphi3 = ::dphi3(h1.phi, h3.phi, h3.phi);
                if (-0.01 < dphi3 && dphi3 < 0.05
                    && std::abs(h3.z + h1.z - 2 * h2.z) < 2.5) {
                    std::vector<hit> newc = pair;
                    newc.push_back(h1);
                    new_candidates.push_back(newc);
                }
            }
        }
        std::swap(candidates, new_candidates);
        new_candidates.clear();

        std::cout << " " << candidates.size() << " candidates." << std::endl;

        std::cout << "Matching to layer 1..." << std::flush;
        for (const auto &pair : candidates) {
            const hit &h3 = pair[0];
            const hit &h2 = pair[1];
            const hit &h1 = pair[2];

            for (const hit &h0 : pb_hits_per_layer[0]) {
                float dphi3 = ::dphi3(h0.phi, h1.phi, h2.phi);
                if (-0.03 < dphi3 && dphi3 < 0.05
                    && std::abs(h3.z + h0.z - 3 * h2.z) < 1.5) {
                    std::vector<hit> newc = pair;
                    newc.push_back(h0);
                    new_candidates.push_back(newc);
                }
            }
        }
        std::swap(candidates, new_candidates);
        new_candidates.clear();

        std::cout << " " << candidates.size() << " candidates." << std::endl;
        std::cout << "Plotting..." << std::endl;

        for (const auto &seed : candidates) {
            seed_phi2_phi3.Fill(seed[2].phi - seed[3].phi, seed[1].phi - seed[3].phi);
            seed_phi2_phi4.Fill(seed[2].phi - seed[3].phi, seed[0].phi - seed[3].phi);
            seed_phi3_phi4.Fill(seed[1].phi - seed[3].phi, seed[0].phi - seed[3].phi);
        }
    }

    out.cd();
    out.Write();
}
