#include <array>
#include <chrono>
#include <iostream>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TPie.h>
#include <TTree.h>

#include "doublet_finder.h"
#include "eventreader.h"
#include "geometry.h"
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

double extrapolated_xi(const compact_beam_spot &bs,
                       const compact_hit &inner,
                       const compact_hit &outer)
{
    const constexpr int layer_1_r = length_to_compact<int>(3);
    const constexpr int layer_2_r = length_to_compact<int>(6.8);

    int inner_r = layer_1_r + inner.dr;
    int outer_r = layer_2_r + outer.dr;

    int dr = outer_r - inner_r;

    // Need a float to compute the cos
    float inner_phi = compact_to_radians(inner.phi);

    int rb_proj = length_to_compact<int>(bs.r * std::cos(bs.phi - inner_phi));

    int num = inner_r - rb_proj;
    return -double(num) / dr;
}

/**
 * \brief Finds transverse component of the impact parameter (wrt the beam spot)
 */
int extrapolated_dr(const compact_beam_spot &bs,
                    const compact_hit &inner,
                    const compact_hit &outer)
{
    const constexpr int layer_1_r = length_to_compact<int>(3);
    const constexpr int layer_2_r = length_to_compact<int>(6.8);

    int inner_r = layer_1_r + inner.dr;

    // Need a float to compute the cos and sin
    float inner_phi = compact_to_radians(inner.phi);

    int rb_proj_x = length_to_compact<int>(bs.r * std::cos(bs.phi - inner_phi));
    int rb_proj_y = length_to_compact<int>(bs.r * std::sin(bs.phi - inner_phi));

    std::int16_t dphi = outer.phi - inner.phi;

    return std::abs(((inner_r - rb_proj_x) * dphi) / radians_to_compact(1) + rb_proj_y);
}

/**
 * \brief Finds z component of the impact parameter (wrt the beam spot)
 *
 * The error on the computed value is about 1.4mm.
 */
int extrapolated_dz(const compact_beam_spot &bs,
                    const compact_hit &inner,
                    const compact_hit &outer)
{
    double xi = extrapolated_xi(bs, inner, outer);
    return inner.z + (outer.z - inner.z) * xi - bs.z;
}

int main(int, char **)
{
    std::chrono::duration<double> formatting, formatting_acc;
    long long formatted_hits = 0;

    std::chrono::duration<double> sorting, sorting_acc;
    long long sorted_hits = 0;

    std::chrono::duration<double> finding, finding_acc;
    long long doublets_found = 0;
   
    bool do_validation = 0;

    TFile out("doublets.root", "RECREATE");
    TTree tree("doublets", "doublets");

    std::vector<int> doublets_inner;
    tree.Branch("inner", &doublets_inner);

    std::vector<int> doublets_outer;
    tree.Branch("outer", &doublets_outer);

    double formatting_seconds, sorting_seconds, finding_seconds, total_seconds;
    tree.Branch("formatting_seconds", &formatting_seconds);
    tree.Branch("sorting_seconds", &sorting_seconds);
    tree.Branch("finding_seconds", &finding_seconds);
    tree.Branch("total_seconds", &total_seconds);

    TH1D doublet_phi1("doublet_phi1", ";phi1;count", 50, -pi, pi);
    TH1D doublet_phi2("doublet_phi2", ";phi2;count", 50, -pi, pi);
    TH1D doublet_phi2_phi1("doublet_phi2_phi1", ";phi2 - phi1;count", 50, -0.05, 0.05);
    TH1D doublet_z1("doublet_z1", ";z1;count", 50, -30, 30);
    TH1D doublet_z2("doublet_z2", ";z2;count", 50, -30, 30);
    TH1D doublet_z0("doublet_z0", ";z0;count", 50, -15, 15);
    TH1D doublet_b0("doublet_b0", ";b0;count", 50, 0, 0.25);

    TH1D hit_count_1("hit_count_1", ";Hits in layer 1;Events", 50, 0, 1500);
    TH1D hit_count_2("hit_count_2", ";Hits in layer 2;Events", 50, 0, 1500);
    TH1D hit_count_12("hit_count_12", ";Number of naive doublets;Events", 50, 0, 2e6);
    TH1D doublet_count("doublet_count", ";Number of doublets;Events", 50, 0, 7000);

    TH1D duration("duration", ";Duration (us);Events", 50, 0, 2500);
    TH2D duration_vs_nvtx("duration_vs_nvtx", ";#vtx;Duration (us)", 50, 0, 75, 50, 0, 2500);

    TH1D trk_pt("trk_pt", "trk_pt", 20, 0, 3.0); trk_pt.Sumw2();
    TH1D trk_eta("trk_eta", "trk_eta", 20, -2.5, 2.5); trk_eta.Sumw2();
    TH1D trk_phi("trk_phi", "trk_phi", 20, -3.14, 3.14); trk_phi.Sumw2();
   
    const int doublet_eff_pt_nbins = 7;
    float doublet_eff_pt_xbins[doublet_eff_pt_nbins+1] = {0.7,0.8,0.9,1.0,1.25,1.5,2.0,3.0};
    TH1D doublet_all_pt("doublet_all_pt", "doublet_all_pt", doublet_eff_pt_nbins, doublet_eff_pt_xbins);
    TH1D doublet_pass_pt("doublet_pass_pt", "doublet_pass_pt", doublet_eff_pt_nbins, doublet_eff_pt_xbins);
    doublet_all_pt.Sumw2();
    doublet_pass_pt.Sumw2();

    const int doublet_eff_eta_nbins = 10;
    float doublet_eff_eta_xbins[doublet_eff_eta_nbins+1] = {-2.5,-2.0,-1.5,-1.0,-0.5,0.,0.5,1.0,1.5,2.0,2.5};
    TH1D doublet_all_eta("doublet_all_eta", "doublet_all_eta", doublet_eff_eta_nbins, doublet_eff_eta_xbins);
    TH1D doublet_pass_eta("doublet_pass_eta", "doublet_pass_eta", doublet_eff_eta_nbins, doublet_eff_eta_xbins);
    doublet_all_eta.Sumw2();
    doublet_pass_eta.Sumw2();

    const int doublet_eff_phi_nbins = 12;
    float doublet_eff_phi_xbins[doublet_eff_phi_nbins+1] = {-3.14,-2.5,-2.0,-1.5,-1.0,-0.5,0.,0.5,1.0,1.5,2.0,2.5,3.14};
    TH1D doublet_all_phi("doublet_all_phi", "doublet_all_phi", doublet_eff_phi_nbins, doublet_eff_phi_xbins);
    TH1D doublet_pass_phi("doublet_pass_phi", "doublet_pass_phi", doublet_eff_phi_nbins, doublet_eff_phi_xbins);
    doublet_all_phi.Sumw2();
    doublet_pass_phi.Sumw2();
   
//     event_reader in("output.root");

    event_reader in("~lmoureau/data/v3.root");

    long long i = 0;
    float n_doub_to_track = 0;
    float n_track = 0;
    while (in.next()) {
        i++;
        std::cout << "==== Next event ====" << std::endl;

        std::unique_ptr<event> e = in.get();

        std::vector<track> interesting_tracks;
        std::copy_if(e->tracks.begin(),
                     e->tracks.end(),
                     std::back_inserter(interesting_tracks),
                     [](const track &trk) {
                        if (trk.pt < 0.7) return false;

                        bool has_hit_in_layer_1 = false;
                        bool has_hit_in_layer_2 = false;

                        for (const hit &hh : trk.hits) {
                            if (hit_is_pixel_barrel(hh)) {
                                has_hit_in_layer_1 |= (hit_pixel_barrel_layer(hh) == 0);
                                has_hit_in_layer_2 |= (hit_pixel_barrel_layer(hh) == 1);

                                if (has_hit_in_layer_1 && has_hit_in_layer_2) {
                                    break;
                                }
                            }
                        }
                        return has_hit_in_layer_1 && has_hit_in_layer_2;
                     }
                    );

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

        std::vector<compact_hit> layer1;
        layer1.reserve(pb_hits_per_layer[0].size());
        for (const hit &h : pb_hits_per_layer[0]) {
            layer1.emplace_back(h, 0);
        }

        std::vector<compact_hit> layer2;
        layer2.reserve(pb_hits_per_layer[1].size());
        for (const hit &h : pb_hits_per_layer[1]) {
            layer2.emplace_back(h, 1);
        }

        cpu_doublet_finder finder;

        formatting = std::chrono::high_resolution_clock::now() - start;
        formatting_acc += formatting;
        sorted_hits += layer1.size();
        sorted_hits += layer2.size();
        auto sorting_start = std::chrono::high_resolution_clock::now();

        finder.sort_hits(layer1, layer2);

        compact_beam_spot bs{
            length_to_compact<std::int32_t>(e->bs.r),
            length_to_compact<std::int32_t>(e->bs.z),
            radians_to_compact(e->bs.phi)
        };

        finder.set_beam_spot(bs);
        finder.set_hits(layer1, layer2);

        sorting = std::chrono::high_resolution_clock::now() - sorting_start;
        sorting_acc += sorting;
        formatted_hits += layer1.size();
        formatted_hits += layer2.size();
        auto finding_start = std::chrono::high_resolution_clock::now();

        finder.start();

        std::vector<cpu_doublet_finder::doublet> doublets;
        finder.get_doublets(doublets);

        auto end = std::chrono::high_resolution_clock::now();
        finding = end - finding_start;
        finding_acc += finding;

        std::chrono::duration<double> event_duration = end - start;
        duration_vs_nvtx.Fill(e->nvtx, 1e6 * event_duration.count());
        duration.Fill(1e6 * event_duration.count());

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

        doublets_inner.clear();
        doublets_outer.clear();

        formatting_seconds = formatting.count();
        sorting_seconds = sorting.count();
        finding_seconds = finding.count();
        total_seconds = event_duration.count();

       int ndoub = 0;
        for (const auto &doublet : doublets) {
            doublets_inner.push_back(doublet.first);
            doublets_outer.push_back(doublet.second);

            const auto &h1 = layer1.at(doublet.first);
            const auto &h2 = layer2.at(doublet.second);

            doublet_phi1.Fill(compact_to_radians(h1.phi));
            doublet_phi2.Fill(compact_to_radians(h2.phi));
            doublet_phi2_phi1.Fill(compact_to_radians(h2.phi - h1.phi));
            doublet_z1.Fill(compact_to_length(h1.z));
            doublet_z2.Fill(compact_to_length(h2.z));

            doublet_z0.Fill(compact_to_length(extrapolated_dz(bs, h1, h2)));
            doublet_b0.Fill(compact_to_length(extrapolated_dr(bs, h1, h2)));
	    
	    if( do_validation ) {
		
                for (const track &t : interesting_tracks) {
		  bool foundh1 = 0;
		  bool foundh2 = 0;
		  
		  for (const hit &hh : t.hits) {
		     
		     if( ! hit_is_pixel_barrel(hh) ) continue;

                     int layer = hit_pixel_barrel_layer(hh);
                     if (layer == 0) {
                        bool pass_phi = (radians_to_compact(hh.phi) == h1.phi);
                        bool pass_z = (length_to_compact<std::int32_t>(hh.z) == h1.z);
                        bool pass_dr = (length_to_compact<std::int16_t>(hh.r - 3) == h1.dr);

                        foundh1 |= (pass_phi && pass_z && pass_dr);
                     } else if (layer == 1) {
                        bool pass_phi = (radians_to_compact(hh.phi) == h2.phi);
                        bool pass_z = (length_to_compact<std::int32_t>(hh.z) == h2.z);
                        bool pass_dr = (length_to_compact<std::int16_t>(hh.r - 6.8) == h2.dr);

                        foundh2 |= (pass_phi && pass_z && pass_dr);
                     }

                     if (foundh1 && foundh2) {
                         break;
                     }
		  }
		  
		  if( foundh1 && foundh2 ) {
		     
		     doublet_pass_pt.Fill(t.pt);
		     doublet_pass_eta.Fill(t.eta);
		     doublet_pass_phi.Fill(t.phi);
		     
		     n_doub_to_track++;
		     break;
		     
		  }
		  
	       }	   
	    }	   
	}

        tree.Fill();
       
       if( do_validation ) {
	  
            for (const track &t : interesting_tracks) {
	     doublet_all_pt.Fill(t.pt);
	     doublet_all_eta.Fill(t.eta);
	     doublet_all_phi.Fill(t.phi);
	     
	     trk_pt.Fill(t.pt);
	     trk_eta.Fill(t.eta);
	     trk_phi.Fill(t.phi);
	     
	     n_track++;
	  }
        }

        hit_count_1.Fill(layer1.size());
        hit_count_2.Fill(layer2.size());
        hit_count_12.Fill(layer1.size() * layer2.size());
        doublet_count.Fill(doublets.size());
    }

    if( do_validation ) {
       
       const int nvals = 2;
       float vals[nvals] = {n_doub_to_track,n_track-n_doub_to_track};
       int colors[nvals] = {3,2};
       TPie *pie = new TPie("pie","pie",nvals,vals,colors);
       pie->SetRadius(.2);
       pie->SetLabelsOffset(.01);
       pie->SetLabelFormat("#splitline{%perc}{%txt}");
       pie->SetEntryLabel(0,"Found");
       pie->SetEntryLabel(1,"Missed");
       
       doublet_pass_pt.Divide(&doublet_pass_pt,&doublet_all_pt,1,1,"B");
       doublet_pass_eta.Divide(&doublet_pass_eta,&doublet_all_eta,1,1,"B");
       doublet_pass_phi.Divide(&doublet_pass_phi,&doublet_all_phi,1,1,"B");

       out.cd();
       if( do_validation ) pie->Write();
    }   
   
    std::cout << "==== Performance info ====" << std::endl;
    std::cout << "Formatted " << formatted_hits
              << " hits in " << formatting_acc.count()
              << " s (" << (1e6 * formatting_acc.count() / i)
              << " us/event)" << std::endl;
    std::cout << "Sorted    " << sorted_hits
              << " hits in " << sorting_acc.count()
              << " s (" << (1e6 * sorting_acc.count() / i)
              << " us/event)" << std::endl;
    std::cout << "Found " << doublets_found
              << " doublets in " << finding_acc.count()
              << " s (" << (1e6 * finding_acc.count() / i)
              << " us/event)" << std::endl;
   
    if( do_validation )
     std::cout << n_doub_to_track << " of doublets are found in " << n_track << " tracks " << std::endl;

    out.cd();
    out.Write();
}
