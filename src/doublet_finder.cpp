#include "doublet_finder.h"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "fast_sincos.h"

std::vector<cpu_doublet_finder::hit_type> cpu_doublet_finder::convert(
        const std::vector<hit> &hits, int layer) const
{
    std::vector<hit_type> res;
    res.reserve(hits.size());
    for (const hit &h : hits) {
        res.emplace_back(h, 1);
    }
    return res;
}

std::size_t cpu_doublet_finder::get_doublets(
    std::vector<cpu_doublet_finder::doublet_type> &output)
{
    if (output.size() == 0) {
        std::swap(_doublets, output);
        return output.size();
    } else {
        std::size_t count = _doublets.size();
        std::copy(_doublets.begin(), _doublets.end(),
                  std::back_inserter(output));
        _doublets.clear();
        return count;
    }
}

void cpu_doublet_finder::sort_hits(
        std::vector<cpu_doublet_finder::hit_type> &layer1,
        std::vector<cpu_doublet_finder::hit_type> &layer2)
{
    std::sort(layer1.begin(),
              layer1.end(),
              [](const hit_type &a, const hit_type &b) {
                  return a.phi < b.phi;
              });
    std::sort(layer2.begin(),
              layer2.end(),
              [](const hit_type &a, const hit_type &b) {
                  return a.phi < b.phi;
              });
}

namespace /* anonymous */
{
    /**
     * \brief Checks that the z component of the impact parameter is within the
     *        beam spot
     *
     * The error on the computed value is about 1.4mm.
     */
    bool check_dz(const compact_hit &inner,
                  const compact_hit &outer,
                  int rb_proj,
                  int b_dz)
    {
        const constexpr int layer_1_r = length_to_compact<int>(3);
        const constexpr int layer_2_r = length_to_compact<int>(6.8);

        int inner_r = layer_1_r + inner.dr;
        int outer_r = layer_2_r + outer.dr;

        int num_xi = inner_r - rb_proj;
        int dz = outer.z - inner.z;

        int dr = outer_r - inner_r;

        num_xi >>= 8;
        dz >>= 8;
        dr >>= 8;

        int dz_times_dr = dr * b_dz - dz * num_xi;

        int bound = length_to_compact<int>(11) * std::abs(dr) >> 8;

        return std::abs(dz_times_dr) < bound;
    }
} // namespace anonymous

void cpu_doublet_finder::find(
        const cpu_doublet_finder::beam_spot_type &bs,
        const std::vector<cpu_doublet_finder::hit_type> &layer1,
        const std::vector<cpu_doublet_finder::hit_type> &layer2)
{
    if (layer1.empty() || layer2.empty()) {
        return;
    }

    std::size_t index = 0;

    _doublets.resize(layer1.size() * layer2.size() / 64);

    const std::int16_t window_width = radians_to_compact(0.04);

    fast_sincos sincos(bs.phi - layer1.front().phi);
    std::size_t iterations = 0;

    auto range_begin = layer2.begin();
    auto range_end = range_begin;

    for (auto it1 = layer1.begin(); it1 != layer1.end(); ++it1) {
        const auto &inner = *it1;

        sincos.step(bs.phi - inner.phi);
        if (iterations % 64 == 0) {
            sincos.sync(bs.phi - inner.phi);
        }
        iterations++;

        int rb_proj = sincos.cos_times(bs.r);
        int b_dz = (inner.z - bs.z) >> 8;

        // We can't use an int16 here, else it wraps around in the first
        // iteration, gets negative and the condition in the while loop is
        // always false
        const int phi_low = inner.phi - window_width;
        while (range_begin->phi < phi_low && range_begin != layer2.end()) {
            ++range_begin;
        }

        // We can't use an int16 here, else it wraps around and the break
        // below happens too early.
        const int phi_high = inner.phi + window_width;
        while (range_end->phi <= phi_high && range_end != layer2.end()) {
            ++range_end;
        }

        for (auto it2 = range_begin; it2 != range_end; ++it2) {
            if (check_dz(inner, *it2, rb_proj, b_dz)) {
                _doublets[index].first = std::distance(layer1.begin(), it1);
                _doublets[index].second = std::distance(layer2.begin(), it2);
                ++index;
            }
        }
    }

    /* Edge cases
     * ----------
     *
     * When the first hit is close to -pi or pi, the search window in the
     * second layer extends to both the beginning and the end of the list.
     * The loop above only takes care of the case where both signs are
     * equal, so we compensate below.
     */

    // Recover efficiency near -pi
    for (auto it1 = layer1.begin(); it1 != layer1.end(); ++it1) {
        const auto &inner = *it1;

        // Here we want to check for the wraparound, so we need int16
        const std::int16_t phi_low = inner.phi - window_width;
        if (phi_low < 0) {
            // Wrapped around
            break;
        }

        // Need a float to compute the cos
        float cos = std::cos(compact_to_radians(bs.phi - inner.phi));
        int rb_proj = bs.r * cos;
        int b_dz = (inner.z - bs.z) >> 8;

        for (auto it2 = layer2.rbegin(); it2 != layer2.rend(); ++it2) {
            if (it2->phi < phi_low) {
                break;
            }
            if (check_dz(inner, *it2, rb_proj, b_dz)) {
                _doublets[index].first = std::distance(layer1.begin(), it1);
                _doublets[index].second = std::distance(it2, layer2.rend()) - 1;
                ++index;
            }
        }
    }

    // Recover efficiency near +pi
    for (auto it1 = layer1.rbegin(); it1 != layer1.rend(); ++it1) {
        const auto &inner = *it1;

        // Here we want to check for the wraparound, so we need int16
        const std::int16_t phi_high = inner.phi + window_width;
        if (phi_high > 0) {
            // Wrapped around
            break;
        }

        // Need a float to compute the cos
        float cos = std::cos(compact_to_radians(bs.phi - inner.phi));
        int rb_proj = bs.r * cos;
        int b_dz = (inner.z - bs.z) >> 8;

        for (auto it2 = layer2.begin(); it2 != layer2.end(); ++it2) {
            if (it2->phi > phi_high) {
                break;
            }
            if (check_dz(inner, *it2, rb_proj, b_dz)) {
                _doublets[index].first = std::distance(it1, layer1.rend()) - 1;
                _doublets[index].second = std::distance(layer2.begin(), it2);
                ++index;
            }
        }
    }

    _doublets.resize(index);
}

////////////////////////////////////////////////////////////////////////////////

std::size_t float_doublet_finder::get_doublets(
    std::vector<float_doublet_finder::doublet_type> &output)
{
    if (output.size() == 0) {
        std::swap(_doublets, output);
        return output.size();
    } else {
        std::size_t count = _doublets.size();
        std::copy(_doublets.begin(), _doublets.end(),
                  std::back_inserter(output));
        _doublets.clear();
        return count;
    }
}

void float_doublet_finder::sort_hits(
        std::vector<float_doublet_finder::hit_type> &layer1,
        std::vector<float_doublet_finder::hit_type> &layer2)
{
    std::sort(layer1.begin(),
              layer1.end(),
              [](const hit_type &a, const hit_type &b) {
                  return a.phi < b.phi;
              });
    std::sort(layer2.begin(),
              layer2.end(),
              [](const hit_type &a, const hit_type &b) {
                  return a.phi < b.phi;
              });
}

namespace /* anonymous */
{
    /**
     * \brief Checks that the z component of the impact parameter is within the
     *        beam spot
     *
     * The error on the computed value is about 1.4mm.
     */
    bool float_check_dz(const hit &inner,
                  const hit &outer,
                  float rb_proj,
                  float b_dz)
    {
        float num_xi = inner.r - rb_proj;
        float dz = outer.z - inner.z;

        float dr = outer.r - inner.r;

        float dz_times_dr = dr * b_dz - dz * num_xi;

        float bound = 11 * std::abs(dr);

        return std::abs(dz_times_dr) < bound;
    }
} // namespace anonymous

void float_doublet_finder::find(
        const float_doublet_finder::beam_spot_type &bs,
        const std::vector<float_doublet_finder::hit_type> &layer1,
        const std::vector<float_doublet_finder::hit_type> &layer2)
{
    if (layer1.empty() || layer2.empty()) {
        return;
    }

    std::size_t index = 0;

    _doublets.resize(layer1.size() * layer2.size() / 64);

    const float window_width = 0.04f;

    fast_float_sincos sincos(bs.phi - layer1.front().phi);
    std::size_t iterations = 0;

    auto range_begin = layer2.begin();
    auto range_end = range_begin;

    for (auto it1 = layer1.begin(); it1 != layer1.end(); ++it1) {
        const auto &inner = *it1;

        sincos.step(bs.phi - inner.phi);
        if (iterations % 64 == 0) {
            sincos.sync(bs.phi - inner.phi);
        }
        iterations++;

        float rb_proj = sincos.cos_times(bs.r);
        float b_dz = inner.z - bs.z;

        float phi_low = inner.phi - window_width;
        while (range_begin->phi < phi_low && range_begin != layer2.end()) {
            ++range_begin;
        }

        // We can't use an int16 here, else it wraps around and the break
        // below happens too early.
        float phi_high = inner.phi + window_width;
        while (range_end->phi <= phi_high && range_end != layer2.end()) {
            ++range_end;
        }

        for (auto it2 = range_begin; it2 != range_end; ++it2) {
            if (float_check_dz(inner, *it2, rb_proj, b_dz)) {
                _doublets[index].first = std::distance(layer1.begin(), it1);
                _doublets[index].second = std::distance(layer2.begin(), it2);
                ++index;
            }
        }
    }

    /* Edge cases
     * ----------
     *
     * When the first hit is close to -pi or pi, the search window in the
     * second layer extends to both the beginning and the end of the list.
     * The loop above only takes care of the case where both signs are
     * equal, so we compensate below.
     */

    // Recover efficiency near -pi
    for (auto it1 = layer1.begin(); it1 != layer1.end(); ++it1) {
        const auto &inner = *it1;

        const float phi_low = inner.phi - window_width;
        if (phi_low < 0) {
            // Wrapped around
            break;
        }

        // Need a float to compute the cos
        float rb_proj = bs.r * std::cos(bs.phi - inner.phi);
        float b_dz = inner.z - bs.z;

        for (auto it2 = layer2.rbegin(); it2 != layer2.rend(); ++it2) {
            if (it2->phi < phi_low) {
                break;
            }
            if (float_check_dz(inner, *it2, rb_proj, b_dz)) {
                _doublets[index].first = std::distance(layer1.begin(), it1);
                _doublets[index].second = std::distance(it2, layer2.rend()) - 1;
                ++index;
            }
        }
    }

    // Recover efficiency near +pi
    for (auto it1 = layer1.rbegin(); it1 != layer1.rend(); ++it1) {
        const auto &inner = *it1;

        const float phi_high = inner.phi + window_width;
        if (phi_high > 0) {
            // Wrapped around
            break;
        }

        // Need a float to compute the cos
        float rb_proj = bs.r * std::cos(bs.phi - inner.phi);
        float b_dz = inner.z - bs.z;

        for (auto it2 = layer2.begin(); it2 != layer2.end(); ++it2) {
            if (it2->phi > phi_high) {
                break;
            }
            if (float_check_dz(inner, *it2, rb_proj, b_dz)) {
                _doublets[index].first = std::distance(it1, layer1.rend()) - 1;
                _doublets[index].second = std::distance(layer2.begin(), it2);
                ++index;
            }
        }
    }

    _doublets.resize(index);
}
