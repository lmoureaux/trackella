#include "doublet_finder.h"

#include <algorithm>
#include <cassert>
#include <cmath>

std::size_t pb_doublet_finder::get_doublets(
    std::vector<pb_doublet_finder::doublet> &output)
{
    assert(get_state() & state_out_of_memory || get_state() & state_finished);

#ifdef HARDWARE_ACCELERATOR
    // TODO
    send_command(command::done_reading);
    return 0;
#else // HARDWARE_ACCELERATOR
    _state |= state_ready;
    _layer1 = nullptr;
    _layer2 = nullptr;

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
#endif // HARDWARE_ACCELERATOR
}

pb_doublet_finder::state pb_doublet_finder::get_state() const
{
#ifdef HARDWARE_ACCELERATOR
    // TODO
#else // HARDWARE_ACCELERATOR
    return _state;
#endif // HARDWARE_ACCELERATOR
}

void pb_doublet_finder::set_hits(const std::vector<compact_pb_hit> &layer1,
                                 const std::vector<compact_pb_hit> &layer2)
{
    assert(get_state() & state_ready);
    assert(layer1.size() <= max_hit_count);
    assert(layer2.size() <= max_hit_count);

#ifdef HARDWARE_ACCELERATOR
    // TODO
#else // HARDWARE_ACCELERATOR
    _layer1 = &layer1;
    _layer2 = &layer2;
#endif // HARDWARE_ACCELERATOR
}

void pb_doublet_finder::set_beam_spot(const compact_beam_spot &bs)
{
    assert(get_state() == state_ready);

#ifdef HARDWARE_ACCELERATOR
    // TODO
#else // HARDWARE_ACCELERATOR
    _bs = bs;
#endif // HARDWARE_ACCELERATOR
}

void pb_doublet_finder::start()
{
    send_command(command::start);
}

#ifndef HARDWARE_ACCELERATOR
namespace /* anonymous */
{
    /**
     * \brief Checks that the z component of the impact parameter is within the
     *        beam spot
     *
     * The error on the computed value is about 1.4mm.
     */
    bool check_dz(const compact_beam_spot &bs,
                  const compact_pb_hit &inner,
                  const compact_pb_hit &outer,
                  int rb_proj)
    {
        const constexpr int layer_1_r = length_to_compact<int>(3);
        const constexpr int layer_2_r = length_to_compact<int>(6.8);

        int inner_r = layer_1_r + inner.dr;
        int outer_r = layer_2_r + outer.dr;

        int num_xi = inner_r - rb_proj;
        int dz = outer.z - inner.z;

        int dr = outer_r - inner_r;
        int a = inner.z - bs.z;

        num_xi >>= 8;
        dz >>= 8;
        dr >>= 8;
        a >>= 8;

        int dz_times_dr = dr * a - dz * num_xi;

        int bound = length_to_compact<int>(11) * std::abs(dr) >> 8;

        return std::abs(dz_times_dr) < bound;
    }
} // namespace anonymous
#endif // HARDWARE_ACCELERATOR

void pb_doublet_finder::send_command(pb_doublet_finder::command cmd)
{
    if (cmd == command::start) {
        assert(get_state() & state_ready);
    } else if (cmd == command::done_reading) {
        assert(get_state() & state_out_of_memory || get_state() & state_finished);
    }

#ifdef HARDWARE_ACCELERATOR
    // TODO
#else // HARDWARE_ACCELERATOR
    if (cmd == command::start) {
        _state ^= state_ready;
        _state ^= state_out_of_memory;
        _state |= state_processing;

        _doublets.reserve(_layer1->size() * _layer2->size() / 128);

        const std::int16_t window_width = radians_to_compact(0.04);

        auto range_begin = _layer2->begin();

        for (auto it1 = _layer1->begin(); it1 != _layer1->end(); ++it1) {
            const auto &inner = *it1;

            // We can't use an int16 here, else it wraps around in the first
            // iteration, gets negative and the condition in the while loop is
            // always false
            const int phi_low = inner.phi - window_width;
            while (range_begin != _layer2->end() && range_begin->phi < phi_low) {
                ++range_begin;
            }
            // Need a float to compute the cos
            float inner_phi = compact_to_radians(inner.phi);
            int rb_proj = length_to_compact<int>(_bs.r * std::cos(_bs.phi - inner_phi));

            // We can't use an int16 here, else it wraps around and the break
            // below happens too early.
            const int phi_high = inner.phi + window_width;
            for (auto it2 = range_begin; it2 != _layer2->end(); ++it2) {
                if (it2->phi > phi_high) {
                    break;
                }
                if (check_dz(_bs, inner, *it2, rb_proj)) {
                    _doublets.push_back({
                        std::distance(_layer1->begin(), it1),
                        std::distance(_layer2->begin(), it2)
                    });
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
        for (auto it1 = _layer1->begin(); it1 != _layer1->end(); ++it1) {
            const auto &inner = *it1;

            // Here we want to check for the wraparound, so we need int16
            const std::int16_t phi_low = inner.phi - window_width;
            if (phi_low < 0) {
                // Wrapped around
                break;
            }

            // Need a float to compute the cos
            float inner_phi = compact_to_radians(inner.phi);
            int rb_proj = length_to_compact<int>(_bs.r * std::cos(_bs.phi - inner_phi));

            for (auto it2 = _layer2->rbegin(); it2 != _layer2->rend(); ++it2) {
                if (it2->phi < phi_low) {
                    break;
                }
                if (check_dz(_bs, inner, *it2, rb_proj)) {
                    _doublets.push_back({
                        std::distance(_layer1->begin(), it1),
                        std::distance(it2, _layer2->rend()) - 1
                    });
                }
            }
        }

        // Recover efficiency near +pi
        for (auto it1 = _layer1->rbegin(); it1 != _layer1->rend(); ++it1) {
            const auto &inner = *it1;

            // Here we want to check for the wraparound, so we need int16
            const std::int16_t phi_high = inner.phi + window_width;
            if (phi_high > 0) {
                // Wrapped around
                break;
            }

            // Need a float to compute the cos
            float inner_phi = compact_to_radians(inner.phi);
            int rb_proj = length_to_compact<int>(_bs.r * std::cos(_bs.phi - inner_phi));

            for (auto it2 = _layer2->begin(); it2 != _layer2->end(); ++it2) {
                if (it2->phi > phi_high) {
                    break;
                }
                if (check_dz(_bs, inner, *it2, rb_proj)) {
                    _doublets.push_back({
                        std::distance(it1, _layer1->rend()) - 1,
                        std::distance(_layer2->begin(), it2)
                    });
                }
            }
        }

        _state ^= state_processing;
        _state |= state_finished;
    }
#endif // HARDWARE_ACCELERATOR
}
