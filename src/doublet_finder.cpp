#include "doublet_finder.h"

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
    _layer1 = nullptr;
    _layer2 = nullptr;
    _state |= state_ready;

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
     * \brief Finds z component of the impact parameter (wrt the beam spot)
     *
     * The error on the computed value is about 1.4mm.
     */
    int extrapolated_dz(const compact_beam_spot &bs,
                        const compact_pb_hit &inner,
                        const compact_pb_hit &outer)
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
        int xi = -(num << 12) / dr;

        return inner.z + (((outer.z - inner.z) * xi) >> 12) - bs.z;
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

        _doublets.reserve(_layer1->size() * _layer2->size());
        for (std::uint16_t i1 = 0; i1 < _layer1->size(); ++i1) {
            for (std::uint16_t i2 = 0; i2 < _layer2->size(); ++i2) {
                const auto &inner = (*_layer1)[i1];
                const auto &outer = (*_layer2)[i2];
                if (std::abs(inner.phi - outer.phi) < radians_to_compact(0.04)
                    && std::abs(extrapolated_dz(_bs, inner, outer)) < length_to_compact<int>(11)) {
                    _doublets.push_back({ i1, i2 });
                }
            }
        }

        _state ^= state_processing;
        _state |= state_finished;
    }
#endif // HARDWARE_ACCELERATOR
}
