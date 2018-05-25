#include <cmath>
#include <atomic>
#include <cstdint>
#include <cstring>

#include <e-lib.h>

#include "epiphany_protocol.h"
#include "fast_sincos.h"

bool data_is_there()
{
    return *device_addr::status != 0;
}

void reset_flags()
{
    *device_addr::status = 0;
}

namespace /* anonymous */
{
    /**
     * \brief Checks that the z component of the impact parameter is within the
     *        beam spot
     *
     * The error on the computed value is about 1.4mm.
     */
    inline bool check_dz(const compact_hit inner,
                         const compact_hit outer,
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

        int bound = length_to_compact<int>(11) * dr;
        bound >>= 8;

        return std::abs(dz_times_dr) < bound;
    }
} // namespace anonymous

int process(const metadata &m)
{
    compact_beam_spot bs = m.bs;

    int layer1_count = m.layer1_count;
    int layer2_count = m.layer2_count;

    const compact_hit *layer1 = device_addr::data;
    const compact_hit *layer2 = device_addr::data + layer1_count;
    doublet *doublets = (doublet *)(0x8e000000 + shm_addr::doublets);

    std::size_t index = 0;

    const std::int16_t window_width = radians_to_compact(0.04);

    fast_sincos sincos(bs.phi - layer1[0].phi);
    std::size_t iterations = 0;

    auto range_begin = 0;
    auto range_end = range_begin;

    for (int it1 = 0; it1 < layer1_count; ++it1) {
        const auto &inner = layer1[it1];

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
        while (layer2[range_begin].phi < phi_low && range_begin < layer2_count) {
            ++range_begin;
        }

        // We can't use an int16 here, else it wraps around and the break
        // below happens too early.
        const int phi_high = inner.phi + window_width;
        while (layer2[range_end].phi <= phi_high && range_end < layer2_count) {
            ++range_end;
        }

        for (int it2 = range_begin; it2 < range_end; ++it2) {
            bool passes = check_dz(inner, layer2[it2], rb_proj, b_dz);
            if (passes) {
                doublet d = { it1, it2 };
                doublets[index] = d;
                ++index;
            }
        }
    }
    return index;
}

int main()
{
    while (true) {
        if (data_is_there()) {
            reset_flags();

            metadata m = *device_addr::meta;

            int written = process(m);

            doublets_status s {
                doublets_status_code::ready,
                written
            };
            *(doublets_status *)(0x8e000000 + shm_addr::status) = s;
        }
    }

    return 0;
}
