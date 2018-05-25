#ifndef EPIPHANY_PROTOCOL_H
#define EPIPHANY_PROTOCOL_H

#include "compact.h"

#ifdef ON_DEVICE
#   define CONST_ON_DEVICE const
#   define CONST_ON_CPU
#else
#   define CONST_ON_DEVICE
#   define CONST_ON_CPU const
#endif // ON_DEVICE

#pragma pack(push,2)
struct metadata
{
    std::int32_t layer1_count, layer2_count;
    compact_beam_spot bs;
};

enum class doublets_status_code : std::int32_t
{
    invalid, working, ready,
};

struct doublets_status
{
    doublets_status_code code;
    std::int32_t count;
};
#pragma pack(pop)

using doublet = std::pair<std::int16_t, std::int16_t>;
static_assert(sizeof(doublet) == 4, "doublet size");

namespace device_addr
{
    volatile int * const status = (volatile int *) 0x2000;
    CONST_ON_DEVICE metadata * const meta = (metadata *)(status + 1);
    CONST_ON_DEVICE int * const start_flag = (int *)(meta + 1);
    CONST_ON_DEVICE compact_hit * const data = (compact_hit *)(start_flag + 1);
}

namespace shm_addr
{
    const constexpr std::ptrdiff_t status   = 0x00;
    const constexpr std::ptrdiff_t doublets = 0x10;
}

#endif // EPIPHANY_PROTOCOL_H
