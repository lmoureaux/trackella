#include "doublet_finder.h"

#include <cassert>

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

void pb_doublet_finder::start()
{
    send_command(command::start);
}

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
                _doublets.push_back({ i1, i2 });
            }
        }

        _state ^= state_processing;
        _state |= state_finished;
    }
#endif // HARDWARE_ACCELERATOR
}
