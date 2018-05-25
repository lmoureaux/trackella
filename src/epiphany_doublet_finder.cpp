#include "epiphany_doublet_finder.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

#include <e-hal.h>
#include <e-loader.h>

#include "epiphany_protocol.h"
#include "fast_sincos.h"

extern "C" e_platform_t e_platform;

epiphany_doublet_finder::epiphany_doublet_finder()
{
    e_set_loader_verbosity(L_D0);
    e_set_host_verbosity(H_D0);

    e_init(NULL);
    e_reset_system();

    e_alloc(&_eram, 0x0, 0x80000);
    e_open(&_device, 0, 0, e_platform.rows, e_platform.cols);

    e_load_group("doublet_finder_device", &_device, 0, 0, 1, 1, E_TRUE);
}

epiphany_doublet_finder::~epiphany_doublet_finder()
{
    e_close(&_device);
    e_free(&_eram);
    e_finalize();
}

std::vector<epiphany_doublet_finder::hit_type> epiphany_doublet_finder::convert(
        const std::vector<hit> &hits, int layer) const
{
    std::vector<hit_type> res;
    res.reserve(hits.size());
    for (const hit &h : hits) {
        res.emplace_back(h, 1);
    }
    return res;
}

std::size_t epiphany_doublet_finder::get_doublets(
    std::vector<epiphany_doublet_finder::doublet_type> &output)
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

void epiphany_doublet_finder::sort_hits(
        std::vector<epiphany_doublet_finder::hit_type> &layer1,
        std::vector<epiphany_doublet_finder::hit_type> &layer2)
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

void epiphany_doublet_finder::upload(
        const epiphany_doublet_finder::beam_spot_type &bs,
        const std::vector<epiphany_doublet_finder::hit_type> &layer1,
        const std::vector<epiphany_doublet_finder::hit_type> &layer2)
{
    // Send the data to the device
    metadata m {
        std::int32_t(layer1.size()), std::int32_t(layer2.size()), bs
    };
    std::cout << "metadata" << std::endl;
    e_write(&_device, 0, 0, (off_t) device_addr::meta, &m, sizeof(m));

    std::cout << "L1 "
        << device_addr::data
        << " .. "
        << (device_addr::data + layer1.size())
        << std::endl;
    e_write(&_device, 0, 0, (off_t) device_addr::data,
            layer1.data(), layer1.size() * sizeof(compact_hit));
    std::cout << "L2 "
        << (device_addr::data + layer1.size())
        << " .. "
        << (device_addr::data + layer1.size() + layer2.size() * sizeof(compact_hit))
        << std::endl;
    e_write(&_device, 0, 0,
            (off_t) (device_addr::data + layer1.size()),
            layer2.data(), layer2.size() * sizeof(compact_hit));

    // Notify the device so it can start
    int data_is_there = 1;
    std::cout << "starting" << std::endl;
    e_write(&_device, 0, 0, (off_t) device_addr::status,
            &data_is_there, sizeof(int));
}

void epiphany_doublet_finder::find(
        const epiphany_doublet_finder::beam_spot_type &bs,
        const std::vector<epiphany_doublet_finder::hit_type> &layer1,
        const std::vector<epiphany_doublet_finder::hit_type> &layer2)
{
    if (layer1.empty() || layer2.empty()) {
        return;
    }

    // Reset the flags in shared memory
    doublets_status s { doublets_status_code::invalid, 0 };
    e_write(&_eram, 0, 0, (off_t) shm_addr::status, &s, sizeof(s));

    // Start remote computing
    upload(bs, layer1, layer2);

    // Wait for the results to appear
    std::cout << "waiting" << std::endl;
    do {
        e_read(&_eram, 0, 0, (off_t) shm_addr::status, &s, sizeof(s));
//         std::cout << (int) s.code << " " << s.count << std::endl;
    } while (s.code != doublets_status_code::ready);

    // The results are there
    _doublets.resize(s.count);
    std::cout << "loading" << std::endl;
    e_read(&_eram, 0, 0, (off_t) shm_addr::doublets,
           _doublets.data(), s.count * sizeof(doublet));

    std::cout << "---" << std::endl;
//     for (auto &doublet : _doublets) {
//         std::cout << doublet.first << " " << doublet.second << std::endl;
//     }
}
