#ifndef EPIPHANY_DOUBLET_FINDER_H
#define EPIPHANY_DOUBLET_FINDER_H

#include <cstdint>
#include <vector>
#include <utility>

#include <e-hal.h>

#include "compact.h"

class epiphany_doublet_finder
{
public:
    /// \brief A doublet, represented as indices within the two layers
    using doublet_type = std::pair<std::uint16_t, std::uint16_t>;

    /// \brief The type to use for hits
    using hit_type = compact_hit;

    /// \brief The type to use for the beam spot
    using beam_spot_type = compact_beam_spot;

    /// \brief Number of cores to use
    const constexpr static int CORES = 16;

    /// \brief Size of each output buffer
    const constexpr static int OUTPUT_BUFFER_SIZE
        = (1024 + 512) * sizeof(doublet_type) / CORES;

    /// \brief Constructor
    explicit epiphany_doublet_finder();

    /// \brief Destructor
    virtual ~epiphany_doublet_finder();

    /// \brief Convert hits to the correct representation
    std::vector<hit_type> convert(const std::vector<hit> &hits, int layer) const;

    /// \brief Convert beam spot info to the correct representation
    beam_spot_type convert(const beam_spot &bs) const
    {
        return {
            length_to_compact<std::int32_t>(bs.r),
            length_to_compact<std::int32_t>(bs.z),
            radians_to_compact(bs.phi)
        };
    }

    /**
     * \brief Gets back the produced doublets.
     *
     * Production will be resumed if the producer was out of memory.
     *
     * \return The number of doublets added to \c output.
     */
    std::size_t get_doublets(std::vector<doublet_type> &output);

    /**
     * \brief Pushes hits to the machine.
     *
     * You are responsible for passing back the vectors to \ref set_hits.
     */
    void sort_hits(std::vector<hit_type> &layer1,
                   std::vector<hit_type> &layer2);

    /**
     * \brief Finds doublets.
     */
    void find(const beam_spot_type &bs,
              const std::vector<hit_type> &layer1,
              const std::vector<hit_type> &layer2);

private:
    /// \brief Uploads data to the device
    void upload(const epiphany_doublet_finder::beam_spot_type &bs,
                const std::vector<hit_type> &layer1,
                const std::vector<hit_type> &layer2);

    e_epiphany_t _device;
    e_mem_t _eram;
    std::array<e_mem_t, CORES> _results_buffers;
    std::vector<doublet_type> _doublets;
};

#endif // EPIPHANY_DOUBLET_FINDER_H
