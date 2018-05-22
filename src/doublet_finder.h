#ifndef DOUBLET_FINDER_H
#define DOUBLET_FINDER_H

#include <cstdint>
#include <vector>
#include <utility>

#include "compact.h"

class cpu_doublet_finder
{
public:
    /// \brief A doublet, represented as indices within the two layers
    using doublet = std::pair<std::uint16_t, std::uint16_t>;

    /// \brief The type to use for hits
    using hit_type = compact_hit;

    /// \brief The type to use for the beam spot
    using beam_spot_type = compact_beam_spot;

    /**
     * \brief Gets back the produced doublets.
     *
     * Production will be resumed if the producer was out of memory.
     *
     * \return The number of doublets added to \c output.
     */
    std::size_t get_doublets(std::vector<doublet> &output);

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
    std::vector<doublet> _doublets;
};

#endif // DOUBLET_FINDER_H
