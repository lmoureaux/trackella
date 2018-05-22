#ifndef DOUBLET_FINDER_H
#define DOUBLET_FINDER_H

#include <cstdint>
#include <vector>
#include <utility>

#include "hit.h"

class pb_doublet_finder
{
public:
    /// \brief A doublet, represented as indices within the two layers
    using doublet = std::pair<std::uint16_t, std::uint16_t>;

    /**
     * \brief Gets back the produced doublets.
     *
     * Production will be resumed if the producer was out of memory.
     *
     * \return The number of doublets added to \c output.
     */
    std::size_t get_doublets(std::vector<doublet> &output);

    /**
     * \brief Sets the beam spot properties.
     */
    void set_beam_spot(const compact_beam_spot &bs);

    /**
     * \brief Pushes hits to the machine.
     *
     * The vectors have to remain valid until all results have been read.
     */
    void set_hits(const std::vector<compact_pb_hit> &layer1,
                  const std::vector<compact_pb_hit> &layer2);

    /// \brief Starts producing doublets.
    void start();

private:
    compact_beam_spot _bs;
    const std::vector<compact_pb_hit> *_layer1;
    const std::vector<compact_pb_hit> *_layer2;
    std::vector<doublet> _doublets;
};

#endif // DOUBLET_FINDER_H
