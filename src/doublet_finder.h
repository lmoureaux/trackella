#ifndef DOUBLET_FINDER_H
#define DOUBLET_FINDER_H

#include <chrono>
#include <cstdint>
#include <vector>
#include <utility>

#include "compact.h"

template<class FinderImpl>
class doublet_finder_wrapper
{
public:
    using finder_type = FinderImpl;
    using doublet_type = typename finder_type::doublet_type;
    using hit_type = typename finder_type::hit_type;

    using clock_type = std::chrono::high_resolution_clock;
    using duration_type = std::chrono::duration<double>;

    struct finding_results {
        duration_type formatting, sorting, finding, total;
        std::vector<doublet_type> doublets;
    };

    std::vector<hit_type> layer1;
    std::vector<hit_type> layer2;

    finding_results find(const beam_spot &bs,
                         std::array<std::vector<hit>, 4> &hits_per_layer);
};

template<class FinderImpl>
typename doublet_finder_wrapper<FinderImpl>::finding_results
    doublet_finder_wrapper<FinderImpl>::find(
        const beam_spot &bs,
        std::array<std::vector<hit>, 4> &hits_per_layer)
{
    finder_type finder;
    finding_results r;

    auto start = clock_type::now();

    auto converted_bs = finder.convert(bs);
    layer1 = finder.convert(hits_per_layer[0], 0);
    layer2 = finder.convert(hits_per_layer[1], 1);

    r.formatting = clock_type::now() - start;
    auto sorting_start = clock_type::now();

    finder.sort_hits(layer1, layer2);

    r.sorting = clock_type::now() - sorting_start;
    auto finding_start = clock_type::now();

    finder.find(converted_bs, layer1, layer2);

    finder.get_doublets(r.doublets);

    auto end = clock_type::now();
    r.finding = end - finding_start;
    r.total = end - start;

    return r;
}

class cpu_doublet_finder
{
public:
    /// \brief A doublet, represented as indices within the two layers
    using doublet_type = std::pair<std::uint16_t, std::uint16_t>;

    /// \brief The type to use for hits
    using hit_type = compact_hit;

    /// \brief The type to use for the beam spot
    using beam_spot_type = compact_beam_spot;

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
    std::vector<doublet_type> _doublets;
};

#endif // DOUBLET_FINDER_H
