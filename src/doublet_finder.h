#ifndef DOUBLET_FINDER_H
#define DOUBLET_FINDER_H

#include <cstdint>
#include <vector>
#include <utility>

#ifdef HARDWARE_ACCELERATOR
#else // HARDWARE_ACCELERATOR
#   include <atomic>
#endif // HARDWARE_ACCELERATOR

#include "hit.h"

class pb_doublet_finder
{
public:
    /// \brief States for the doublet finder state machine.
    using state = std::uint32_t;

    /**
     * \brief The finder is ready. One may upload hits and start it.
     *
     * Expected command: \ref command_start
     *
     * Next state: \ref state_processing
     */
    const constexpr static state state_ready = 1u << 31;

    /**
     * \brief The finder is processing hits.
     *
     * Expected command: none
     *
     * Next state: \ref state_out_of_memory or \ref state_finished
     */
    const constexpr static state state_processing = 1u << 30;

    /**
     * \brief The finder has run out of memory.
     *
     * Expected command: \ref command_done_reading
     *
     * Next state: \ref state_processing
     */
    const constexpr static state state_out_of_memory = 1u << 29;

    /**
     * \brief The finder has finished.
     *
     * Expected command: \ref command_done_reading
     *
     * Next state: \ref state_ready
     */
    const constexpr static state state_finished = 1u << 28;

    /// \brief A doublet, represented as indices within the two layers
    using doublet = std::pair<std::uint16_t, std::uint16_t>;

    /// \brief Maximum number of hits passed to \ref set_hits.
    const constexpr static std::size_t max_hit_count = 2048;

    /**
     * \brief Gets back the produced doublets.
     *
     * Production will be resumed if the producer was out of memory.
     *
     * \return The number of doublets added to \c output.
     */
    std::size_t get_doublets(std::vector<doublet> &output);

    /// \brief Retrieves the state of the machine.
    state get_state() const;

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

protected:
    /**
     * \brief Commands for the doublet finder state machine.
     */
    enum class command : std::uint32_t
    {
        /**
         * \brief Start processing.
         *
         * Allowed in state: \c ready
         *
         * Next state: \c processing
         */
        start = 1u << 31,

        /**
         * \brief Notify that the host has run all the generated doublets.
         *
         * Allowed in state: \c out_of_memory or \c finished
         *
         * Next state: \c processing or \c ready
         */
        done_reading = 1u << 30,
    };

    /// \brief Sends a command to the state machine.
    void send_command(command cmd);

private:
#ifdef HARDWARE_ACCELERATOR
#else // HARDWARE_ACCELERATOR
    compact_beam_spot _bs;
    const std::vector<compact_pb_hit> *_layer1;
    const std::vector<compact_pb_hit> *_layer2;
    std::vector<doublet> _doublets;
    std::atomic<state> _state = state_ready;
#endif // HARDWARE_ACCELERATOR
};

#endif // DOUBLET_FINDER_H
