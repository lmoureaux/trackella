#ifndef EVENT_READER_H
#define EVENT_READER_H

#include <memory>
#include <string>

#include "event.h"

class TFile;

class event_reader final
{
    struct data;
    std::unique_ptr<data> _d;

public:
    explicit event_reader(const std::string &filename);
    ~event_reader();

    bool next();
    std::unique_ptr<event> get();
};

#endif // EVENT_READER_H
