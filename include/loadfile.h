#ifndef LOADFILE_H
#define LOADFILE_H
#include <string>
#include "trucks.h"
#include "graph.h"

class LoadFile {
public:
    bool load(const std::string& filename, Graph& graph, Truck& truck); 
};

#endif
