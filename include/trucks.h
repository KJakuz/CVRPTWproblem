#ifndef TRUCKS_H
#define TRUCKS_H
#include "graph.h"
#include <fstream>
#include <vector>

class Node;

class Truck
{
public:
    int trucks_number, capacity, id;
    int cargo;           // ile jeszcze ma na pace
    int which_node;      // gdzie jest
    double current_time; // aktualny czas
    std::vector<int> route;

    Truck(int trucks_number, int id, int capacity, int cargo, int which_node, int current_time) : trucks_number(trucks_number), id(id),
                                                                                                  capacity(capacity), cargo(cargo), which_node(which_node), current_time(current_time) {}

    Truck() : id(-1), capacity(-1), cargo(-1), which_node(-1), current_time(-1) {}

    void show_trucks_info();
};
#endif