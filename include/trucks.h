#ifndef TRUCKS_H
#define TRUCKS_H
#include "graph.h"

class Node;

class Truck{
    public:
        int trucks_number, capacity;
        int cargo; //ile jeszcze ma na pace
        int which_node; // gdzie jest
        int current_time; //aktualny czas/ może float?


        bool check_time(const Node& node, const std::vector<std::vector<float>>& distances);
        bool check_cargo(const Node& node);
        void show_trucks_info();


};
#endif