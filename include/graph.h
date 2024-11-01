#ifndef GRAPH_H
#define GRAPH_H
#include <vector>
#include "trucks.h"

class Truck;

class Node{
    public:
        int id; //numer węzła (0 to magazyn)
        int xcord, ycord, demand, readytime, duedate, servicetime; //koordynaty, ile potrzeba, kiedy sie otwiera, kiedy zamyka, ile czasu na rozladunek
        bool check_if_done; //czy juz dostal towar

};

class Graph{
    public:
        std::vector<Node> Nodes;
        int number_nodes;
        std::vector<std::vector<float>> distances; //tablica z dlugosciami pomiedzy węzłami
        std::vector<Node> unvisited;
        std::vector<Truck> trucksvector;

        void init_trucks(Truck& truck);
        void measure_distances();
        void show_distances_matrix();
        void show_nodes_values();
        void show_number_nodes();
        void unvisitedmap();
        void GRASP(int alfa,int beta,int gamma);
        bool all_visited();
        void create_trucks();
        void makeunvisitedvector();
        void show_one_node_values(Node& node);
        //...
        
};

#endif