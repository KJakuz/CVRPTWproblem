#ifndef GRAPH_H
#define GRAPH_H
#include <vector>
#include "trucks.h"


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


        void measure_distances();
        void show_distances_matrix();
        void show_nodes_values();
        void show_number_nodes();
        //...
        
};

#endif