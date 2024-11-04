#include <iostream>
#include <string>
#include "graph.h"
#include "trucks.h"
#include "loadfile.h"

int main(int argc, char* argv[]){
    Graph graph;
    Truck truck;
    LoadFile loader;

    if (argc < 2) {
        std::cerr << "Użycie: " << argv[0] << " <nazwa_pliku>" << std::endl;
        return 1;
    }

    // Wczytywanie danych z pliku z argumentu programu
    std::string filename = argv[1];
    if (!loader.load(filename, graph, truck)) {
        std::cerr << "Nie wczytano pliku: " << filename << std::endl;
        return 1;
    }

    graph.init_trucks(truck);
    //liczenie odleglosci
    graph.measure_distances();

    // wyswietlanie pomocniczych wartosci
    //truck.show_trucks_info();
    //graph.show_nodes_values();
    //graph.show_distances_matrix();
    //...
    //graph.trucksvector[0].show_trucks_info();

    graph.rungrasp();


    return 0;
}