#include <iostream>
#include <string>
#include "graph.h"
#include "trucks.h"
#include "loadfile.h"
#include "tabu.h"

int main(int argc, char *argv[])
{
    Graph graph;
    Truck truck;
    Tabu tabu;
    LoadFile loader;

    if (argc < 2)
    {
        std::cerr << "Użycie: " << argv[0] << " <nazwa_pliku>" << std::endl;
        return 1;
    }

    // Wczytywanie danych z pliku z argumentu programu
    std::string filename = argv[1];
    if (!loader.load(filename, graph, truck))
    {
        std::cerr << "Nie wczytano pliku: " << filename << std::endl;
        return 1;
    }

    graph.init_trucks(truck);
    graph.measure_distances();
    //graph.rungrasp();
    tabu.Tabu_search(graph);

    return 0;
}