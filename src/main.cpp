#include <iostream>
#include <string>
#include "graph.h"
#include "trucks.h"
#include "loadfile.h"
#include "tabu.h"

extern Parameters defaultParametersfortabu;

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
    Graph graf = graph;

    double temperature[4]={10000,500,1000,5000};
    double coolingfactor[3]={0.99,0.945,0.93};

    for(int i=0;i<2;i++){
        for(int j=1;j<4;j++){
            for(int s=0;s<3;s++){
                graph = graf;
                //defaultParametersfortabu.temperature=temperature[0];
                //defaultParametersfortabu.cooling_factor=coolingfactor[s];
                //graph.rungrasp();
                tabu.Tabu_search(graph);

            }
        }
        

    }



    return 0;
}