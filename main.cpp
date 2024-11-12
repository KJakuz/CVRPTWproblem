#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include "graph.h"
#include "trucks.h"
#include "loadfile.h"
#include "parameters.h"

int main(int argc, char* argv[]){
    Graph graph;
    Truck truck;
    LoadFile loader;

    extern Parameters parameters;

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


    //grid search
    int distanceparams[5]={1,5,10,15,20};
    int waitingtimeparams[5]={1,5,10,15,20};
    int demandparams[5]={1,5,10,15,20};

    int trucks=0;
    double distance=0;


    //graph.rungrasp();
for(int s=0;s<3;s++){
    for(int i = 0; i < 5; i++) {
            for(int k = 0; k < 5; k++) {
                for(int j=0;j < 5;j++){
                    parameters.distance_cost_param = distanceparams[i];
                    parameters.waiting_time_param = waitingtimeparams[k];
                    parameters.demand_param = demandparams[j];
                    trucks, distance = graph.rungrasp();
                    //outputFile <<parameters.distance_cost_param<<";"<< parameters.window_time_param<<";"<< parameters.waiting_time_param<<";"<< parameters.RCLpercent <<";";
                    //outputFile << std::fixed << std::setprecision(5)<< trucks + 1 << ";" << distance << std::endl;
                    //outputFile <<"siema";
            }
            }
    }
}

    return 0;
}