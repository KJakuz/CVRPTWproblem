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
    int distanceparams[4]={1,15,25,35};
    int windowtimeparams[4]={1,15,25,35};
    int waitingtimeparams[4]={1,15,25,35};
    int RLCpercentparams[4]={1,2,4,5};
    int trucks=0;
    float distance=0;


    graph.rungrasp();
    /*
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            for(int k = 0; k < 4; k++) {
                for(int l = 0; l < 4; l++) {
                    parameters.distance_cost_param = distanceparams[i];
                    parameters.window_time_param = windowtimeparams[j];
                    parameters.waiting_time_param = waitingtimeparams[k];
                    parameters.RCLpercent = RLCpercentparams[l];
                    trucks, distance = graph.rungrasp();
                    //outputFile <<parameters.distance_cost_param<<";"<< parameters.window_time_param<<";"<< parameters.waiting_time_param<<";"<< parameters.RCLpercent <<";";
                    //outputFile << std::fixed << std::setprecision(5)<< trucks + 1 << ";" << distance << std::endl;
                    //outputFile <<"siema";

                }
            }
        }
    }
    */
    

    return 0;
}