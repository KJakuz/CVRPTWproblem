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
    int distanceparams[4]={1,20,40,60};
    int waitingtimeparams[4]={1,20,40,60};
    int demandparams[4]={20,60,120,300};
    int windowparams[4]={1,20,40,60};
    int timeparam[7]={10,30,50,100,150,300,500};

    int trucks=0;
    double distance=0;
    for(int i=0;i<3;i++){
        for(int s=5;s<7;s++){
            parameters.time_limit_in_seconds = timeparam[s];
            graph.rungrasp();
        }
    }


    /*
        //graph.rungrasp();
    for(int s=0;s<3;s++){
        for(int i = 0; i < 4; i++) {
                for(int k = 0; k < 4; k++) {
                    for(int j=0;j < 4;j++){
                        for(int m=0;m < 4;m++){
                        parameters.distance_cost_param = distanceparams[i];
                        parameters.waiting_time_param = waitingtimeparams[k];
                        parameters.window_time_param = windowparams[j];
                        parameters.demand_param = demandparams[m];
                        trucks, distance = graph.rungrasp();
                        //outputFile <<parameters.distance_cost_param<<";"<< parameters.window_time_param<<";"<< parameters.waiting_time_param<<";"<< parameters.RCLpercent <<";";
                        //outputFile << std::fixed << std::setprecision(5)<< trucks + 1 << ";" << distance << std::endl;
                        //outputFile <<"siema";
                    }
                    }
                }
        }
    }*/

    return 0;
}