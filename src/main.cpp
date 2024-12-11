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

    double temperature[5]={10000,5000,1000,500,100};
    double coolingfactor[5]={0.99,0.975,0.96,0.945,0.93};
    int wagi[3]={1,3,5};
/*
    for(int i=0;i<2;i++){
        for(int j=1;j<3;j++){
            for(int s=0;s<3;s++){
                for(int k=0;k<3;k++){
                    for(int l=0;l<3;l++){
                        for(int o=0;o<3;o++){
                            graph = graf;
                            //defaultParametersfortabu.op1 = wagi[l];
                            //defaultParametersfortabu.op2 = wagi[k];
                            //defaultParametersfortabu.op3 = wagi[j];
                            //defaultParametersfortabu.op4 = wagi[s];
                            //defaultParametersfortabu.op5 = wagi[o];



                //defaultParametersfortabu.temperature=temperature[0];
                //defaultParametersfortabu.cooling_factor=coolingfactor[s];
                //graph.rungrasp();
                tabu.Tabu_search(graph);
                        }}}
            }
        }
  

    }
    
    */
  /*
                for(int k=0;k<5;k++){
                    for(int l=0;l<5;l++){
                        for(int o=0;o<5;o++){
                            graph = graf;
                            defaultParametersfortabu.temperature=temperature[l];
                            defaultParametersfortabu.cooling_factor=coolingfactor[o];
                            tabu.Tabu_search(graph);
                        }}}
*/



    tabu.Tabu_search(graph);
    return 0;
}