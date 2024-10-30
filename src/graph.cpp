#include <iostream>
#include <cmath>
#include "graph.h"

void Graph::show_number_nodes(){
            std::cout<<"123";
        std::cout<<"numer węzła: "<<number_nodes<<" | "<<Nodes[1].ycord<<std::endl;
}


    //liczenie macierzy odleglosci od węzłów
void Graph::measure_distances(){
    distances.resize(number_nodes, std::vector<float>(number_nodes, 0.0f)); //alokacja pamiec, wielkosc macierzy number_nodes i zapisanej zerami
    for(int i=0;i<number_nodes;i++){
        for(int j=0;j<number_nodes;j++){
            //wzor na odleglosc na plaszczyznie euklidesowej
            distances[i][j]=sqrt( pow(Nodes[j].xcord - Nodes[i].xcord,2) + pow(Nodes[j].ycord - Nodes[i].ycord,2) );
        }
    }
}


void Graph::show_distances_matrix(){
    std::cout<<"macierz odleglosci: \n   ";
    for(int k=0;k<number_nodes;k++){
        std::cout<<k<<"  ";
    }
    std::cout<<std::endl;
    for(int i=0;i<number_nodes;i++){
        std::cout<<i<<": ";
    for(int j=0;j<number_nodes;j++){
        std::cout<<distances[i][j]<<" ";
        //std::cout<<"distances["<<i<<"]"<<"["<<j<<"]: "<<distances[i][j]<<"\n";
    }
    std::cout<<"\n";
    }
}


void Graph::show_nodes_values(){
    for (int i = 0; i < number_nodes; i++) {
        Node& node = Nodes[i];
        std::cout << "ID: " << node.id << ", X: " << node.xcord
            << ", Y: " << node.ycord
            << ", Zapotrzebowanie: " << node.demand
            << ", Czas otwarcia: " << node.readytime
            << ", Czas zamknięcia: " << node.duedate
            << ", Czas obsługi: " << node.servicetime << std::endl;
    }
}