#include <iostream>
#include <cmath>
#include <vector>
#include "graph.h"

class Truck;

//tworzymy wektor obiektow ciezarowek
void Graph::init_trucks(Truck& truckinfo){
    for(int i=0;i<truckinfo.trucks_number+1;i++){
        trucksvector.push_back(Truck(truckinfo.trucks_number,i,truckinfo.capacity,truckinfo.capacity,0,0));
    }

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

void Graph::show_one_node_values(Node& node){
    std::cout << "ID: " << node.id << ", X: " << node.xcord
            << ", Y: " << node.ycord
            << ", Zapotrzebowanie: " << node.demand
            << ", Czas otwarcia: " << node.readytime
            << ", Czas zamknięcia: " << node.duedate
            << ", Czas obsługi: " << node.servicetime << std::endl;
};


bool Graph::all_visited(){
    if(unvisited.size() == 0){
        return true;
    }
    return false;
}

void Graph::makeunvisitedvector(){
    for(int i = 1; i < number_nodes; i++){
        if(Nodes[i].check_if_done==false){
            unvisited.push_back(Nodes[i]);
        }
    }
}


///*
void Graph::GRASP(){
    std::vector<std::vector<int>> solution;
    int counter=0;
    makeunvisitedvector();
    while(!all_visited()){
        std::vector<Node> candidates; // lista kandydatow w formie <<id,odleglosc,czas>,<id,odleglosc,czas>> 
        for(int i=0;i<unvisited.size();i++){
            if(trucksvector[counter].which_node!=unvisited[i].id){
                if(trucksvector[counter].cargo>=unvisited[i].demand){
                    if(trucksvector[counter].current_time + distances[trucksvector[counter].which_node][unvisited[i].id] <= unvisited[i].duedate){
                        candidates.push_back(unvisited[i]);
                        //logika usuwania z unvisited visited wezlow (chyba wydajniejsze niz zawsze robienie nowej listy, tj. przegladania wezlow if done)
                        unvisited.erase(unvisited.begin() + i);
                        i--;
                    }
                }
            }
            else{
                continue;
            }
            //w for
        }
        //po for
        for(int k=0;k<candidates.size();k++){
            show_one_node_values(candidates[k]);
        }

    }

}
