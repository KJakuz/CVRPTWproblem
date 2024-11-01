#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
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

bool compareCandidates( std::pair<Node,float>& a, std::pair<Node,float>& b) {
    return a.second < b.second; // rosnaco wedlug kosztow
}

///*
void Graph::GRASP(int alfa,int beta,int gamma){
    std::ofstream outputFile("program.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Nie można otworzyć pliku" << std::endl;
        return;
    }
    //zapewnienie losowosci 
    srand(static_cast<unsigned int>(time(nullptr)));

    std::vector<std::vector<int>> solution;
    int counter = 0;
    makeunvisitedvector();
    while(!all_visited()){

        if(counter >= trucksvector.size()) {
            trucksvector.push_back(Truck(0,counter,trucksvector[0].capacity,trucksvector[0].capacity,0,0));
        }

        std::vector<std::pair<Node,float>> candidates; // lista kandydatow w formie <<wezel>,<koszt>> 
        for(int i=0;i<unvisited.size();i++){
            if((trucksvector[counter].capacity>=unvisited[i].demand) && (distances[0][unvisited[i].id]<unvisited[i].duedate)){
                if(trucksvector[counter].which_node!=unvisited[i].id){
                    if(trucksvector[counter].cargo>=unvisited[i].demand){
                        if(trucksvector[counter].current_time + distances[trucksvector[counter].which_node][unvisited[i].id] < unvisited[i].duedate){
                            //obliczanie kosztow
                            float waiting_time = std::max(0.0f,Nodes[i].readytime-(trucksvector[counter].current_time + distances[trucksvector[counter].which_node][unvisited[i].id]));
                            float window_time = std::max(0.0f,Nodes[i].duedate-(trucksvector[counter].current_time + distances[trucksvector[counter].which_node][unvisited[i].id]));
                            float cost = alfa*distances[trucksvector[counter].which_node][unvisited[i].id] + beta*window_time + gamma*waiting_time;
                            candidates.push_back(std::make_pair(unvisited[i],cost));

                        }
                    }
                }
            }
            else{
                std::cout<<distances[trucksvector[counter].which_node][unvisited[i].id]<<" "<<unvisited[i].duedate;
                int failed=-1;
                outputFile<<failed;
                outputFile.close();
                exit(1);
                
            }
            //w for
        }
        //po for kandydatow
        if(candidates.size()==0){
            counter++;
            continue;
        }
        //sortowanie kandydatow rosnaco pod wzgledem kosztu
        std::sort(candidates.begin(),candidates.end(),compareCandidates);
        
        //wybieramy losowo z pierwszych 10% kandydatow
        int limit = std::max(static_cast<int>(candidates.size()/10),1);
        int next_node_index = std::rand() % limit;

        //OCZEKWIANIE !

        trucksvector[counter].route.push_back(candidates[next_node_index].first.id);
        trucksvector[counter].cargo -= candidates[next_node_index].first.demand;
        trucksvector[counter].current_time += distances[trucksvector[counter].which_node][candidates[next_node_index].first.id] + candidates[next_node_index].first.servicetime;
        trucksvector[counter].which_node = candidates[next_node_index].first.id;
        candidates[next_node_index].first.check_if_done=true;
        
        //usuniecie z unvisited wezla przez ktory przejezdzamy / nie wiem co to
        unvisited.erase(std::remove_if(unvisited.begin(), unvisited.end(),[&](Node& node){ return node.id == candidates[next_node_index].first.id; }),unvisited.end());

    }

    float alltime = 0;
    for(int i=0;i<trucksvector.size();i++){
        alltime += (trucksvector[i].current_time + distances[trucksvector[i].which_node][0]);
    }
    outputFile << counter + 1 << " " << alltime << std::endl;
    for(int l=0; l<trucksvector.size(); l++) {
        trucksvector[l].show_route(outputFile);
    }

    outputFile.close();

}
