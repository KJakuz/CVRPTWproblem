#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <cmath>
#include <random>
#include <iomanip>
#include "graph.h"
#include "tabu.h"
#include "trucks.h"
#include "parameters.h"

Parameters defaultParametersfortabu;

bool compare_Candidates_for_first_grasp(std::pair<Node, double> &a, std::pair<Node, double> &b)
{
    return a.second < b.second; // rosnaco wedlug kosztow
}

int Tabu::create_first_solution_with_grasp(Graph& graph){

    std::ofstream outputFile("wyniki.txt");
    if (!outputFile.is_open())
    {
        std::cerr << "Nie można otworzyć pliku" << std::endl;
        return -1;
    }


    graph.measure_distances();

    std::vector<std::vector<int>> solution;
    int idx_of_truck = 0;
    graph.makeunvisitedvector();

    while (!graph.all_visited())
    {

        if (idx_of_truck >= graph.trucksvector.size())
        {
            graph.trucksvector.push_back(Truck (0, idx_of_truck, graph.trucksvector[0].capacity, graph.trucksvector[0].capacity, 0, 0));
        }

        std::vector<std::pair<Node, double>> candidates; // lista kandydatow w formie <<wezel>,<koszt>>

        // sprawdzanie mozliwych kandydatow
        for (int i = 0; i < graph.unvisited.size(); i++)
        {
            if ((graph.trucksvector[idx_of_truck].capacity >= graph.unvisited[i].demand) && (graph.distances[0][graph.unvisited[i].id] < graph.unvisited[i].duedate) && (graph.unvisited[i].readytime + graph.unvisited[i].servicetime + graph.distances[0][graph.unvisited[i].id] <= graph.Nodes[0].duedate))
            {
                if (graph.trucksvector[idx_of_truck].cargo >= graph.unvisited[i].demand)
                {
                    if ((graph.trucksvector[idx_of_truck].current_time + graph.distances[graph.trucksvector[idx_of_truck].which_node][graph.unvisited[i].id] < graph.unvisited[i].duedate) && (graph.trucksvector[idx_of_truck].current_time + graph.distances[graph.trucksvector[idx_of_truck].which_node][graph.unvisited[i].id] + graph.unvisited[i].servicetime + graph.distances[0][graph.unvisited[i].id] <= graph.Nodes[0].duedate))
                    {
                        // obliczanie kosztow kandydatow
                        double waiting_time_costs = std::max(0.0, graph.unvisited[i].readytime - (graph.trucksvector[idx_of_truck].current_time + graph.distances[graph.trucksvector[idx_of_truck].which_node][graph.unvisited[i].id]));
                        double window_time_costs = std::max(0.0, graph.unvisited[i].duedate - (graph.trucksvector[idx_of_truck].current_time + graph.distances[graph.trucksvector[idx_of_truck].which_node][graph.unvisited[i].id]));
                        double demand_cost = (graph.trucksvector[idx_of_truck].capacity - graph.trucksvector[idx_of_truck].cargo) - graph.unvisited[i].demand;
                        double cost = defaultParametersfortabu.distance_cost_param * graph.distances[graph.trucksvector[idx_of_truck].which_node][graph.unvisited[i].id] +
                                      +defaultParametersfortabu.demand_param * demand_cost + defaultParametersfortabu.window_time_param * window_time_costs + defaultParametersfortabu.waiting_time_param * waiting_time_costs;
                        candidates.push_back(std::make_pair(graph.unvisited[i], cost));
                    }
                }
            }
            else
            {
                // jesli niemozliwy transport
                int failed = -1;
                outputFile << failed;
                outputFile.close();
                exit(1);
            }
        }
        // po sprawdzaniu kandydatiw
        if (candidates.size() == 0)
        {
            idx_of_truck++;
            continue;
        }

        // sortowanie kandydatow rosnaco pod wzgledem kosztu
        std::sort(candidates.begin(), candidates.end(), compare_Candidates_for_first_grasp);

        // wybieramy losowo wezel z pierwszych iles procent kandydatow
        int limit = std::max(static_cast<int>(candidates.size() * (defaultParametersfortabu.RCLpercent / 100.0)), 1);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, limit - 1);
        int next_node_index = dis(gen);

        // aktualizujemy atrybuty ciezarowki
        graph.trucksvector[idx_of_truck].route.push_back(candidates[next_node_index].first.id);
        graph.trucksvector[idx_of_truck].cargo -= candidates[next_node_index].first.demand;
        double waiting_time = std::max(0.0, candidates[next_node_index].first.readytime - (graph.trucksvector[idx_of_truck].current_time + graph.distances[graph.trucksvector[idx_of_truck].which_node][candidates[next_node_index].first.id]));
        graph.trucksvector[idx_of_truck].current_time += graph.distances[graph.trucksvector[idx_of_truck].which_node][candidates[next_node_index].first.id] + candidates[next_node_index].first.servicetime + waiting_time;
        graph.trucksvector[idx_of_truck].which_node = candidates[next_node_index].first.id;
        candidates[next_node_index].first.check_if_done = true;

        // usuniecie z graph.unvisited wezla przez ktory przejezdzamy (znajdujemy id do usuniecia w wektorze)
        graph.unvisited.erase(std::remove_if(graph.unvisited.begin(), graph.unvisited.end(), [&](Node &node)
                                       { return node.id == candidates[next_node_index].first.id; }),
                        graph.unvisited.end());
    }
    // zwracamy liczbe tirow
    return idx_of_truck;
}


bool Tabu::can_be_swaped(int idx_of_node1, int idx_of_node2,int idx_of_truck1,int idx_of_truck2,Graph& graph){
    if(idx_of_node1 !=idx_of_node2){
        Truck truck1 = graph.trucksvector[idx_of_truck1];
        Truck truck2 = graph.trucksvector[idx_of_truck2];
        Node node1 = graph.Nodes[truck1.route[idx_of_node1]];
        Node node2 = graph.Nodes[truck2.route[idx_of_node2]];
        
        truck1.route[idx_of_node1] = node2.id;
        truck2.route[idx_of_node2] = node1.id;

        double time_for_truck1 = 0, time_for_truck2 = 0,waiting_time_first_node_truck1=0,waiting_time_first_node_truck2=0;
        double cargo_for_truck1 = truck1.capacity, cargo_for_truck2 = truck2.capacity;

        //sciezka1
        time_for_truck1 = is_single_route_possible(truck1.route,graph);
        time_for_truck2 = is_single_route_possible(truck2.route,graph);

        if(time_for_truck1 == -1 || time_for_truck2 == -1){

            return false;
        
        }
        else{

            graph.trucksvector[idx_of_truck2].current_time = time_for_truck2; 
            graph.trucksvector[idx_of_truck1].current_time = time_for_truck1; 

        } 
        
    return true;
    }
    else{
        return false;
        }
}

double Tabu::is_single_route_possible(std::vector<int>& route,Graph& graph){


    Truck infotruck = graph.trucksvector[0];

    double cost = 0, waiting_time_first_node=0;
    int cargo = infotruck.capacity;


    waiting_time_first_node = std::max(0.0,(graph.Nodes[route[0]].readytime-graph.distances[0][route[0]]));


    cost += graph.distances[0][route[0]] + graph.Nodes[route[0]].servicetime + waiting_time_first_node;
    cargo -= graph.Nodes[route[0]].demand;

    for(int i=1;i<route.size();i++){
        cargo -= graph.Nodes[route[i]].demand;
        cost += graph.distances[route[i-1]][route[i]];
        cost += std::max(0.0,(graph.Nodes[route[i]].readytime-cost));
        cost += graph.Nodes[route[i]].servicetime;
        if((cost <= graph.Nodes[route[i]].duedate) && (cargo >= graph.Nodes[route[i]].demand)){

            continue;
        }
        else{
            return -1;
        }
    }

    cost += graph.distances[0][route[route.size()-1]];
    //cost += graph.Nodes[route[route.size()-1]].servicetime;
    if(cost <= graph.Nodes[0].duedate && cargo>=0){
        cost -=graph.distances[0][route[route.size()-1]];
        return cost;
    }else{
        return -1;
    }
return cost;
}


void Tabu::generate_neighbour(double& current_cost,int& current_used_trucks,std::vector<std::vector<int>>& current_routes,Graph& graph){
    for(int i=0;i<current_routes.size();i++){
        for(int j=0;j<current_routes.size();j++){
            if(i != j && current_routes[i].size()>1 && current_routes[j].size()>1 ){
                for(int k=0;k<=current_routes[i].size()-1;k++){
                    for(int l=0;l<=current_routes[j].size()-1;l++){
                        if(can_be_swaped(k,l,i,j,graph)){
                            int temp_x = current_routes[i][k];
                            current_routes[i][k] = current_routes[j][l];
                            current_routes[j][l] = temp_x;


                            //obliczanie nowego czasu trasy 'i' i 'j'
                            current_cost = 0;
                            for (int i = 0; i < graph.trucksvector.size(); i++)
                            {
                                if (!graph.trucksvector[i].route.empty()){
                                    current_cost += graph.trucksvector[i].current_time + graph.distances[(current_routes[i][current_routes[i].size()-1])][0];

                                }
                            } 
                            return;
                        }
                    }
                }
            }
        }
    }
}


void Tabu::Tabu_search(Graph& graph){
    //inicjalizacja zmiennych
    double current_solution_cost = 0;
    int current_solution_used_trucks = 0;
    std::vector<std::vector<int>> current_solution_routes;


    // TWORZENIE PIERWSZEGO ROZWIAZANIA
    std::vector<std::vector<int>> first_solution_routes;
    double first_solution_cost = 0;
    int first_solution_trucks_number = create_first_solution_with_grasp(graph); 

    for (int i = 0; i < graph.trucksvector.size(); i++)
        {
            first_solution_cost += graph.trucksvector[i].current_time + graph.distances[graph.trucksvector[i].which_node][0];
        }

    for (int l = 0; l < graph.trucksvector.size(); l++)
        {
            first_solution_routes.push_back(graph.trucksvector[l].route);
        }

    //pierwsze rozwiazanie: first_solution_trucks_number, first_solution_cost, first_solution_routes
    current_solution_cost = first_solution_cost;
    current_solution_used_trucks = first_solution_trucks_number;
    current_solution_routes = first_solution_routes;
    int no_improvement = 0;
    int iteration = 0;
    double neighbour_solution_cost = current_solution_cost;
    int neighbour_solution_used_trucks = current_solution_used_trucks;
    std::vector<std::vector<int>> neighbour_solution_routes = current_solution_routes;
    //SZUKAMY LEPSZYCH SASIEDNICH ROZWIAZAN (NIE POZWALAJAC NA 2 TAKIE SAME JESLI POPRZEDNIE JESZCZE W TABU LIST)

    while(iteration < Max_iterations && no_improvement < no_improvement_limit){
        //tworzenie sąsiednich rozwiązań

        generate_neighbour(neighbour_solution_cost, neighbour_solution_used_trucks, neighbour_solution_routes, graph);





        iteration ++;
        
    } 













// WYPISANIE 
    std::ofstream outputFile("wyniki.txt");
    if (!outputFile.is_open())
    {
        std::cerr << "Nie można otworzyć pliku" << std::endl;
        return;
    }
/*
    outputFile << std::fixed << std::setprecision(5) << first_solution_trucks_number + 1 << " " << first_solution_cost << std::endl;
    for (int l = 0; l < first_solution_routes.size(); l++)
    {
        if (!first_solution_routes[l].empty())
        {
            for (int k = 0; k < first_solution_routes[l].size(); k++)
            {
                outputFile << first_solution_routes[l][k] << " ";
            }
            outputFile << std::endl;
        }
    }
    */
    outputFile << std::fixed << std::setprecision(5) << neighbour_solution_used_trucks + 1 << " " << neighbour_solution_cost << std::endl;
    for (int l = 0; l < neighbour_solution_routes.size(); l++)
    {
        if (!neighbour_solution_routes[l].empty())
        {
            for (int k = 0; k < neighbour_solution_routes[l].size(); k++)
            {
                outputFile << neighbour_solution_routes[l][k] << " ";
            }
            outputFile << std::endl;
        }
    }




    outputFile.close();
}