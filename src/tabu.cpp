#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <cmath>
#include <random>
#include <iomanip>
#include <functional>
#include "graph.h"
#include "tabu.h"
#include "trucks.h"
#include "parameters.h"

Parameters defaultParametersfortabu;

//potrzebne do graspa
bool compare_Candidates_for_first_grasp(std::pair<Node, double> &a, std::pair<Node, double> &b)
{
    return a.second < b.second; // rosnaco wedlug kosztow
}

//grasp z metody pierwszej
int Tabu::create_first_solution_with_grasp(Graph& graph){

    std::ofstream outputFile("wyniki.txt");
    if (!outputFile.is_open())
    {
       //std::cerr << "Nie można otworzyć pliku" << std::endl;
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


bool Tabu::not_in_Tabu(const int* operation){
    //std::cout<<"tabu: "<< operation[0] <<" "<< operation[1] <<" " <<operation[2] <<" "<< operation[3] <<" \n";
    for(int i=0;i<default_parameters.Tabu_list_size;i++){
        if( (Tabu_list[i][0] == operation[0] &&
            Tabu_list[i][1] == operation[1] &&
            Tabu_list[i][2] == operation[2] &&
            Tabu_list[i][3] == operation[3] &&
            Tabu_list[i][4] == operation[4]) ||
            //inverted
            (Tabu_list[i][0] == operation[2] &&
            Tabu_list[i][1] == operation[3] &&
            Tabu_list[i][2] == operation[0] &&
            Tabu_list[i][3] == operation[1] &&
            Tabu_list[i][4] == operation[4]) ) 
            {
                return false;
            }

            }

    return true;
}


void Tabu::add_to_Tabu(const int* operation,int iteration){
    //operation = [row1,column1,row2,column2,how_moved] | how_moved = 0 - swap_2_nodes | how_moved = 1 - multi route two opt? | how_moved = 2 - twoopt? ...
    if(not_in_Tabu(operation)){
        Tabu_list[iteration % defaultParametersfortabu.Tabu_list_size][0]= operation[0];
        Tabu_list[iteration % defaultParametersfortabu.Tabu_list_size][1]= operation[1];
        Tabu_list[iteration % defaultParametersfortabu.Tabu_list_size][2]= operation[2];
        Tabu_list[iteration % defaultParametersfortabu.Tabu_list_size][3]= operation[3];
        Tabu_list[iteration % defaultParametersfortabu.Tabu_list_size][4]= operation[4];
    }

    return;
}



bool Tabu::can_be_swaped(int idx_of_node1, int idx_of_node2, int idx_of_truck1, int idx_of_truck2, Graph& graph, std::vector<std::vector<int>>& current_routes) {
    if(idx_of_node1 == idx_of_node2 || idx_of_truck1 == idx_of_truck2) {
        return false;
    }

    Truck truck1 = graph.trucksvector[idx_of_truck1];
    Truck truck2 = graph.trucksvector[idx_of_truck2];
    
    std::vector<int> test_swap_route1 = truck1.route;
    std::vector<int> test_swap_route2 = truck2.route;

    std::swap(test_swap_route1[idx_of_node1], test_swap_route2[idx_of_node2]);


    double time_for_truck1 = is_single_route_possible(test_swap_route1, graph);
    double time_for_truck2 = is_single_route_possible(test_swap_route2, graph);

    if(time_for_truck1 == -1 || time_for_truck2 == -1) {
        return false;
    }

    graph.trucksvector[idx_of_truck2].current_time = time_for_truck2; 
    graph.trucksvector[idx_of_truck1].current_time = time_for_truck1;
    graph.trucksvector[idx_of_truck1].route = test_swap_route1;
    graph.trucksvector[idx_of_truck2].route = test_swap_route2;

    current_routes[idx_of_truck1] = test_swap_route1;
    current_routes[idx_of_truck2] = test_swap_route2;

    return true;
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
    if(cost <= graph.Nodes[0].duedate && cargo>=0){
        cost -=graph.distances[0][route[route.size()-1]];
        return cost;
    }else{
        return -1;
    }
return cost;
}


void Tabu::swap_two_nodes(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, int* used_ops, Graph& graph) {
    for (int i = 0; i < current_routes.size(); ++i) {
        for (int j = 0; j < current_routes.size(); ++j) {
            if (i == j || current_routes[i].empty() || current_routes[j].empty()) {
                continue;
            }
            for (int k = 0; k < current_routes[i].size(); ++k) {
                for (int l = 0; l < current_routes[j].size(); ++l) {
                    // Test routes after swapping
                    std::vector<int> test_route1 = current_routes[i];
                    std::vector<int> test_route2 = current_routes[j];
                    std::swap(test_route1[k], test_route2[l]);

                    // Check if the new routes are valid
                    double time_for_route1 = is_single_route_possible(test_route1, graph);
                    double time_for_route2 = is_single_route_possible(test_route2, graph);

                    if (time_for_route1 != -1 && time_for_route2 != -1) {
                        used_ops[0] = i; // Truck 1 index
                        used_ops[1] = k; // Node index in truck 1
                        used_ops[2] = j; // Truck 2 index
                        used_ops[3] = l; // Node index in truck 2
                        used_ops[4] = 0; // Operation type (swap_two_nodes)

                        // Check if the operation is not Tabu
                        if (not_in_Tabu(used_ops)) {
                            // Apply the swap
                            current_routes[i] = test_route1;
                            current_routes[j] = test_route2;

                            graph.trucksvector[i].route = test_route1;
                            graph.trucksvector[j].route = test_route2;

                            graph.trucksvector[i].current_time = time_for_route1;
                            graph.trucksvector[j].current_time = time_for_route2;

                            return;
                        }
                    }
                }
            }
        }
    }
}



void Tabu::two_opt_swap(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, int* used_ops, Graph& graph) {
    for(int i = 0; i < current_routes.size(); i++) {
        for(int j = i + 1; j < current_routes.size(); j++) {

            if(current_routes[i].size() < 2 || current_routes[j].size() < 2) continue;

            for(int k = 1; k < current_routes[i].size(); k++) {
                for(int l = 1; l < current_routes[j].size(); l++) {
                    
                    std::vector<int> test_route1 = current_routes[i];
                    std::vector<int> test_route2 = current_routes[j];


                    //wycinanie fragmentow
                    std::vector<int> segment1(test_route1.begin() + k, test_route1.end());
                    std::vector<int> segment2(test_route2.begin() + l, test_route2.end());

                    //zlaczenie w trasy
                    test_route1.erase(test_route1.begin() + k, test_route1.end());
                    test_route1.insert(test_route1.end(), segment2.begin(), segment2.end());

                    test_route2.erase(test_route2.begin() + l, test_route2.end());
                    test_route2.insert(test_route2.end(), segment1.begin(), segment1.end());
                    double time_for_route1 = is_single_route_possible(test_route1, graph);
                    double time_for_route2 = is_single_route_possible(test_route2, graph);          
                    if(time_for_route1 != -1 && time_for_route2 != -1) {


                        used_ops[0] = i;
                        used_ops[1] = k;
                        used_ops[2] = j;    
                        used_ops[3] = l;
                        used_ops[4] = 1;

                        // Check Tabu list
                        if(not_in_Tabu(used_ops)) {
                            current_routes[i] = test_route1;
                            current_routes[j] = test_route2;

                            graph.trucksvector[i].route = test_route1;
                            graph.trucksvector[j].route = test_route2;

                            graph.trucksvector[i].current_time = time_for_route1;
                            graph.trucksvector[j].current_time = time_for_route2;

                            return;
                        }
                    }
                }
            }
        }
    }
}


void Tabu::generate_neighbour(double& current_cost,int& current_used_trucks,std::vector<std::vector<int>>& current_routes,int* used_ops, Graph& graph){

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(1, 2);

    int which_operation = dist(gen);
//jak to polaczyc?
    if(which_operation == 1){
        two_opt_swap(current_cost, current_used_trucks, current_routes, used_ops, graph);
    }
    else if(which_operation == 2){
        swap_two_nodes(current_cost,current_used_trucks,current_routes,used_ops,graph);
    }

    current_cost = 0;
    current_used_trucks = 0;
    for (int i = 0; i < graph.trucksvector.size(); i++) {
        if (!graph.trucksvector[i].route.empty()){
            current_cost += graph.trucksvector[i].current_time + graph.distances[(current_routes[i][current_routes[i].size()-1])][0];
            current_used_trucks ++;
        }
    } //std::cout<<"\n";
    return;
}

void update_graph_from_routes(Graph& graph, const std::vector<std::vector<int>>& routes) {
    for (int i = 0; i < routes.size(); ++i) {
        graph.trucksvector[i].route = routes[i];
    }
}

void Tabu::Tabu_search(Graph& graph){
    auto start=std::chrono::high_resolution_clock::now();
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

   std::cout<<"starting cost: "<<first_solution_cost<<" \n ";

    //pierwsze rozwiazanie: first_solution_trucks_number, first_solution_cost, first_solution_routes
    current_solution_cost = first_solution_cost;
    current_solution_used_trucks = first_solution_trucks_number +1;
    current_solution_routes = first_solution_routes;
    int no_improvement = 0;
    int iteration = 0;
    double neighbour_solution_cost = current_solution_cost;
    int neighbour_solution_used_trucks = current_solution_used_trucks;
    std::vector<std::vector<int>> neighbour_solution_routes = current_solution_routes;
    int used_ops[5];
    Graph graphcopy = graph;
    //SZUKAMY LEPSZYCH SASIEDNICH ROZWIAZAN (NIE POZWALAJAC NA 2 TAKIE SAME JESLI POPRZEDNIE JESZCZE W TABU LIST)

    while(iteration < Max_iterations && no_improvement < no_improvement_limit){  
        generate_neighbour(neighbour_solution_cost, neighbour_solution_used_trucks, neighbour_solution_routes, used_ops, graph);

        if (neighbour_solution_cost < current_solution_cost){
            current_solution_cost = neighbour_solution_cost;
            current_solution_routes = neighbour_solution_routes;
            current_solution_used_trucks = neighbour_solution_used_trucks;
            std::cout<<"iter: "<<iteration<<" better -> "<<std::fixed << std::setprecision(5) <<current_solution_cost<<" used_operation "<<used_ops[4]<<" ikjl values: "<<used_ops[0]<<" "<<used_ops[1]<<" "<<used_ops[2]<<" "<<used_ops[3]<<" \n";
            update_graph_from_routes(graph, current_solution_routes);
            graphcopy.trucksvector = graph.trucksvector;
        }
        else{
            no_improvement ++;
            neighbour_solution_cost = current_solution_cost;
            neighbour_solution_routes = current_solution_routes;
            neighbour_solution_used_trucks = current_solution_used_trucks;
            graph.trucksvector = graphcopy.trucksvector;

        }

        add_to_Tabu(used_ops,iteration);
        for (int i = 0; i < current_solution_routes.size(); i++) {
            graph.trucksvector[i].route = current_solution_routes[i];
        }
        
        if(no_improvement >= no_improvement_limit){
            break;
        }
        iteration ++;
        
    } 









auto stop=std::chrono::high_resolution_clock::now();  // Koniec pomiaru
auto duration=std::chrono::duration_cast<std::chrono::seconds>(stop-start);
std::cout<<"\nCzas wykonania funkcji: "<<duration.count()<<" sekund"<<std::endl;

// WYPISANIE 
    std::ofstream outputFile("wyniki.txt");
    if (!outputFile.is_open())
    {
       //std::cerr << "Nie można otworzyć pliku" << std::endl;
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
    outputFile << std::fixed << std::setprecision(5) << current_solution_used_trucks << " " << current_solution_cost << std::endl;
    for (int l = 0; l < current_solution_routes.size(); l++)
    {
        if (!current_solution_routes[l].empty())
        {
            for (int k = 0; k < current_solution_routes[l].size(); k++)
            {
                outputFile << current_solution_routes[l][k] << " ";
            }
            outputFile << std::endl;
        }
    }




    outputFile.close();
}