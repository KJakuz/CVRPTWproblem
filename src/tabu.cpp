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
int iteration =0;

//potrzebne do graspa
bool compare_Candidates_for_first_grasp(std::pair<Node, double> &a, std::pair<Node, double> &b)
{
    return a.second < b.second; // rosnaco wedlug kosztow
}

//grasp z metody pierwszej
int Tabu::create_first_solution(Graph& graph){

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
            Tabu_list[i][4] == operation[4])) 
            {
                return false;
            }

            }

    return true;
}


void Tabu::add_to_Tabu(const int* operation){
    // operation = [row1,column1,row2,column2,how_moved] how_moved = 0 - swap_2_nodes, 1 - multi route two opt, 2 - twoopt, etc.
        for(int i = defaultParametersfortabu.Tabu_list_size - 1; i > 0; i--){
            Tabu_list[i][0]= Tabu_list[i-1][0];
           Tabu_list[i][1]= Tabu_list[i-1][1];
            Tabu_list[i][2]= Tabu_list[i-1][2];
            Tabu_list[i][3]= Tabu_list[i-1][3];
            Tabu_list[i][4]= Tabu_list[i-1][4];
        }
        Tabu_list[0][0]= operation[0];
        Tabu_list[0][1]= operation[1];
        Tabu_list[0][2]= operation[2];
        Tabu_list[0][3]= operation[3];
        Tabu_list[0][4]= operation[4];

    
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


void update_graph_from_routes(Graph& graph, const std::vector<std::vector<int>>& routes) {
    for (int i = 0; i < routes.size(); ++i) {
        graph.trucksvector[i].route = routes[i];
    }
}


void Tabu::restart_solution(std::vector<std::vector<int>>& restart_solution_routes, double& restart_solution_cost, int& restart_solution_trucks_number, Graph& graph){
    restart_solution_routes.clear();
    restart_solution_cost = 0;
    restart_solution_trucks_number = create_first_solution(graph) +1; 

    for (int i = 0; i < graph.trucksvector.size(); i++)
        {
            restart_solution_cost += graph.trucksvector[i].current_time + graph.distances[graph.trucksvector[i].which_node][0];
        }

    for (int l = 0; l < graph.trucksvector.size(); l++)
        {
            restart_solution_routes.push_back(graph.trucksvector[l].route);
        }

   std::cout<<"\nrestarting cost: "<<restart_solution_trucks_number<<" "<<restart_solution_cost<<" \n ";

}


bool accept_worse_solution(double delta_cost, double temperature) {
    if (delta_cost < 0) return true; // Akceptacja lepszych rozwiązań
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0, 1);
    return dist(gen) < std::exp(-delta_cost / temperature);
}


void Tabu::Tabu_search(Graph& graph){
    
    auto start=std::chrono::high_resolution_clock::now();
    Graph restartgraph = graph;


    // first solution 
    std::vector<std::vector<int>> first_solution_routes;
    double first_solution_cost = 0;
    int first_solution_trucks_number = create_first_solution(graph); 

    for (int i = 0; i < graph.trucksvector.size(); i++)
        {
            first_solution_cost += graph.trucksvector[i].current_time + graph.distances[graph.trucksvector[i].which_node][0];
        }
    
    for (int l = 0; l < graph.trucksvector.size(); l++)
        {
            first_solution_routes.push_back(graph.trucksvector[l].route);
        }
 
   std::cout<<"starting cost: "<<first_solution_trucks_number+1<<" "<<first_solution_cost<<" \n ";

    //pierwsze rozwiazanie: first_solution_trucks_number, first_solution_cost, first_solution_routes
    //current solution sprawdza czy tworzone neighbour solution jest lepsze od aktualnego(current), best sprawdza w razie resetów przy no improvement

    double current_solution_cost = first_solution_cost;
    int current_solution_used_trucks = first_solution_trucks_number +1;
    std::vector<std::vector<int>> current_solution_routes = first_solution_routes;
    
    double neighbour_solution_cost = current_solution_cost;
    int neighbour_solution_used_trucks = current_solution_used_trucks;
    std::vector<std::vector<int>> neighbour_solution_routes = current_solution_routes;

    double best_solution_cost = current_solution_cost;
    int best_solution_used_trucks = current_solution_used_trucks;
    std::vector<std::vector<int>> best_solution_routes = current_solution_routes;
    

    int used_ops[5] ={0,0,0,0,0};
    Graph graphcopy = graph;
    int no_improvement = 0;

    //parametry simulated annealing 
    double temperature = 300.0;
    double cooling_rate = 0.998;
    double min_temperature = 0.1; 
    
    
    auto start_time = std::chrono::steady_clock::now();
    auto end_time = start_time + std::chrono::seconds(defaultParametersfortabu.time_limit_in_seconds);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    
    //local search z neighbourhood searchem
    while(std::chrono::steady_clock::now() < end_time){  

        generate_neighbour(neighbour_solution_cost, neighbour_solution_used_trucks, neighbour_solution_routes, used_ops, graph);

        //no improvement = reset
        if(no_improvement >= defaultParametersfortabu.no_improvement_limit){
            //zmien tabu list size albo cos w tym stylu
            Graph new_restart = restartgraph;

            if(current_solution_cost <= best_solution_cost && current_solution_used_trucks <= best_solution_used_trucks){
                best_solution_cost = current_solution_cost;
                best_solution_used_trucks = current_solution_used_trucks;
                best_solution_routes = current_solution_routes;
            }

            restart_solution(current_solution_routes, current_solution_cost,current_solution_used_trucks,new_restart);
            
            graph = new_restart;
            graphcopy = new_restart;

            neighbour_solution_cost = current_solution_cost;
            neighbour_solution_used_trucks = current_solution_used_trucks;
            neighbour_solution_routes = current_solution_routes;

            no_improvement = 0;
        }

    //delta do simulated annealing
    double delta_cost = neighbour_solution_cost - current_solution_cost;

    //simulated annealing warunek przyciecia rozwiazania
    if (accept_worse_solution(delta_cost, defaultParametersfortabu.temperature)) {
        current_solution_cost = neighbour_solution_cost;
        current_solution_routes = neighbour_solution_routes;
        current_solution_used_trucks = neighbour_solution_used_trucks;
        //std::cout << "iter: " << iteration << " -> accepted with delta: " << delta_cost << " at temp: " << temperature << "\n";
        update_graph_from_routes(graph, current_solution_routes);
        graphcopy.trucksvector = graph.trucksvector;
        no_improvement = 0;
/*
        if (current_solution_cost > neighbour_solution_cost && current_solution_used_trucks >= neighbour_solution_used_trucks){
            
            current_solution_cost = neighbour_solution_cost;
            current_solution_routes = neighbour_solution_routes;
            current_solution_used_trucks = neighbour_solution_used_trucks;
            std::cout<<"iter: "<<iteration<<" better -> "<<current_solution_used_trucks<<" "<<std::fixed << std::setprecision(5) <<current_solution_cost<<" used_operation "<<used_ops[4]<<" ikjl values: "<<used_ops[0]<<" "<<used_ops[1]<<" "<<used_ops[2]<<" "<<used_ops[3]<<" \n";
            update_graph_from_routes(graph, current_solution_routes);
            graphcopy.trucksvector = graph.trucksvector;
            no_improvement =0;
*/
        }
        else{
            //std::cout<<"iter: "<<iteration<<" worse -> "<<current_solution_used_trucks<<" "<<std::fixed << std::setprecision(5) <<neighbour_solution_cost<<" used_operation "<<used_ops[4]<<" ikjl values: "<<used_ops[0]<<" "<<used_ops[1]<<" "<<used_ops[2]<<" "<<used_ops[3]<<" \n";
            no_improvement ++;
            neighbour_solution_cost = current_solution_cost;
            neighbour_solution_routes = current_solution_routes;
            neighbour_solution_used_trucks = current_solution_used_trucks;
            graph.trucksvector = graphcopy.trucksvector;

        }

        defaultParametersfortabu.temperature = std::max(defaultParametersfortabu.temperature * defaultParametersfortabu.cooling_factor, defaultParametersfortabu.min_temperature);
        
        //aktualizacja sciezek w truckvector
        for (int i = 0; i < current_solution_routes.size(); i++) {
            graph.trucksvector[i].route = current_solution_routes[i];
        }
        
        iteration ++;
        
    } 

    //przy koncu czasu
    if(current_solution_cost < best_solution_cost || current_solution_used_trucks < best_solution_used_trucks){
        best_solution_cost = current_solution_cost;
        best_solution_used_trucks = current_solution_used_trucks;
        best_solution_routes = current_solution_routes;
    }




    //std::cout<<"ending iter: "<<iteration<<"\n";
    auto stop=std::chrono::high_resolution_clock::now();  // Koniec pomiaru
    auto duration=std::chrono::duration_cast<std::chrono::seconds>(stop-start);
    //std::cout<<"\nCzas wykonania funkcji: "<<duration.count()<<" sekund"<<std::endl;

    std::ofstream outputFile("wyniki.txt");
    if (!outputFile.is_open())
    {
        return;
    }

    outputFile << std::fixed << std::setprecision(5) << best_solution_used_trucks << " " << best_solution_cost << std::endl;
    for (int l = 0; l < best_solution_routes.size(); l++)
    {
        if (!best_solution_routes[l].empty())
        {
            for (int k = 0; k < best_solution_routes[l].size(); k++)
            {
                outputFile << best_solution_routes[l][k] << " ";
            }
            outputFile << std::endl;
        }
    }


    outputFile.close();
}


