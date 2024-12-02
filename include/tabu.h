#ifndef TABU_H
#define TABU_H
#include <vector>
#include "trucks.h"
#include "parameters.h"
#include "graph.h"

extern Parameters default_parameters;

class Node;
class Truck;

class Tabu
{
    public:
        int Tabu_list[default_parameters.Tabu_list_size][5];
        int no_improvement_limit = 1000;
        int Max_iterations = 20;

        int create_first_solution_with_grasp(Graph& graph);
        void Tabu_search(Graph& graph);
        double calculate_cost(const std::vector<std::vector<int>>& routes);
        void generate_neighbour(double& current_cost,int& current_used_trucks,std::vector<std::vector<int>>& current_routes,int* used_ops, Graph& graph);
        bool can_be_swaped(int index_of_node1, int index_of_node2,int idx_of_truck1,int idx_of_truck2,Graph& graph,std::vector<std::vector<int>>& current_routes);
        double is_single_route_possible(std::vector<int>& route,Graph& graph);
        bool not_in_Tabu(const int operation[5]);
        void add_to_Tabu(const int operation[5],int iteration);
        void swap_two_nodes(double& current_cost,int& current_used_trucks,std::vector<std::vector<int>>& current_routes,int* used_ops, Graph& graph);
        void two_opt_swap(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, int* used_ops, Graph& graph);

};

#endif