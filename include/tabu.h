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
        double Tabu_list[default_parameters.Tabu_list_size][6] = {0.0};

        //main functions
        int create_first_solution(Graph& graph);
        void Tabu_search(Graph& graph);

        double is_single_route_possible(std::vector<int>& route,Graph& graph);
        void restart_solution(std::vector<std::vector<int>>& first_solution_routes, double& first_solution_cost, int& first_solution_trucks_number, Graph& graph);

        //tabu list 
        bool not_in_Tabu(const double operation[6]);
        void add_to_Tabu(const double operation[6]);
        
        //creating neighbourhood
        void generate_neighbour(double& current_cost,int& current_used_trucks,std::vector<std::vector<int>>& current_routes,double* used_ops, Graph& graph);
        void adaptive_large_neighborhood_search(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, double* used_ops, Graph& graph);
        
        bool accept_worse_solution(double delta_cost, double temperature);


        //operators
        void swap_two_nodes(double& current_cost,int& current_used_trucks,std::vector<std::vector<int>>& current_routes,double* used_ops, Graph& graph);
        void two_opt_swap(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, double* used_ops, Graph& graph);
        void advanced_node_redistribution(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, double* used_ops, Graph& graph);
        void route_splitting_and_merging(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, double* used_ops, Graph& graph);
        void relocate_customers(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, double* used_ops, Graph& graph);
        void cross_route_balanced_swap(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, double* used_ops, Graph& graph);
        void cluster_swap(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, double* used_ops, Graph& graph);
        void global_swap(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, double* used_ops, Graph& graph);
        void reduce_truck_count(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, double* used_ops, Graph& graph);
        void route_reconstruction(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, double* used_ops, Graph& graph);
        void aggressive_route_consolidation(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, double* used_ops, Graph& graph);
                                     
};

#endif