#ifndef TABU_H
#define TABU_H
#include <vector>
#include "trucks.h"
#include "parameters.h"
#include "graph.h"

class Parameters;
class Node;
class Truck;

class Tabu
{
    public:
        int Tabu_list_size = 50;
        std::vector<std::pair<int, int>> tabu_list;
        int no_improvement_limit = 5;
        int Max_iterations = 1;

        int create_first_solution_with_grasp(Graph& graph);
        void Tabu_search(Graph& graph);
        double calculate_cost(const std::vector<std::vector<int>>& routes);
        void generate_neighbour(double& current_cost,int& current_used_trucks,std::vector<std::vector<int>>& current_routes, Graph& graph);
        bool can_be_swaped(int index_of_node1, int index_of_node2,int idx_of_truck1,int idx_of_truck2,Graph& graph);
        double is_single_route_possible(std::vector<int>& route,Graph& graph);


};

#endif