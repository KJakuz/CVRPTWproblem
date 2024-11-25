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
        int Max_iterations = 10;

        int create_first_solution_with_grasp(Graph& graph);
        void Tabu_search(Graph& graph);
        double calculate_cost(const std::vector<std::vector<int>>& routes, Graph& graph);
        void update_tabu_list();
        bool is_tabu(int node1, int node2);

};

#endif