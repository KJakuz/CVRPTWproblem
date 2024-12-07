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
#include "trucks.h"
#include "parameters.h"
#include "tabu.h"


void Tabu::generate_neighbour(double& current_cost,int& current_used_trucks, std::vector<std::vector<int>>& current_routes, int* used_ops, Graph& graph) {

    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<double> operation_weights = {100, 0, 0, 0, 0, 0, 0};
    std::discrete_distribution<> dist(operation_weights.begin(), operation_weights.end());

    int which_operation = dist(gen);

    switch(which_operation) {

        case 0:
            adaptive_large_neighborhood_search(current_cost, current_used_trucks, current_routes, used_ops, graph);
            break;
        case 1:
            route_splitting_and_merging(current_cost, current_used_trucks, current_routes, used_ops, graph);
            break;
        case 2:
            swap_two_nodes(current_cost, current_used_trucks, current_routes, used_ops, graph);
            break;
        case 3:
            two_opt_swap(current_cost, current_used_trucks, current_routes, used_ops, graph);
            break;
        case 4:
            advanced_node_redistribution(current_cost, current_used_trucks, current_routes, used_ops, graph);
            break;
        case 5:
            relocate_customers(current_cost, current_used_trucks, current_routes, used_ops, graph);
            break;
        case 6:
            cross_route_balanced_swap(current_cost, current_used_trucks, current_routes, used_ops, graph);
            break;
    }

    current_cost = 0;
    current_used_trucks = 0;
    for (int i = 0; i < graph.trucksvector.size(); i++) {
        if (!graph.trucksvector[i].route.empty()) {
            current_cost += graph.trucksvector[i].current_time + graph.distances[(current_routes[i][current_routes[i].size()-1])][0];
            current_used_trucks++;
        }
    }
}


void Tabu::adaptive_large_neighborhood_search(double& current_cost, int& current_used_trucks,std::vector<std::vector<int>>& current_routes, int* used_ops, Graph& graph) {

    const int num_suboperators = 6;
    std::vector<double> weights(num_suboperators, 1.0);
    std::vector<int> scores(num_suboperators, 0);
    std::vector<int> usages(num_suboperators, 0);

    std::random_device rd;
    std::mt19937 gen(rd());

    double best_cost = current_cost;
    std::vector<std::vector<int>> best_routes = current_routes;

    int iter = 0;
    const int max_iterations = 200;
    const double reaction_factor = 0.7; 
    const double exploration_factor = 0.1;

    while (iter < max_iterations) {
        std::discrete_distribution<> dist(weights.begin(), weights.end());
        int selected_suboperator = dist(gen);

        switch(selected_suboperator) {
            case 0:
                swap_two_nodes(current_cost, current_used_trucks, current_routes, used_ops, graph);
                break;
            case 1:
                two_opt_swap(current_cost, current_used_trucks, current_routes, used_ops, graph);
                break;
            case 2:
                advanced_node_redistribution(current_cost, current_used_trucks, current_routes, used_ops, graph);
                break;
            case 3:
                route_splitting_and_merging(current_cost, current_used_trucks, current_routes, used_ops, graph);
                break;
            case 4:
                relocate_customers(current_cost, current_used_trucks, current_routes, used_ops, graph);
                break;
            case 5:
                cross_route_balanced_swap(current_cost, current_used_trucks, current_routes, used_ops, graph);
                break;


        }

        usages[selected_suboperator]++;

        if (current_cost < best_cost) {
            best_cost = current_cost;
            best_routes = current_routes;
            scores[selected_suboperator] += 10; 
        } else {
            scores[selected_suboperator] += (current_cost < best_cost * 1.02) ? 4 : 2;
        }

        if (iter % 200 == 0) {
            for (int i = 0; i < num_suboperators; ++i) {
                if (usages[i] > 0) {
                    double performance = static_cast<double>(scores[i]) / usages[i];
                    weights[i] = std::pow(weights[i], 0.9) * std::pow(performance, 1.1);
                }
            }

            double total_weight = std::accumulate(weights.begin(), weights.end(), 0.0);
            for (auto& w : weights) {
                w /= total_weight;
            }

            scores.assign(num_suboperators, 0);
            usages.assign(num_suboperators, 0);
        }

        iter++;
    }
}


void Tabu::swap_two_nodes(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, int* used_ops, Graph& graph) {
    for (int i = 0; i < current_routes.size(); ++i) {
        for (int j = 0; j < current_routes.size(); ++j) {
            if (i == j || current_routes[i].empty() || current_routes[j].empty()) continue;
            
            for (int k = 0; k < current_routes[i].size(); ++k) {
                for (int l = 0; l < current_routes[j].size(); ++l) {
                        std::vector<int> test_route1 = current_routes[i];
                        std::vector<int> test_route2 = current_routes[j];
                        std::swap(test_route1[k], test_route2[l]);

                        double time_for_route1 = is_single_route_possible(test_route1, graph);
                        double time_for_route2 = is_single_route_possible(test_route2, graph);

                        if (time_for_route1 != -1 && time_for_route2 != -1) {
                            used_ops[0] = i;
                        used_ops[1] = k;
                        used_ops[2] = j;
                        used_ops[3] = l;
                        used_ops[4] = 0;

                    if (not_in_Tabu(used_ops)) {

                            current_routes[i] = test_route1;
                            current_routes[j] = test_route2;

                            graph.trucksvector[i].route = test_route1;
                            graph.trucksvector[j].route = test_route2;

                            graph.trucksvector[i].current_time = time_for_route1;
                            graph.trucksvector[j].current_time = time_for_route2;

                            add_to_Tabu(used_ops);
                            return;
                        }
                    }
                }
            }
        }
    }
}



void Tabu::two_opt_swap(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, int* used_ops, Graph& graph) {
    for (int i = 0; i < current_routes.size(); i++) {
        for (int j = i + 1; j < current_routes.size(); j++) {
            if (current_routes[i].size() < 2 || current_routes[j].size() < 2) continue;

            for (int k = 1; k < current_routes[i].size(); k++) {
                for (int l = 1; l < current_routes[j].size(); l++) {
                
                        std::vector<int> test_route1(current_routes[i].begin(), current_routes[i].begin() + k);
                        test_route1.insert(test_route1.end(), current_routes[j].begin() + l, current_routes[j].end());

                        std::vector<int> test_route2(current_routes[j].begin(), current_routes[j].begin() + l);
                        test_route2.insert(test_route2.end(), current_routes[i].begin() + k, current_routes[i].end());

                        double time_for_route1 = is_single_route_possible(test_route1, graph);
                        double time_for_route2 = is_single_route_possible(test_route2, graph);

                        if (time_for_route1 != -1 && time_for_route2 != -1) {
                            used_ops[0] = i;
                            used_ops[1] = k;
                            used_ops[2] = j;
                            used_ops[3] = l;
                            used_ops[4] = 1; 

                            if (not_in_Tabu(used_ops)) {
                            current_routes[i] = test_route1;
                            current_routes[j] = test_route2;

                            graph.trucksvector[i].route = test_route1;
                            graph.trucksvector[j].route = test_route2;

                            graph.trucksvector[i].current_time = time_for_route1;
                            graph.trucksvector[j].current_time = time_for_route2;

                            add_to_Tabu(used_ops);
                            return;
                        }
                    }
                }
            }
        }
    }
}




void Tabu::advanced_node_redistribution(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, int* used_ops, Graph& graph) {
    for (int i = 0; i < current_routes.size(); ++i) {
        for (int j = 0; j < current_routes.size(); ++j) {
            if (i == j || current_routes[i].empty() || current_routes[j].empty()) continue;

            for (int k = 0; k < current_routes[i].size(); ++k) {
                for (int insert_pos = 0; insert_pos <= current_routes[j].size(); ++insert_pos) {
        
                        std::vector<int> test_route1 = current_routes[i];
                        std::vector<int> test_route2 = current_routes[j];

                        int node_to_move = test_route1[k];
                        test_route1.erase(test_route1.begin() + k);
                        test_route2.insert(test_route2.begin() + insert_pos, node_to_move);

                        double time_for_route1 = is_single_route_possible(test_route1, graph);
                        double time_for_route2 = is_single_route_possible(test_route2, graph);

                        if (time_for_route1 != -1 && time_for_route2 != -1) {
                           
                             used_ops[0] = i;
                            used_ops[1] = k;
                            used_ops[2] = j;
                            used_ops[3] = insert_pos;
                            used_ops[4] = 4; 

                            if (not_in_Tabu(used_ops)) {
                            current_routes[i] = test_route1;
                            current_routes[j] = test_route2;

                            graph.trucksvector[i].route = test_route1;
                            graph.trucksvector[j].route = test_route2;

                            graph.trucksvector[i].current_time = time_for_route1;
                            graph.trucksvector[j].current_time = time_for_route2;

                          
                            add_to_Tabu(used_ops);
                            return;
                        }
                    }
                }
            }
        }
    }
}



void Tabu::relocate_customers(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, int* used_ops, Graph& graph) {
    for (int source_route = 0; source_route < current_routes.size(); ++source_route) {
        for (int target_route = 0; target_route < current_routes.size(); ++target_route) {
            if (source_route == target_route) continue;

            for (int node_idx = 0; node_idx < current_routes[source_route].size(); ++node_idx) {
                for (int insert_pos = 0; insert_pos <= current_routes[target_route].size(); ++insert_pos) {
                        std::vector<int> test_source_route = current_routes[source_route];
                        std::vector<int> test_target_route = current_routes[target_route];

                        int customer = test_source_route[node_idx];
                        test_source_route.erase(test_source_route.begin() + node_idx);
                        test_target_route.insert(test_target_route.begin() + insert_pos, customer);

                        double source_route_time = is_single_route_possible(test_source_route, graph);
                        double target_route_time = is_single_route_possible(test_target_route, graph);

                        if (source_route_time != -1 && target_route_time != -1) {
                            used_ops[0] = source_route;
                            used_ops[1] = node_idx;
                            used_ops[2] = target_route;
                            used_ops[3] = insert_pos;
                            used_ops[4] = 6; 

                            if (not_in_Tabu(used_ops)){
                                current_routes[source_route] = test_source_route;
                                current_routes[target_route] = test_target_route;

                                graph.trucksvector[source_route].route = test_source_route;
                                graph.trucksvector[target_route].route = test_target_route;

                                graph.trucksvector[source_route].current_time = source_route_time;
                                graph.trucksvector[target_route].current_time = target_route_time;

                                add_to_Tabu(used_ops);
                                return;
                        }
                    }
                }
            }
        }
    }
}

void Tabu::cross_route_balanced_swap(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, int* used_ops, Graph& graph) {
        for (int route1 = 0; route1 < current_routes.size(); ++route1) {
            for (int route2 = route1 + 1; route2 < current_routes.size(); ++route2) {
                for (int node1_idx = 0; node1_idx < current_routes[route1].size(); ++node1_idx) {
                    for (int node2_idx = 0; node2_idx < current_routes[route2].size(); ++node2_idx) {
                        std::vector<int> test_route1 = current_routes[route1];
                        std::vector<int> test_route2 = current_routes[route2];

                        std::swap(test_route1[node1_idx], test_route2[node2_idx]);

                        double route1_time = is_single_route_possible(test_route1, graph);
                        double route2_time = is_single_route_possible(test_route2, graph);

                        if (route1_time != -1 && route2_time != -1) {
                            used_ops[0] = route1;
                            used_ops[1] = node1_idx;
                            used_ops[2] = route2;
                            used_ops[3] = node2_idx;
                            used_ops[4] = 7;
                            if (not_in_Tabu(used_ops)){

                                current_routes[route1] = test_route1;
                                current_routes[route2] = test_route2;

                                graph.trucksvector[route1].route = test_route1;
                                graph.trucksvector[route2].route = test_route2;

                                graph.trucksvector[route1].current_time = route1_time;
                                graph.trucksvector[route2].current_time = route2_time;

                                add_to_Tabu(used_ops);
                                return;
                        }
                    }
                }
            }
        }
    }
}

void Tabu::route_splitting_and_merging(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, int* used_ops, Graph& graph) {
    
    std::vector<int> route_lengths(current_routes.size());
    std::vector<double> route_times(current_routes.size());
    
    for (int i = 0; i < current_routes.size(); ++i) {
        route_lengths[i] = current_routes[i].size();
        route_times[i] = graph.trucksvector[i].current_time;
    }

    std::vector<int> route_indices(current_routes.size());
    std::iota(route_indices.begin(), route_indices.end(), 0);
    
    for (int idx : route_indices) {
        if (current_routes[idx].size() <= 2 || idx < 0 || idx > current_routes.size()) continue;

        int split_point = current_routes[idx].size() / 2;
        
        std::vector<int> first_part(current_routes[idx].begin(), current_routes[idx].begin() + split_point);
        std::vector<int> second_part(current_routes[idx].begin() + split_point, current_routes[idx].end());

        double first_route_time = is_single_route_possible(first_part, graph);
        double second_route_time = is_single_route_possible(second_part, graph);

        if (first_route_time != -1 && second_route_time != -1) {
            for (int merge_idx = 0; merge_idx < current_routes.size(); ++merge_idx) {
                if (merge_idx == idx) continue;

                for (int insert_pos = 0; insert_pos <= current_routes[merge_idx].size(); ++insert_pos) {
                    std::vector<int> test_merged_route = current_routes[merge_idx];
                    test_merged_route.insert(
                        test_merged_route.begin() + insert_pos, 
                        second_part.begin(), 
                        second_part.end()
                    );

                    double merged_route_time = is_single_route_possible(test_merged_route, graph);
                    if (not_in_Tabu(used_ops)) {
                        if (merged_route_time != -1) {

                            used_ops[0] = idx;
                            used_ops[1] = split_point;
                            used_ops[2] = merge_idx;
                            used_ops[3] = insert_pos;
                            used_ops[4] = 5; 

                                current_routes[idx] = first_part;
                                current_routes[merge_idx] = test_merged_route;

                                graph.trucksvector[idx].route = first_part;
                                graph.trucksvector[merge_idx].route = test_merged_route;

                                graph.trucksvector[idx].current_time = first_route_time;
                                graph.trucksvector[merge_idx].current_time = merged_route_time;
                                add_to_Tabu(used_ops);
                                return;
                        }
                    }
                }
            }
        }
    }
}

