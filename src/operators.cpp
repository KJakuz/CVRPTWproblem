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

extern double temperature;

void Tabu::generate_neighbour(double& current_cost,int& current_used_trucks, std::vector<std::vector<int>>& current_routes, int* used_ops, Graph& graph) {

    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<double> operation_weights = { 3, 3, 5, 4, 3, 25, 12};
    //std::vector<double> operation_weights = { 1,1,1,1,1,1,1};

    std::discrete_distribution<> dist(operation_weights.begin(), operation_weights.end());

    int which_operation = dist(gen);

    switch(which_operation) {

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
            relocate_customers(current_cost, current_used_trucks, current_routes, used_ops, graph);
            break;
        case 4:
            cross_route_balanced_swap(current_cost, current_used_trucks, current_routes, used_ops, graph);
            break;
        case 5:
            route_splitting_and_merging(current_cost, current_used_trucks, current_routes, used_ops, graph);
            break;
        case 6:
            reduce_truck_count(current_cost, current_used_trucks, current_routes, used_ops, graph);
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



void Tabu::swap_two_nodes(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, int* used_ops, Graph& graph) {
    double best_absolute_delta = 0; 
    double best_cost = current_cost;
    int best_i = -1, best_j = -1, best_k = -1, best_l = -1;
    std::vector<std::vector<int>> best_routes = current_routes;
    double best_time_for_route1 = -1, best_time_for_route2 = -1, best_delta;

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
                            double old_time_for_route1 = graph.trucksvector[i].current_time;
                            double old_time_for_route2 = graph.trucksvector[j].current_time;

                            double new_cost = current_cost - old_time_for_route1 - old_time_for_route2 
                                              + time_for_route1 + time_for_route2;

                            double delta_cost = new_cost - current_cost;
                            double absolute_delta = std::abs(delta_cost); 
                            if (delta_cost > 0 && !accept_worse_solution(delta_cost, temperature)) continue;
                            if (absolute_delta > best_absolute_delta) {
                                best_absolute_delta = absolute_delta;
                                best_cost = new_cost;
                                best_i = i;
                                best_j = j;
                                best_k = k;
                                best_l = l;
                                best_delta = delta_cost;
                                best_routes = {test_route1, test_route2};
                                best_time_for_route1 = time_for_route1;
                                best_time_for_route2 = time_for_route2;
                            }
                        }
                    }
                }
            }
        }
    }

    if (best_i != -1 && best_j != -1) {
        current_routes[best_i] = best_routes[0];
        current_routes[best_j] = best_routes[1];

        graph.trucksvector[best_i].route = best_routes[0];
        graph.trucksvector[best_j].route = best_routes[1];

        graph.trucksvector[best_i].current_time = best_time_for_route1;
        graph.trucksvector[best_j].current_time = best_time_for_route2;

        current_cost = best_cost;

        used_ops[0] = best_i;
        used_ops[1] = best_k;
        used_ops[2] = best_j;
        used_ops[3] = best_l;
        used_ops[4] = 0;
        used_ops[5] = best_delta;
        add_to_Tabu(used_ops);
    }
}



void Tabu::two_opt_swap(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, int* used_ops, Graph& graph) {
    double best_absolute_delta = 0;
    double best_cost = current_cost;
    int best_i = -1, best_j = -1, best_k = -1, best_l = -1;
    std::vector<std::vector<int>> best_routes = current_routes;
    double best_time_for_route1 = -1, best_time_for_route2 = -1;

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
                            double old_time_for_route1 = graph.trucksvector[i].current_time;
                            double old_time_for_route2 = graph.trucksvector[j].current_time;

                            double new_cost = current_cost - old_time_for_route1 - old_time_for_route2
                                              + time_for_route1 + time_for_route2;
                            double delta_cost = new_cost - current_cost;

                            if (delta_cost < 0 || accept_worse_solution(delta_cost, temperature)) {
                                if (std::abs(delta_cost) > best_absolute_delta) {
                                    best_absolute_delta = std::abs(delta_cost);
                                    best_cost = new_cost;
                                    best_i = i;
                                    best_j = j;
                                    best_k = k;
                                    best_l = l;
                                    best_routes = {test_route1, test_route2};
                                    best_time_for_route1 = time_for_route1;
                                    best_time_for_route2 = time_for_route2;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (best_i != -1 && best_j != -1) {
        current_routes[best_i] = best_routes[0];
        current_routes[best_j] = best_routes[1];

        graph.trucksvector[best_i].route = best_routes[0];
        graph.trucksvector[best_j].route = best_routes[1];

        graph.trucksvector[best_i].current_time = best_time_for_route1;
        graph.trucksvector[best_j].current_time = best_time_for_route2;

        current_cost = best_cost;

        used_ops[0] = best_i;
        used_ops[1] = best_k;
        used_ops[2] = best_j;
        used_ops[3] = best_l;
        used_ops[4] = 1;

        add_to_Tabu(used_ops);
    }
}



void Tabu::advanced_node_redistribution(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, int* used_ops, Graph& graph) {
    double best_absolute_delta = 0;
    double best_cost = current_cost;
    int best_i = -1, best_j = -1, best_k = -1, best_insert_pos = -1;
    std::vector<std::vector<int>> best_routes = current_routes;
    double best_time_for_route1 = -1, best_time_for_route2 = -1;

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
                        used_ops[4] = 2;

                        if (not_in_Tabu(used_ops)) {
                            double old_time_for_route1 = graph.trucksvector[i].current_time;
                            double old_time_for_route2 = graph.trucksvector[j].current_time;

                            double new_cost = current_cost - old_time_for_route1 - old_time_for_route2
                                              + time_for_route1 + time_for_route2;
                            double delta_cost = new_cost - current_cost;

                            if (delta_cost < 0 || accept_worse_solution(delta_cost, temperature)) {
                                if (std::abs(delta_cost) > best_absolute_delta) {
                                    best_absolute_delta = std::abs(delta_cost);
                                    best_cost = new_cost;
                                    best_i = i;
                                    best_j = j;
                                    best_k = k;
                                    best_insert_pos = insert_pos;
                                    best_routes = {test_route1, test_route2};
                                    best_time_for_route1 = time_for_route1;
                                    best_time_for_route2 = time_for_route2;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (best_i != -1 && best_j != -1) {
        current_routes[best_i] = best_routes[0];
        current_routes[best_j] = best_routes[1];

        graph.trucksvector[best_i].route = best_routes[0];
        graph.trucksvector[best_j].route = best_routes[1];

        graph.trucksvector[best_i].current_time = best_time_for_route1;
        graph.trucksvector[best_j].current_time = best_time_for_route2;

        current_cost = best_cost;

        used_ops[0] = best_i;
        used_ops[1] = best_k;
        used_ops[2] = best_j;
        used_ops[3] = best_insert_pos;
        used_ops[4] = 2;

        add_to_Tabu(used_ops);
    }
}


void Tabu::relocate_customers(double& current_cost,int& current_used_trucks,std::vector<std::vector<int>>& current_routes,int* used_ops,Graph& graph){
    double best_cost=current_cost;
    int best_source_route=-1,best_target_route=-1,best_node_idx=-1,best_insert_pos=-1;
    std::vector<std::vector<int>> best_routes=current_routes;
    double best_source_route_time=-1,best_target_route_time=-1;

    for(int source_route=0;source_route<current_routes.size();++source_route){
        for(int target_route=0;target_route<current_routes.size();++target_route){
            if(source_route==target_route) continue;

            for(int node_idx=0;node_idx<current_routes[source_route].size();++node_idx){
                for(int insert_pos=0;insert_pos<=current_routes[target_route].size();++insert_pos){
                    std::vector<int> test_source_route=current_routes[source_route];
                    std::vector<int> test_target_route=current_routes[target_route];

                    int customer=test_source_route[node_idx];
                    test_source_route.erase(test_source_route.begin()+node_idx);
                    test_target_route.insert(test_target_route.begin()+insert_pos,customer);

                    double source_route_time=is_single_route_possible(test_source_route,graph);
                    double target_route_time=is_single_route_possible(test_target_route,graph);

                    if(source_route_time!=-1&&target_route_time!=-1){
                        used_ops[0]=source_route;
                        used_ops[1]=node_idx;
                        used_ops[2]=target_route;
                        used_ops[3]=insert_pos;
                        used_ops[4]=3;

                        if(not_in_Tabu(used_ops)){
                            double old_time_for_source_route=graph.trucksvector[source_route].current_time;
                            double old_time_for_target_route=graph.trucksvector[target_route].current_time;

                            double new_cost=current_cost-old_time_for_source_route-old_time_for_target_route+source_route_time+target_route_time;
                            double delta_cost=new_cost-current_cost;
                            if(delta_cost<0||accept_worse_solution(delta_cost,temperature)){
                                if(new_cost<best_cost){
                                    best_cost=new_cost;
                                    best_source_route=source_route;
                                    best_target_route=target_route;
                                    best_node_idx=node_idx;
                                    best_insert_pos=insert_pos;
                                    best_routes={test_source_route,test_target_route};
                                    best_source_route_time=source_route_time;
                                    best_target_route_time=target_route_time;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if(best_source_route!=-1&&best_target_route!=-1){
        current_routes[best_source_route]=best_routes[0];
        current_routes[best_target_route]=best_routes[1];

        graph.trucksvector[best_source_route].route=best_routes[0];
        graph.trucksvector[best_target_route].route=best_routes[1];

        graph.trucksvector[best_source_route].current_time=best_source_route_time;
        graph.trucksvector[best_target_route].current_time=best_target_route_time;

        current_cost=best_cost;
        used_ops[0]=best_source_route;
        used_ops[1]=best_node_idx;
        used_ops[2]=best_target_route;
        used_ops[3]=best_insert_pos;
        used_ops[4]=3;
        add_to_Tabu(used_ops);
    }
}


void Tabu::cross_route_balanced_swap(double& current_cost,int& current_used_trucks,std::vector<std::vector<int>>& current_routes,int* used_ops,Graph& graph){
    double best_cost=current_cost;
    int best_route1=-1,best_route2=-1;
    int best_node1_idx=-1,best_node2_idx=-1;
    std::vector<int> best_route1_vec,best_route2_vec;
    double best_route1_time=-1,best_route2_time=-1;

    for(int route1=0;route1<current_routes.size();++route1){
        if(current_routes[route1].empty()) continue;

        for(int route2=route1+1;route2<current_routes.size();++route2){
            if(current_routes[route2].empty()) continue;

            for(int node1_idx=0;node1_idx<current_routes[route1].size();++node1_idx){
                for(int node2_idx=0;node2_idx<current_routes[route2].size();++node2_idx){
                    std::vector<int> test_route1=current_routes[route1];
                    std::vector<int> test_route2=current_routes[route2];
                    std::swap(test_route1[node1_idx],test_route2[node2_idx]);

                    double route1_time=is_single_route_possible(test_route1,graph);
                    double route2_time=is_single_route_possible(test_route2,graph);

                    if(route1_time!=-1&&route2_time!=-1){
                        used_ops[0]=route1;
                        used_ops[1]=node1_idx;
                        used_ops[2]=route2;
                        used_ops[3]=node2_idx;
                        used_ops[4]=4;

                        if(not_in_Tabu(used_ops)){
                            double old_route1_time=graph.trucksvector[route1].current_time;
                            double old_route2_time=graph.trucksvector[route2].current_time;

                            double new_cost=current_cost-old_route1_time-old_route2_time+route1_time+route2_time;
                            double delta_cost=new_cost-current_cost;
                            if(delta_cost<0||accept_worse_solution(delta_cost,temperature)){
                                if(new_cost<best_cost){
                                    best_cost=new_cost;
                                    best_route1=route1;
                                    best_route2=route2;
                                    best_node1_idx=node1_idx;
                                    best_node2_idx=node2_idx;
                                    best_route1_vec=test_route1;
                                    best_route2_vec=test_route2;
                                    best_route1_time=route1_time;
                                    best_route2_time=route2_time;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if(best_route1!=-1&&best_route2!=-1){
        current_routes[best_route1]=best_route1_vec;
        current_routes[best_route2]=best_route2_vec;

        graph.trucksvector[best_route1].route=best_route1_vec;
        graph.trucksvector[best_route2].route=best_route2_vec;

        graph.trucksvector[best_route1].current_time=best_route1_time;
        graph.trucksvector[best_route2].current_time=best_route2_time;

        current_cost=best_cost;

        used_ops[0]=best_route1;
        used_ops[1]=best_node1_idx;
        used_ops[2]=best_route2;
        used_ops[3]=best_node2_idx;
        used_ops[4]=4;

        add_to_Tabu(used_ops);
    }
}


void Tabu::route_splitting_and_merging(double& current_cost,int& current_used_trucks,std::vector<std::vector<int>>& current_routes,int* used_ops,Graph& graph) {
    for (int idx=0;idx<current_routes.size();++idx) {
        if (current_routes[idx].size()<=2) continue;

        int split_point= rand() % (current_routes[idx].size() - 1) + 1;


        std::vector<int> first_part(current_routes[idx].begin(),current_routes[idx].begin()+split_point);
        std::vector<int> second_part(current_routes[idx].begin()+split_point,current_routes[idx].end());

        double first_route_time=is_single_route_possible(first_part,graph);
        double second_route_time=is_single_route_possible(second_part,graph);

        if (first_route_time==-1||second_route_time==-1) continue;

        bool improvement_found=false;
        for (int merge_idx=0;merge_idx<current_routes.size();++merge_idx) {
            if (merge_idx==idx) continue;

            for (int insert_pos=0;insert_pos<=current_routes[merge_idx].size();++insert_pos) {
                std::vector<int> test_merged_route=current_routes[merge_idx];
                test_merged_route.insert(test_merged_route.begin()+insert_pos,second_part.begin(),second_part.end());

                double merged_route_time=is_single_route_possible(test_merged_route,graph);

                if (merged_route_time!=-1) {
                    used_ops[0]=idx;
                    used_ops[1]=split_point;
                    used_ops[2]=merge_idx;
                    used_ops[3]=insert_pos;
                    used_ops[4]=5;

                    if (not_in_Tabu(used_ops)) {
                        double old_merge_time=graph.trucksvector[merge_idx].current_time;
                        double old_idx_time=graph.trucksvector[idx].current_time;
                        double new_cost=current_cost-old_merge_time-old_idx_time+first_route_time+merged_route_time;

                        double delta_cost=new_cost-current_cost;
                        if (delta_cost<0||accept_worse_solution(delta_cost,temperature)) {
                            if (new_cost<current_cost) {
                                current_cost=new_cost;

                                current_routes[idx]=first_part;
                                current_routes[merge_idx]=test_merged_route;

                                graph.trucksvector[idx].route=first_part;
                                graph.trucksvector[merge_idx].route=test_merged_route;

                                graph.trucksvector[idx].current_time=first_route_time;
                                graph.trucksvector[merge_idx].current_time=merged_route_time;

                                add_to_Tabu(used_ops);

                                improvement_found=true;
                                break;
                            }
                        }
                    }
                }
            }
            if (improvement_found) break;
        }
    }
}


void Tabu::reduce_truck_count(double& current_cost, int& current_used_trucks, std::vector<std::vector<int>>& current_routes, int* used_ops, Graph& graph) {
    double best_cost = current_cost;
    int best_source_route = -1, best_target_route = -1, best_insert_pos =-1;
    std::vector<std::vector<int>> best_routes = current_routes;
    double best_target_route_time = -1, best_source_route_time = -1;
    bool truck_reduction_possible = false;

    std::vector<std::pair<int, double>> route_utilization;
    for (int i = 0; i < current_routes.size(); ++i) {
        if (!current_routes[i].empty()) {
            route_utilization.push_back({i, graph.trucksvector[i].current_time});
        }
    }
    std::sort(route_utilization.begin(), route_utilization.end(), 
              [](const auto& a, const auto& b) { return a.second < b.second; });

    for (auto& [source_route_idx, _] : route_utilization) {
        if (current_routes[source_route_idx].empty()) continue;

        for (int target_route = 0; target_route < current_routes.size(); ++target_route) {
            if (target_route == source_route_idx || current_routes[target_route].empty()) continue;

            for (int insert_pos = 0; insert_pos <= current_routes[target_route].size(); ++insert_pos) {
                std::vector<int> test_source_route = current_routes[source_route_idx];
                std::vector<int> test_target_route = current_routes[target_route];

                test_target_route.insert(test_target_route.begin() + insert_pos, 
                                         test_source_route.begin(), 
                                         test_source_route.end());
                test_source_route.clear();

                double target_route_time = is_single_route_possible(test_target_route, graph);
                if (target_route_time != -1) {
                    used_ops[0] = source_route_idx;
                    used_ops[1] = 0; 
                    used_ops[2] = target_route;
                    used_ops[3] = insert_pos;
                    used_ops[4] = 6;
                    if (not_in_Tabu(used_ops)) {
                        truck_reduction_possible = true;

                        double old_time_for_source_route = graph.trucksvector[source_route_idx].current_time;
                        double old_time_for_target_route = graph.trucksvector[target_route].current_time;

                        double new_cost = current_cost - old_time_for_source_route - old_time_for_target_route 
                                          + target_route_time;
double delta_cost = new_cost - current_cost;  
                            if (delta_cost < 0 || accept_worse_solution(delta_cost, temperature)) {
                        if (new_cost < best_cost) {
                            best_cost = new_cost;
                            best_source_route = source_route_idx;
                            best_target_route = target_route;
                            best_insert_pos = insert_pos;
                            best_routes = current_routes;
                            best_routes[source_route_idx].clear(); 
                            best_routes[target_route] = test_target_route;
                            best_target_route_time = target_route_time;
                        }
                    }}
                }
            }
        }
    }

    if (truck_reduction_possible && best_source_route != -1 && best_target_route != -1) {
        current_routes = best_routes;

        graph.trucksvector[best_source_route].route.clear();
        graph.trucksvector[best_target_route].route = best_routes[best_target_route];

        graph.trucksvector[best_target_route].current_time = best_target_route_time;

        current_cost = best_cost;
        current_used_trucks--;

        used_ops[0] = best_source_route;
        used_ops[1] = static_cast<int>(best_target_route_time);
        used_ops[2] = best_target_route;
        used_ops[3] = best_insert_pos;
        used_ops[4] = 6;
        add_to_Tabu(used_ops);
    }
}


