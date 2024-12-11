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
#include "simulated_annealing.h"

extern double temperature;
extern double iteration;

void Sa::generate_neighbour(double &current_cost, int &current_used_trucks, std::vector<std::vector<int>> &current_routes, double *used_ops, Graph &graph)
{

    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<double> operation_weights = {1, 1, 1, 1, 1, 1, 3};

    std::discrete_distribution<> dist(operation_weights.begin(), operation_weights.end());

    int which_operation = dist(gen);

    switch (which_operation)
    {

    case 0:
        swap_two_nodes(current_cost, current_used_trucks, current_routes, used_ops, graph);
        break;
    case 1:
        multi_route_two_opt_swap(current_cost, current_used_trucks, current_routes, used_ops, graph);
        break;
    case 2:
        move_node(current_cost, current_used_trucks, current_routes, used_ops, graph);
        break;
    case 3:
        route_splitting_and_merging(current_cost, current_used_trucks, current_routes, used_ops, graph);
        break;
    case 4:
        reduce_truck_count(current_cost, current_used_trucks, current_routes, used_ops, graph);
        break;
    }

    current_cost = 0;
    current_used_trucks = 0;
    for (int i = 0; i < graph.trucksvector.size(); i++)
    {
        if (!graph.trucksvector[i].route.empty())
        {
            current_cost += graph.trucksvector[i].current_time + graph.distances[(current_routes[i][current_routes[i].size() - 1])][0];
            current_used_trucks++;
        }
    }
}

void Sa::swap_two_nodes(double &current_cost, int &current_used_trucks, std::vector<std::vector<int>> &current_routes, double *used_ops, Graph &graph)
{
    double best_absolute_delta = 0;
    double best_cost = current_cost;
    std::vector<std::vector<int>> best_routes = current_routes;

    double best_time_for_route1 = -1, best_time_for_route2 = -1, best_delta = -1;
    int best_i = -1, best_j = -1, best_k = -1, best_l = -1;

    for (int i = 0; i < current_routes.size(); ++i)
    {
        for (int j = 0; j < current_routes.size(); ++j)
        {
            if (i == j || current_routes[i].empty() || current_routes[j].empty())
                continue;
            for (int k = 0; k < current_routes[i].size(); ++k)
            {
                for (int l = 0; l < current_routes[j].size(); ++l)
                {

                    std::vector<int> test_route1 = current_routes[i];
                    std::vector<int> test_route2 = current_routes[j];
                    std::swap(test_route1[k], test_route2[l]);

                    double time_for_route1 = is_single_route_possible(test_route1, graph);
                    double time_for_route2 = is_single_route_possible(test_route2, graph);

                    if (time_for_route1 != -1 && time_for_route2 != -1)
                    {
                        used_ops[0] = 0; // id operacji swap

                        double new_cost = current_cost - graph.trucksvector[i].current_time - graph.trucksvector[j].current_time - graph.distances[0][current_routes[i].back()] - graph.distances[0][current_routes[j].back()] + time_for_route1 + time_for_route2 + graph.distances[0][test_route1.back()] + graph.distances[0][test_route2.back()];
                        double delta_cost = new_cost - current_cost;

                        used_ops[1] = delta_cost;

                        if (not_in_Tabu(used_ops))
                        {

                            double absolute_delta = std::abs(delta_cost);
                            if (delta_cost > 0 && !accept_worse_solution(delta_cost, temperature))
                                continue;

                            if (absolute_delta > best_absolute_delta && static_cast<int>(iteration) % 2 == 0)
                            {
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
                            else if (new_cost > best_cost && static_cast<int>(iteration) % 2 == 1)
                            {
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

    if (best_i != -1 && best_j != -1)
    {
        current_routes[best_i] = best_routes[0];
        current_routes[best_j] = best_routes[1];

        graph.trucksvector[best_i].route = best_routes[0];
        graph.trucksvector[best_j].route = best_routes[1];

        graph.trucksvector[best_i].current_time = best_time_for_route1;
        graph.trucksvector[best_j].current_time = best_time_for_route2;

        if (current_cost != best_cost)
        {
            current_cost = best_cost;

            used_ops[0] = 0;
            used_ops[1] = std::fabs(best_delta);
            add_to_Tabu(used_ops);
        }
    }
}

void Sa::multi_route_two_opt_swap(double &current_cost, int &current_used_trucks, std::vector<std::vector<int>> &current_routes, double *used_ops, Graph &graph)
{
    double best_absolute_delta = 0;
    double best_cost = current_cost;
    std::vector<std::vector<int>> best_routes = current_routes;

    double best_time_for_route1 = -1, best_time_for_route2 = -1, best_delta = -1;
    int best_i = -1, best_j = -1, best_k = -1, best_l = -1;

    for (int i = 0; i < current_routes.size(); i++)
    {
        for (int j = i + 1; j < current_routes.size(); j++)
        {
            if (current_routes[i].size() < 2 || current_routes[j].size() < 2)
                continue;

            for (int k = 1; k < current_routes[i].size(); k++)
            {
                for (int l = 1; l < current_routes[j].size(); l++)
                {
                    // laczenie prefixu jednej trasy z sufixem drugiej
                    std::vector<int> test_route1(current_routes[i].begin(), current_routes[i].begin() + k);
                    test_route1.insert(test_route1.end(), current_routes[j].begin() + l, current_routes[j].end());

                    std::vector<int> test_route2(current_routes[j].begin(), current_routes[j].begin() + l);
                    test_route2.insert(test_route2.end(), current_routes[i].begin() + k, current_routes[i].end());

                    double time_for_route1 = is_single_route_possible(test_route1, graph);
                    double time_for_route2 = is_single_route_possible(test_route2, graph);

                    if (time_for_route1 != -1 && time_for_route2 != -1)
                    {

                        used_ops[0] = 1;

                        double new_cost = current_cost - graph.trucksvector[i].current_time - graph.trucksvector[j].current_time - graph.distances[0][current_routes[i].back()] - graph.distances[0][current_routes[j].back()] + time_for_route1 + time_for_route2 + graph.distances[0][test_route1.back()] + graph.distances[0][test_route2.back()];
                        double delta_cost = new_cost - current_cost;

                        used_ops[1] = std::fabs(delta_cost);

                        if (not_in_Tabu(used_ops))
                        {
                            if (delta_cost < 0 || accept_worse_solution(delta_cost, temperature))
                            {
                                if (std::abs(delta_cost) > best_absolute_delta)
                                {
                                    best_absolute_delta = std::abs(delta_cost);
                                    best_cost = new_cost;
                                    best_i = i;
                                    best_j = j;
                                    best_k = k;
                                    best_l = l;
                                    best_delta = std::fabs(delta_cost);
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

    if (best_i != -1 && best_j != -1)
    {
        current_routes[best_i] = best_routes[0];
        current_routes[best_j] = best_routes[1];

        graph.trucksvector[best_i].route = best_routes[0];
        graph.trucksvector[best_j].route = best_routes[1];

        graph.trucksvector[best_i].current_time = best_time_for_route1;
        graph.trucksvector[best_j].current_time = best_time_for_route2;

        current_cost = best_cost;

        used_ops[0] = 1;
        used_ops[1] = std::fabs(best_delta);

        add_to_Tabu(used_ops);
    }
}

void Sa::move_node(double &current_cost, int &current_used_trucks, std::vector<std::vector<int>> &current_routes, double *used_ops, Graph &graph)
{
    double best_absolute_delta = 0;
    double best_cost = current_cost;
    std::vector<std::vector<int>> best_routes = current_routes;
    double best_time_for_route1 = -1, best_time_for_route2 = -1, best_delta = -1;
    int best_i = -1, best_j = -1, best_from = -1, best_insert_pos = -1;

    for (int i = 0; i < current_routes.size(); ++i)
    {
        for (int j = 0; j < current_routes.size(); ++j)
        {
            if (i == j || current_routes[i].empty() || current_routes[j].empty())
                continue;

            for (int from = 0; from < current_routes[i].size(); ++from)
            {
                for (int insert_pos = 0; insert_pos <= current_routes[j].size(); ++insert_pos)
                {
                    std::vector<int> test_route1 = current_routes[i];
                    std::vector<int> test_route2 = current_routes[j];

                    int node_to_move = test_route1[from];
                    test_route1.erase(test_route1.begin() + from);
                    test_route2.insert(test_route2.begin() + insert_pos, node_to_move);

                    double time_for_route1 = is_single_route_possible(test_route1, graph);
                    double time_for_route2 = is_single_route_possible(test_route2, graph);

                    if (time_for_route1 != -1 && time_for_route2 != -1)
                    {
                        used_ops[0] = 2;

                        double new_cost = current_cost - graph.trucksvector[i].current_time - graph.trucksvector[j].current_time - graph.distances[0][current_routes[i].back()] - graph.distances[0][current_routes[j].back()] + time_for_route1 + time_for_route2 + graph.distances[0][test_route1.back()] + graph.distances[0][test_route2.back()];
                        double delta_cost = new_cost - current_cost;

                        used_ops[1] = delta_cost;
                        if (not_in_Tabu(used_ops))
                        {

                            if (delta_cost < 0 || accept_worse_solution(delta_cost, temperature))
                            {
                                if (std::fabs(delta_cost) > best_absolute_delta && static_cast<int>(iteration) % 2 == 0)
                                {
                                    best_absolute_delta = std::abs(delta_cost);
                                    best_cost = new_cost;
                                    best_i = i;
                                    best_j = j;
                                    best_from = from;
                                    best_insert_pos = insert_pos;
                                    best_delta = delta_cost;
                                    best_routes = {test_route1, test_route2};
                                    best_time_for_route1 = time_for_route1;
                                    best_time_for_route2 = time_for_route2;
                                }
                                else if (new_cost < best_cost && static_cast<int>(iteration) % 2 == 1)
                                {
                                    best_absolute_delta = std::abs(delta_cost);
                                    best_cost = new_cost;
                                    best_i = i;
                                    best_j = j;
                                    best_from = from;
                                    best_insert_pos = insert_pos;
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
    }

    if (best_i != -1 && best_j != -1)
    {
        current_routes[best_i] = best_routes[0];
        current_routes[best_j] = best_routes[1];

        graph.trucksvector[best_i].route = best_routes[0];
        graph.trucksvector[best_j].route = best_routes[1];

        graph.trucksvector[best_i].current_time = best_time_for_route1;
        graph.trucksvector[best_j].current_time = best_time_for_route2;

        current_cost = best_cost;

        used_ops[0] = 2;
        used_ops[1] = std::fabs(best_delta);

        add_to_Tabu(used_ops);
    }
}

void Sa::route_splitting_and_merging(double &current_cost, int &current_used_trucks, std::vector<std::vector<int>> &current_routes, double *used_ops, Graph &graph)
{
    double best_absolute_delta = 0;
    double best_cost = current_cost;
    std::vector<std::vector<int>> best_routes = current_routes;
    int best_split_idx = -1, best_merge_idx = -1, best_split_point = -1, best_insert_pos = -1;
    double best_first_route_time = -1, best_merged_route_time = -1, best_delta = -1;

    for (int idx = 0; idx < current_routes.size(); ++idx)
    {
        if (current_routes[idx].size() <= 2)
            continue;

        int split_point = rand() % (current_routes[idx].size() - 1) + 1;

        std::vector<int> first_part(current_routes[idx].begin(), current_routes[idx].begin() + split_point);
        std::vector<int> second_part(current_routes[idx].begin() + split_point, current_routes[idx].end());

        double first_route_time = is_single_route_possible(first_part, graph);
        double second_route_time = is_single_route_possible(second_part, graph);

        if (first_route_time == -1 || second_route_time == -1)
            continue;

        for (int merge_idx = 0; merge_idx < current_routes.size(); ++merge_idx)
        {
            if (merge_idx == idx)
                continue;

            for (int insert_pos = 0; insert_pos <= current_routes[merge_idx].size(); ++insert_pos)
            {
                std::vector<int> test_merged_route = current_routes[merge_idx];
                test_merged_route.insert(test_merged_route.begin() + insert_pos, second_part.begin(), second_part.end());

                double merged_route_time = is_single_route_possible(test_merged_route, graph);

                if (merged_route_time != -1)
                {
                    used_ops[0] = 5;

                    if (current_routes[merge_idx].empty() || current_routes[idx].empty())
                        continue;

                    double old_merge_time = graph.trucksvector[merge_idx].current_time;
                    double old_idx_time = graph.trucksvector[idx].current_time;
                    double new_cost = current_cost - old_merge_time - old_idx_time - graph.distances[0][current_routes[merge_idx].back()] - graph.distances[0][current_routes[idx].back()] + first_route_time + merged_route_time + graph.distances[0][first_part.back()] + graph.distances[0][test_merged_route.back()];

                    double delta_cost = new_cost - current_cost;
                    used_ops[1] = delta_cost;

                    if (not_in_Tabu(used_ops))
                    {
                        if (delta_cost < 0 || accept_worse_solution(delta_cost, temperature))
                        {
                            if (std::abs(delta_cost) > std::fabs(best_delta))
                            {
                                best_cost = new_cost;
                                best_routes = current_routes;
                                best_routes[idx] = first_part;
                                best_routes[merge_idx] = test_merged_route;

                                best_split_idx = idx;
                                best_merge_idx = merge_idx;
                                best_split_point = split_point;
                                best_insert_pos = insert_pos;
                                best_delta = delta_cost;

                                best_first_route_time = first_route_time;
                                best_merged_route_time = merged_route_time;
                            }
                        }
                    }
                }
            }
        }
    }

    if (best_split_idx != -1 && best_merge_idx != -1)
    {
        current_cost = best_cost;
        current_routes = best_routes;

        graph.trucksvector[best_split_idx].route = best_routes[best_split_idx];
        graph.trucksvector[best_merge_idx].route = best_routes[best_merge_idx];

        graph.trucksvector[best_split_idx].current_time = best_first_route_time;
        graph.trucksvector[best_merge_idx].current_time = best_merged_route_time;

        used_ops[0] = 5;
        used_ops[1] = std::fabs(best_delta);

        add_to_Tabu(used_ops);
    }
}

void Sa::reduce_truck_count(double &current_cost, int &current_used_trucks, std::vector<std::vector<int>> &current_routes, double *used_ops, Graph &graph)
{
    double best_cost = current_cost;
    bool truck_reduction_achieved = false;

    // sort po skomplikowaniu tras
    std::vector<std::pair<int, double>> route_utilization;
    for (int i = 0; i < current_routes.size(); ++i)
    {
        if (!current_routes[i].empty())
        {
            double complexity_score = graph.trucksvector[i].current_time * (1.0 + current_routes[i].size() * 0.1);
            route_utilization.push_back({i, complexity_score});
        }
    }
    std::sort(route_utilization.begin(), route_utilization.end(),
              [](const auto &a, const auto &b)
              { return a.second < b.second; });

    for (int merge_attempts = 0; merge_attempts < current_used_trucks; ++merge_attempts)
    {
        bool reduction_in_this_attempt = false;

        // od poczatku posortowania sa latwiejsze trasy do polaczenai
        for (int source_route_idx = 0; source_route_idx < route_utilization.size(); source_route_idx++)
        {
            if (current_routes[source_route_idx].empty())
                continue;

            int best_target_route = -1;
            int best_insert_pos = -1;
            double best_merge_time = -1;
            double best_delta_cost = std::numeric_limits<double>::max();

            for (int target_route = 0; target_route < current_routes.size(); ++target_route)
            {
                if (target_route == source_route_idx || current_routes[target_route].empty())
                    continue;

                std::vector<int> insertion_points;
                for (int i = 0; i < static_cast<int>(current_routes[target_route].size()); i++)
                {
                    insertion_points.push_back(i);
                }
                for (int insert_pos : insertion_points)
                {
                    std::vector<int> test_source_route = current_routes[source_route_idx];
                    std::vector<int> test_target_route = current_routes[target_route];

                    // wkladanie sciezki do innej w 3 miejscach
                    test_target_route.insert(test_target_route.begin() + insert_pos,
                                             test_source_route.begin(),
                                             test_source_route.end());

                    double target_route_time = is_single_route_possible(test_target_route, graph);

                    if (target_route_time != -1)
                    {
                        double old_source_time = graph.trucksvector[source_route_idx].current_time;
                        double old_target_time = graph.trucksvector[target_route].current_time;

                        double new_cost = current_cost - old_source_time - old_target_time - graph.distances[0][current_routes[source_route_idx].back()] - graph.distances[0][current_routes[target_route].back()] + target_route_time + graph.distances[0][test_target_route.back()];
                        double delta_cost = new_cost - current_cost;

                        double time_balance_penalty = std::abs(old_source_time - old_target_time);
                        double complexity_penalty = test_source_route.size() * 0.5;
                        double merge_score = delta_cost + time_balance_penalty + complexity_penalty;
                        used_ops[0] = 6;
                        used_ops[1] = delta_cost;

                        if (merge_score < best_delta_cost)
                        {
                            best_target_route = target_route;
                            best_insert_pos = insert_pos;
                            best_merge_time = target_route_time;
                            best_delta_cost = merge_score;
                        }
                    }
                }
            }

            if (best_target_route != -1)
            {
                current_routes[best_target_route].insert(
                    current_routes[best_target_route].begin() + best_insert_pos,
                    current_routes[source_route_idx].begin(),
                    current_routes[source_route_idx].end());
                current_routes[source_route_idx].clear();

                graph.trucksvector[best_target_route].route = current_routes[best_target_route];
                graph.trucksvector[best_target_route].current_time = best_merge_time;
                graph.trucksvector[source_route_idx].route.clear();
                graph.trucksvector[source_route_idx].current_time = 0;

                current_cost = current_cost + best_delta_cost;
                current_used_trucks--;
                reduction_in_this_attempt = true;
                truck_reduction_achieved = true;

                route_utilization.erase(std::remove_if(route_utilization.begin(), route_utilization.end(), [source_route_idx](const auto &pair)
                                                       { return pair.first == source_route_idx; }),
                                        route_utilization.end());

                break;
            }
        }

        if (!reduction_in_this_attempt)
            break;
    }

    if (truck_reduction_achieved)
    {
        used_ops[0] = 6;
        used_ops[1] = current_cost - best_cost;
        add_to_Tabu(used_ops);
    }
}
