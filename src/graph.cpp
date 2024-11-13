#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <random>
#include "graph.h"

class Truck;
Parameters defaultParameters;

// tworzymy wektor obiektow ciezarowek
void Graph::init_trucks(Truck &truckinfo)
{
    for (int i = 0; i < truckinfo.trucks_number + 1; i++)
    {
        trucksvector.push_back(Truck(truckinfo.trucks_number, i, truckinfo.capacity, truckinfo.capacity, 0, 0));
    }
}

// liczenie macierzy odleglosci od węzłów
void Graph::measure_distances()
{
    // TODO: Moze lepiej by bylo zastosowac std::array?
    distances.resize(number_nodes, std::vector<double>(number_nodes, 0.0f)); // alokacja pamiec, wielkosc macierzy number_nodes i zapisanej zerami
    for (int i = 0; i < number_nodes; i++)
    {
        for (int j = 0; j < number_nodes; j++)
        {
            // wzor na odleglosc na plaszczyznie euklidesowej
            distances[i][j] = sqrt(pow(Nodes[j].xcord - Nodes[i].xcord, 2) + pow(Nodes[j].ycord - Nodes[i].ycord, 2));
        }
    }
}

// pomocnicze
void Graph::show_distances_matrix()
{
    std::cout << "macierz odleglosci: \n   ";
    for (int k = 0; k < number_nodes; k++)
    {
        std::cout << k << "  ";
    }
    std::cout << std::endl;
    for (int i = 0; i < number_nodes; i++)
    {
        std::cout << i << ": ";
        for (int j = 0; j < number_nodes; j++)
        {
            std::cout << distances[i][j] << " ";
            // std::cout<<"distances["<<i<<"]"<<"["<<j<<"]: "<<distances[i][j]<<"\n";
        }
        std::cout << "\n";
    }
}

// pomocnicze
void Graph::show_nodes_values()
{
    for (int i = 0; i < number_nodes; i++)
    {
        Node &node = Nodes[i];
        std::cout << "ID: " << node.id << ", X: " << node.xcord
                  << ", Y: " << node.ycord
                  << ", Zapotrzebowanie: " << node.demand
                  << ", Czas otwarcia: " << node.readytime
                  << ", Czas zamknięcia: " << node.duedate
                  << ", Czas obsługi: " << node.servicetime << std::endl;
    }
}
// pomocnicze
void Graph::show_one_node_values(Node &node)
{
    std::cout << "ID: " << node.id << ", X: " << node.xcord
              << ", Y: " << node.ycord
              << ", Zapotrzebowanie: " << node.demand
              << ", Czas otwarcia: " << node.readytime
              << ", Czas zamknięcia: " << node.duedate
              << ", Czas obsługi: " << node.servicetime << std::endl;
};

bool Graph::all_visited()
{
    if (unvisited.size() == 0)
    {
        return true;
    }
    return false;
}

// wektor nieodwiedzonych wezlow
void Graph::makeunvisitedvector()
{
    for (int i = 1; i < number_nodes; i++)
    {
        if (Nodes[i].check_if_done == false)
        {
            unvisited.push_back(Nodes[i]);
        }
    }
}

// do sortowania
bool compareCandidates(std::pair<Node, double> &a, std::pair<Node, double> &b)
{
    return a.second < b.second; // rosnaco wedlug kosztow
}

// glowny algorytm
int Graph::GRASP()
{
    std::ofstream outputFile("wyniki.txt");
    if (!outputFile.is_open())
    {
        std::cerr << "Nie można otworzyć pliku" << std::endl;
        return -1;
    }

    std::vector<std::vector<int>> solution;
    int idx_of_truck = 0;
    makeunvisitedvector();

    while (!all_visited())
    {

        if (idx_of_truck >= trucksvector.size())
        {
            trucksvector.push_back(Truck(0, idx_of_truck, trucksvector[0].capacity, trucksvector[0].capacity, 0, 0));
        }

        std::vector<std::pair<Node, double>> candidates; // lista kandydatow w formie <<wezel>,<koszt>>

        // sprawdzanie mozliwych kandydatow
        for (int i = 0; i < unvisited.size(); i++)
        {
            if ((trucksvector[idx_of_truck].capacity >= unvisited[i].demand) && (distances[0][unvisited[i].id] < unvisited[i].duedate))
            {
                if (trucksvector[idx_of_truck].which_node != unvisited[i].id) // TODO: czy to jest potrzebne?
                {
                    if (trucksvector[idx_of_truck].cargo >= unvisited[i].demand)
                    {
                        if (trucksvector[idx_of_truck].current_time + distances[trucksvector[idx_of_truck].which_node][unvisited[i].id] < unvisited[i].duedate)
                        {
                            // obliczanie kosztow kandydatow
                            double waiting_time_costs = std::max(0.0, unvisited[i].readytime - (trucksvector[idx_of_truck].current_time + distances[trucksvector[idx_of_truck].which_node][unvisited[i].id]));
                            double window_time_costs = std::max(0.0, unvisited[i].duedate - (trucksvector[idx_of_truck].current_time + distances[trucksvector[idx_of_truck].which_node][unvisited[i].id]));
                            double demand_cost = (trucksvector[idx_of_truck].capacity - trucksvector[idx_of_truck].cargo) - unvisited[i].demand;
                            double cost = defaultParameters.distance_cost_param * distances[trucksvector[idx_of_truck].which_node][unvisited[i].id] +
                                          +defaultParameters.demand_param * demand_cost + defaultParameters.window_time_param * window_time_costs + defaultParameters.waiting_time_param * waiting_time_costs;
                            candidates.push_back(std::make_pair(unvisited[i], cost));
                        }
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
        std::sort(candidates.begin(), candidates.end(), compareCandidates);

        // wybieramy losowo wezel z pierwszych iles procent kandydatow
        int limit = std::max(static_cast<int>(candidates.size() * (defaultParameters.RCLpercent / 100.0)), 1);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, limit - 1);
        int next_node_index = dis(gen);

        // aktualizujemy atrybuty ciezarowki
        trucksvector[idx_of_truck].route.push_back(candidates[next_node_index].first.id);
        trucksvector[idx_of_truck].cargo -= candidates[next_node_index].first.demand;
        double waiting_time = std::max(0.0, candidates[next_node_index].first.readytime - (trucksvector[idx_of_truck].current_time + distances[trucksvector[idx_of_truck].which_node][candidates[next_node_index].first.id]));
        trucksvector[idx_of_truck].current_time += distances[trucksvector[idx_of_truck].which_node][candidates[next_node_index].first.id] + candidates[next_node_index].first.servicetime + waiting_time;
        trucksvector[idx_of_truck].which_node = candidates[next_node_index].first.id;
        candidates[next_node_index].first.check_if_done = true;

        // usuniecie z unvisited wezla przez ktory przejezdzamy (znajdujemy id do usuniecia w wektorze)
        unvisited.erase(std::remove_if(unvisited.begin(), unvisited.end(), [&](Node &node)
                                       { return node.id == candidates[next_node_index].first.id; }),
                        unvisited.end());
    }
    // zwracamy liczbe tirow
    return idx_of_truck;
}

void Graph::reset_trucks()
{
    for (int i = 0; i < number_nodes; i++)
    {
        Nodes[i].check_if_done = false;
    }

    unvisited.clear();
    makeunvisitedvector();

    Truck truck_info = trucksvector[0];
    trucksvector.clear();
    init_trucks(truck_info);
}

// wielokrotne uruchamianie graspa
void Graph::rungrasp()
{
    double best_cost = 99999999;
    int best_trucks_number = 99999999;
    std::vector<std::vector<int>> best_routes;

    auto start_time = std::chrono::steady_clock::now();
    auto end_time = start_time + std::chrono::seconds(defaultParameters.time_limit_in_seconds);
    // procent do wyboru kandydata z rcl jest adaptacyjny, na poczatku programu 1 -> 6 -> 11 -> 16 aby rozwijac mozliwe odpowiedzi
    auto quarter_duration = (end_time - start_time) / 4;
    auto change_rcl_percent = start_time + quarter_duration;

    while (std::chrono::steady_clock::now() < end_time)
    {
        reset_trucks();
        int number_of_trucks = GRASP();

        if (std::chrono::steady_clock::now() >= change_rcl_percent)
        {
            defaultParameters.RCLpercent += 5;
            change_rcl_percent += quarter_duration;
        }
        double cost = 0;
        for (int i = 0; i < trucksvector.size(); i++)
        {
            cost += trucksvector[i].current_time + distances[trucksvector[i].which_node][0];
        }

        if (number_of_trucks <= best_trucks_number)
        {
            if (cost < best_cost)
            {
                best_cost = cost;
                best_trucks_number = number_of_trucks;

                best_routes.clear();
                for (int l = 0; l < trucksvector.size(); l++)
                {
                    best_routes.push_back(trucksvector[l].route);
                }
            }
        }
    }

    std::ofstream outputFile("wyniki.txt");
    if (!outputFile.is_open())
    {
        std::cerr << "Nie można otworzyć pliku" << std::endl;
        return;
    }

    outputFile << std::fixed << std::setprecision(5) << best_trucks_number + 1 << " " << best_cost << std::endl;
    for (int l = 0; l < best_routes.size(); l++)
    {
        if (!best_routes[l].empty())
        {
            for (int k = 0; k < best_routes[l].size(); k++)
            {
                outputFile << best_routes[l][k] << " ";
            }
            outputFile << std::endl;
        }
    }

    outputFile.close();
}