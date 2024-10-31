#include "loadfile.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

bool LoadFile::load(const std::string& filename, Graph& graph, Truck& truck) {
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "nie można otworzyć pliku " << filename << std::endl;
        return false;
    }

    std::string line;
    int vehicle_number = 0;
    int vehicle_capacity = 0;

    // sekcja VEHICLE
    while (std::getline(file, line)) {
        if (line.find("VEHICLE") != std::string::npos) { //npos to stala no position, czyli tutaj brak w linii
            std::getline(file, line); // Pomijamy linię "NUMBER CAPACITY"
            file >> vehicle_number >> vehicle_capacity;
            truck.trucks_number = vehicle_number;
            truck.capacity = vehicle_capacity;
            break;
        }
    }

    // sekcja CUSTOMER
    std::vector<Node> nodes;
    while (std::getline(file, line)) {
        if (line.find("CUSTOMER") != std::string::npos) {
            std::getline(file, line); 
            int id, x, y, demand, ready_time, due_date, service_time;
            
            while (file >> id >> x >> y >> demand >> ready_time >> due_date >> service_time) {
                Node node;
                node.id = id;
                node.xcord = x;
                node.ycord = y;
                node.demand = demand;
                node.readytime = ready_time;
                node.duedate = due_date;
                node.servicetime = service_time;
                node.check_if_done = false; 
                
                nodes.push_back(node);
            }
            break;
        }
    }

    graph.Nodes = nodes;
    graph.number_nodes = nodes.size();

    file.close();
    return true;
}
