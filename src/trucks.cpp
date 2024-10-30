#include "trucks.h"
#include <iostream>

    //mozna rozladowac w nastepnym węźle?
bool Truck::check_time(const Node& node,const std::vector<std::vector<float>>& distances){ //const Node& node tworzy stala referencje tj. nie mozna modyfikowac
    float distance_to_next = distances[which_node][node.id]; //moze -1
    if((current_time + distance_to_next <= node.duedate) && (node.readytime >= current_time + distance_to_next)){
        return true;
    }
    else{
        return false;
    }
}

    //starczy ładunku w następnym węźle?
bool Truck::check_cargo(const Node& node){
    if(cargo >= node.demand){
        return true;
    }
    else{
        return false;
    }
}

void Truck::show_trucks_info(){
    std::cout << "Liczba pojazdów: " << trucks_number << std::endl;
    std::cout << "Pojemność pojazdu: " << capacity << std::endl;
}