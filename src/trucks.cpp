#include "trucks.h"
#include <iostream>
#include <fstream>

void Truck::show_trucks_info()
{
    std::cout << "Liczba pojazdów: " << trucks_number << std::endl;
    std::cout << "Pojemność pojazdu: " << capacity << std::endl;
    std::cout << "aktualnie towaru: " << cargo << std::endl;
    std::cout << "który węzeł: " << which_node << std::endl;
    std::cout << "jaki czas: " << current_time << std::endl;
}