#include "parameters.h"
#include <iostream>

void Parameters::setAlfa(int value){
    distance_cost_param = value;
}

void Parameters::setBeta(int value){
    window_time_param = value;
}

void Parameters::setGamma(int value){
    waiting_time_param = value;
}

void Parameters::setRCLpercent(int value){
    RCLpercent = value;
}

void Parameters::setTimeLimit(int value){
    time_limit_in_seconds = value;
}