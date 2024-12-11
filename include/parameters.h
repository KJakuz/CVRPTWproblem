#ifndef PARAMETERS_H
#define PARAMETERS_H

class Parameters {
    public:
        //grasp params
        int distance_cost_param = 60;
        int window_time_param = 1;
        int waiting_time_param = 20;
        int demand_param = 20;
        double RCLpercent = 0.5;

        //time param
        int time_limit_in_seconds = 5;

        //tabu search params
        static const int Tabu_list_size = 20;
        int no_improvement_limit = Tabu_list_size;
        //simulated annealing params
        double temperature = 100;
        double cooling_factor = 0.945;
        double min_temperature = 0.1;

        //operators param
        double op1=1;
        double op2=1;
        double op3=1;
        double op4=1;
        double op5=3;




        Parameters() {}

};


#endif