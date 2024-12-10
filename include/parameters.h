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
        int time_limit_in_seconds = 300;

        //tabu search params
        static const int Tabu_list_size = 20;
        int no_improvement_limit = Tabu_list_size;
        //simulated annealing params
        double temperature = 500;
        double cooling_factor = 0.975;
        double min_temperature = 0.1;


        Parameters() {}

};


#endif