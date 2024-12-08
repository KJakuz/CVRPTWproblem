#ifndef PARAMETERS_H
#define PARAMETERS_H

class Parameters {
    public:
        //grasp params
        int distance_cost_param = 60;
        int window_time_param = 1;
        int waiting_time_param = 20;
        int demand_param = 20;
        int RCLpercent = 1;

        //time param
        int time_limit_in_seconds = 300;

        //tabu search param
        int no_improvement_limit = 4000;
        static const int Tabu_list_size = 500;

        //simulated annealing
        double temperature = 1000;
        double cooling_factor = 0.98;
        double min_temperature = 0.1;


        Parameters() {}

};


#endif