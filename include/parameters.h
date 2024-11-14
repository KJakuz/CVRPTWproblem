#ifndef PARAMETERS_H
#define PARAMETERS_H

class Parameters {
    public:
        int distance_cost_param = 60;
        int window_time_param = 1;
        int waiting_time_param = 40;
        int demand_param = 40;
        int RCLpercent = 1;
        int time_limit_in_seconds = 10;

        Parameters() {}

        void setAlfa(int value);
        void setBeta(int value);
        void setGamma(int value);
        void setRCLpercent(int value);
        void setTimeLimit(int value);

};


#endif