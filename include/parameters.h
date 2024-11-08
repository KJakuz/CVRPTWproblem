#ifndef PARAMETERS_H
#define PARAMETERS_H

class Parameters {
    public:
        int distance_cost_param = 35;
        int window_time_param = 1;
        int waiting_time_param = 35;
        int RCLpercent = 2;
        int time_limit_in_seconds = 300;

        Parameters() {}

        void setAlfa(int value);
        void setBeta(int value);
        void setGamma(int value);
        void setRCLpercent(int value);
        void setTimeLimit(int value);

};


#endif