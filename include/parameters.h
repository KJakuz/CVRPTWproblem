#ifndef PARAMETERS_H
#define PARAMETERS_H

class Parameters {
    public:
        int distance_cost_param = 25;
        int window_time_param = 1;
        int waiting_time_param = 25;
        int RCLpercent = 1;
        int time_limit_in_seconds = 180;

        Parameters() {}

        void setAlfa(int value);
        void setBeta(int value);
        void setGamma(int value);
        void setRCLpercent(int value);
        void setTimeLimit(int value);

};


#endif