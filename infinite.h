#ifndef INFINITE_H
#define INFINITE_H

#include "system.h"

class Infinite : public System
{
public:
    Infinite();

    //parameters
    //

    //single particle states
    void generateSP_States();

    //configurations
    void generateConfigurations();
    void printConfigurations();

    //matrix elements
    double V1B(int,int);
    double V2B(int,int,int,int);
};

class SP_Infinite : public SP_State
{
public:
    SP_Infinite(int _p, bool _spin): SP_State(_p-1), p(_p), spin(_spin){}
    int p; //from 1
    bool spin; //0 -- spin up, 1 -- spin down
};

#endif
