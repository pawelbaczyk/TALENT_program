#ifndef PAIRING_H
#define PAIRING_H

#include "system.h"
#include <bitset>

class Pairing : public System
{
public:
    Pairing(int,int,double);

    //parameters
    double g;

    //single particle states
    void generateSP_States(int);

    //configurations
    void generateConfigurations();
    void printConfigurations();

    //matrix elements
    double V1B(int,int);
    double V2B(int,int,int,int);
};

class SP_Pairing : public SP_State
{
public:
    SP_Pairing(int _p, bool _spin): SP_State(_p-1), p(_p), spin(_spin){}
    int p; //from 1
    bool spin; //0 -- spin up, 1 -- spin down
};

#endif
