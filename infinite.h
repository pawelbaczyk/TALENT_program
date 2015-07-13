#ifndef INFINITE_H
#define INFINITE_H

#include "system.h"
#include "help.h"

const double hbar = 197;//TODO
const double m = 938;

class Infinite : public System
{
public:
    Infinite(int, int);

    //parameters
    //

    //single particle states
    void generateSP_States(int);

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
    SP_Infinite(double _kx, double _ky, double _kz, bool _spin, bool _isospin)
        : SP_State(hbar*hbar/2.0/m*(_kx*_kx + _ky*_ky + _kz*_kz)),
          kx(_kx), ky(_ky), kz(_kz), spin(_spin), isospin(_isospin) {}
    double kx;
    double ky;
    double kz;
    bool spin; //0 -- spin up, 1 -- spin down
    bool isospin; //0 -- neutron, 1 -- proton
};

#endif
