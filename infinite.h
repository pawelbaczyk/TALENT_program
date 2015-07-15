#ifndef INFINITE_H
#define INFINITE_H

#include "system.h"
#include "help.h"
#include <iomanip>
#include <map>
#include <cstdio>
#include <string>

const double hbar = 197;//TODO
const double m = 938;

const double V0R = 200;
const double KR = 1.487;
const double V0T = -178;
const double KT = 0.639;
const double V0S = -91.85;
const double KS = 0.465;

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

class Infinite : public System
{
public:
    Infinite(int,int,double,int);

    //parameters
    int g_s;
    double k_F;
    double rho;
    int nMax;
    double L;

    //gauss-legendre points and weights in range [-1,1]
    vector<double> gauss_x;
    vector<double> gauss_w;

    void setRho(double);

    //single particle states
    void generateSP_States(int);
    void printSP_States();

    //configurations
    void generateConfigurations();
    void printConfigurations();

    //matrix elements
    double V1B(int,int);
    double V2B(int,int,int,int);
    double V2B_sym(int,int,int,int);
    bool deltaSpin(SP_Infinite*, SP_Infinite*, SP_Infinite*, SP_Infinite*);
    bool deltaIsospin(SP_Infinite*, SP_Infinite*, SP_Infinite*, SP_Infinite*);
    bool deltaSpinIsospin(SP_Infinite*, SP_Infinite*, SP_Infinite*, SP_Infinite*);

    //HF
    void HF_calculateE0();
    double HF_E0;

    //exact HF
    double HF_exact_f(double r);
    void HF_cal_exact_E0();
    double HF_exact_E0;

    //maps
    void map_generateV2B();
    map<string,double> map_V2B;
};

#endif
