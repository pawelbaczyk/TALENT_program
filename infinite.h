#ifndef INFINITE_H
#define INFINITE_H

#include "system.h"
#include "help.h"
#include <iomanip>
#include <map>
#include <cstdio>
#include <string>


const double hbar = 197.326968;//TODO
const double mn = 939.565;//n
const double mp = 938.27200;
const double m = mn;

const double V0R = 200;
const double KR = 1.487;
const double V0T = -178;
const double KT = 0.639;
const double V0S = -91.85;
const double KS = 0.465;

class SP_Infinite : public SP_State
{
public:
    SP_Infinite(int _nx,int _ny, int _nz, double _kx, double _ky, double _kz, int _spin, int _isospin)
        : nx(_nx),ny(_ny),nz(_nz),SP_State(hbar*hbar/2.0/m*(_kx*_kx + _ky*_ky + _kz*_kz)),
          kx(_kx), ky(_ky), kz(_kz), spin(_spin), isospin(_isospin),
          HF_spEnergy(hbar*hbar/2.0/m*(_kx*_kx + _ky*_ky + _kz*_kz)) {}
    double kx;
    double ky;
    double kz;
    int nx,ny,nz;
    int spin; //1 -- spin up, -1 -- spin down
    int isospin; //1 -- neutron, -1 -- proton
    double HF_spEnergy;//init to spEnergy
};

class TwoBody_Infinite: public TwoBody_State
{
public:
    TwoBody_Infinite(int _p,int _q,int _Nx,int _Ny,int _Nz,int _Sz,int _Tz):TwoBody_State(_p,_q),Nx(_Nx),Ny(_Ny),Nz(_Nz),Sz(_Sz),Tz(_Tz){}
    int Nx;
    int Ny;
    int Nz;
    int Sz;
    int Tz;
};

class TwoBody_compare
{
public:
    bool operator()(TwoBody_State*,TwoBody_State*);
};


class Channel
{
public:
    Channel(int _Nx,int _Ny,int _Nz,int _Sz,int _Tz):Nx(_Nx),Ny(_Ny),Nz(_Nz),Sz(_Sz),Tz(_Tz){}
    int Nx,Ny,Nz,Sz,Tz;
    vector<int> bra,ket;
    MatrixXd mat;
    bool Include(TwoBody_Infinite*);
};


struct Position
{
public:
    Position() : channel(-1), index(0){}
    Position(int _channel, int _index) : channel(_channel), index(_index){}
    int channel;
    int index;
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


    //twobody states
    void generateTwoBody_States();
    void printTwoBody_States();

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

    //CCD
    void CCD_generateBlockMatrices();
    void CCD_BlockMatricesLadders();
    void CCD_BlockMatricesIntermediates();
    vector<Channel> CCD_V_hhhh, CCD_V_hhpp,CCD_V_pppp,CCD_T_hhpp,CCD_e_hhpp,CCD_V_hphp;
    Matrix<Position,Dynamic,Dynamic> CCD_position;
};

#endif
