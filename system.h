#ifndef SYSTEM_H
#define SYSTEM_H

#include <vector>
#include "Eigen/Dense"
#include "help.h"
#include <iostream>

using namespace std;
using namespace Eigen;

class SP_State
{
public:
    SP_State(double _spEnergy) : spEnergy(_spEnergy){}
    double spEnergy; //single particle energy
};

class System
{
public:
    virtual ~System(){}
    System(int,int);
    int A; //number of particles, have to be even for pairing
    int numberSP; //number of single particle orbitals

    //single particle orbitals
    virtual void generateSP_States(int) = 0;
    vector<SP_State*> SP_States;

    //configurations
    virtual void generateConfigurations() = 0;
    virtual void printConfigurations() = 0;
    vector<int> configurations;

    //matrix elements
    virtual double V1B(int,int) = 0;
    virtual double V2B(int,int,int,int) = 0;

    //diagonalization
    MatrixXd H;
    void generateH();
    double getH(int,int);
    void diagonalization();
    void diag_printCoefficients();
    double diag_E_GS;
    double diag_deltaE;
    VectorXd diag_coefficients;

    //MBPT
    void MBPT();
    double MBPT_getCoefficient(int,int,int,int);
    void MBPT_printCoefficients();
    void MBPT_calculateDeltaE();
    vector<double> MBPT_coefficients;
    double MBPT_E_GS;
    double MBPT_deltaE;

    //CCD
    void CCD_generateMatrices();
    void CCD_calculateTau();
    MatrixXd CCD_V_ph, CCD_V_pp, CCD_V_hh, CCD_e_ph, CCD_Tau, CCD_t_m;
    double CCD_t(int,int,int,int);
    double CCD_E_GS;
    double CCD_deltaE;

    //GF
    void GF_generateMatrices();
    void GF_diag();
    double GF_E_GS;
    double GF_deltaE;
    MatrixXd GF_Sigma_static;
    MatrixXd GF_Mat;
    MatrixXd GF_Mstatic;
    MatrixXd GF_Mr;
    MatrixXd GF_Mq;
    MatrixXd GF_Er;
    MatrixXd GF_Eq;
    MatrixXd GF_Z;
    MatrixXd GF_W;
    VectorXd GF_e;
};

#endif
