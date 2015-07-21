#ifndef SYSTEM_H
#define SYSTEM_H

#include <vector>
#include "Eigen/Dense"
#include "help.h"
#include <iostream>
#include "Eigen/SparseCore"

#include <iomanip>

using namespace std;
using namespace Eigen;

class SP_State
{
public:
    SP_State(double _spEnergy) : spEnergy(_spEnergy){}
    double spEnergy; //single particle energy
};

class TwoBody_State
{
 public:
 TwoBody_State(int _p,int _q):p(_p),q(_q){}
  int p,q;
};

class System
{
 public:
  virtual ~System();
  System(int);
  int A; //number of particles, have to be even for pairing
  int numberSP; //number of single particle orbitals

  //single particle orbitals
  virtual void generateSP_States(int) = 0;
  vector<SP_State*> SP_States;


  //two_body states
  virtual void generateTwoBody_States()=0;
  vector<TwoBody_State*> TwoBody_States_hh;
  vector<TwoBody_State*> TwoBody_States_hp;
  vector<TwoBody_State*> TwoBody_States_pp;

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
  void CCD_SparseMatrices();
  void CCD_OnFlight();
  double CCD_OnFlight_t(int,int,int,int);
  MatrixXd CCD_OnFlight_t_m;
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
