#include <iostream>
#include <vector>
#include <cmath>
#include <bitset>
#include "help.h"
#include "Eigen/Dense"
#include <iomanip>

using namespace std;
using namespace Eigen;

//////////////////////////////

struct PairingSP
{
public:
    PairingSP(int, bool);
    int p; //from 1
    bool spin; //0 -- spin up, 1 -- spin down
    double spEnergy; //single particle energy
};

PairingSP::PairingSP(int _p, bool _spin) : p(_p), spin(_spin)
{
    spEnergy = p-1;
}

//////////////////////////////

class Pairing
{
public:
    vector<PairingSP> SP_States;
    void generateSP(int);
    int PMax;
    double g;
};

void Pairing::generateSP(int _PMax)//generate s.p. states (p orbits, spin up/down)
{
    PMax = _PMax;
    for (int p = 1; p <= PMax; p++)
    {
        SP_States.push_back(PairingSP(p,0));//spin up
        SP_States.push_back(PairingSP(p,1));//spin down
    }
}

//////////////////////////////

class System
{
public:
    System(int,int,double);
    int A; //number of particles, have to be even for pairing

    //pairing
    Pairing spOrbitals;
    void generateSPOrbitals(int);
    void set_g(double);

    //configurations
    void generateConfigurations();
    void printConfigurations();
    vector<int> configurations;

    //system
    MatrixXd H;
    double V2B(int,int,int,int);
    double V1B(int,int);
    double getH(int,int);
    void generate_H();

    //diagonalization
    void diagonalization();
    void diag_printCoefficients();
    double diag_E_GS;
    VectorXd diag_coefficients;

    //MBPT
    void MBPT();
    void MBPT_calculateCoefficients();
    double MBPT_getCoefficient(int,int,int,int);
    void MBPT_printCoefficients();
    void MBPT_calculateDeltaE();
    vector<double> MBPT_coefficients;
    double MBPT_deltaE;
};

void System::generateSPOrbitals(int PMax)//generate s.p. orbitals
{
    spOrbitals.generateSP(PMax);
}

void System::set_g(double g)//set pairing g value
{
    spOrbitals.g = g;
}

System::System(int _A, int PMax, double g) : A(_A)//generate system with A particles
{
    generateSPOrbitals(PMax);
    set_g(g);
}

void System::diagonalization()//perform diagonalization
{
  generate_H();
  SelfAdjointEigenSolver<MatrixXd> solver(H);
  if(solver.info()!=Success) abort();
  diag_E_GS=solver.eigenvalues()[0];
  diag_coefficients=solver.eigenvectors().col(0);
}

void System::diag_printCoefficients()//print ground state eigenvector with C_0 = 1
{
    for (int i = 0; i < diag_coefficients.size(); i++)
        cout << setw(8) << diag_coefficients[i]/diag_coefficients[0] << "\t";
    cout << endl;
}

double System::V2B(int i,int j,int a,int b)//calculate 2body matrix element <ij|V_pairing|ab>
{
    if ((spOrbitals.SP_States[a].p == spOrbitals.SP_States[b].p) && //same 'p'
        (spOrbitals.SP_States[a].spin != spOrbitals.SP_States[b].spin) &&//different spin
        (spOrbitals.SP_States[i].p == spOrbitals.SP_States[j].p) &&
        (spOrbitals.SP_States[i].spin != spOrbitals.SP_States[j].spin))
        return -0.5*spOrbitals.g;
    else
        return 0;
}
double System::V1B(int a,int b)
{
  if((spOrbitals.SP_States[a].p == spOrbitals.SP_States[b].p) && 
     (spOrbitals.SP_States[a].spin == spOrbitals.SP_States[b].spin))
    return spOrbitals.SP_States[a].p-1;
}
void System::generate_H()
{
  double g=spOrbitals.g;
  H.resize(6,6);
  for(int i=0;i<6;i++)
    for(int j=0;j<6;j++)
      {
	if(i!=j)
	  H(i,j)=-0.5*g;
	if(i+j==6)
	  H(i,j)=0;
	if(i==j)
	  H(i,j)=2*(i+1)-g;
      }
}
double System::getH(int a,int b)
{
  int bra=configurations[a];
  int ket=configurations[b];
  int diffCount=0;
  vector<int> OccupiedPos;
  vector<int> UnoccupiedPos;
  vector<int> NumOrbitalsKet;
  vector<int> NumOrbitalsBra;
  int numKet(0),numBra(0);
  for(int i=0;i<spOrbitals.SP_States.size();i++)
    {
      if(getBit(bra,i)!=getBit(ket,i))
	{
	  ++diffCount;
	  if(diffCount>4) break;
	  if(getBit(ket,i))
	    {
	      OccupiedPos.push_back(i);
	      NumOrbitalsKet.push_back(numKet);
	    }
	  else
	    {
	      UnoccupiedPos.push_back(i);
	      NumOrbitalsBra.push_back(numBra);
	    }
	}
      if(getBit(ket,i)) ++numKet;
      if(getBit(bra,i)) ++numBra;
    }
  if(diffCount>4) return 0;
  if(diffCount==4)
    {
      int phase=(NumOrbitalsBra[0]+NumOrbitalsBra[1]+NumOrbitalsKet[0]+NumOrbitalsKet[1])%2?-1:1;
      return phase*V2B(UnoccupiedPos[0],UnoccupiedPos[1],OccupiedPos[0],OccupiedPos[1]);
    }
  else if(diffCount==2)
    {
      int phase=(NumOrbitalsBra[0]+NumOrbitalsKet[0])%2?-1:1;
      return phase*V1B(UnoccupiedPos[0],OccupiedPos[0]);
    }
  else if(diffCount==0)
    {
      double temp;
      for(int i=0;i<spOrbitals.SP_States.size();i++)
	{
	  if(getBit(ket,i))
	    {
	      //	      temp+=V1B(i,i)+0.5*
	    }
	}
    }
  else
    return 0;
}
// void System::generate_H()
// {
//   int dim=configurations.size();
//   H.resize(dim,dim);
//   for(int i=0;i<dim;i++)
//     for(int j=i;j<dim;j++)
//       {
// 	H(i,j)=getH(i,j);
// 	H(j,i)=H(i,j);
//       }
// }

void System::MBPT()
{
    MBPT_calculateCoefficients();
    MBPT_calculateDeltaE();
}

void System::MBPT_calculateCoefficients()//calculate the coefficients for MBPT
{
    MBPT_coefficients.push_back(1.0); // C0
    MBPT_coefficients.push_back(MBPT_getCoefficient(2,3,4,5)); // C1
    MBPT_coefficients.push_back(MBPT_getCoefficient(2,3,6,7)); // C2
    MBPT_coefficients.push_back(MBPT_getCoefficient(0,1,4,5)); // C3
    MBPT_coefficients.push_back(MBPT_getCoefficient(0,1,6,7)); // C4
    // there is no C5 -- 4p4h excitation
}

void System::MBPT_printCoefficients()//printing coefficients for MBPT
{
    for (int i = 0; i < MBPT_coefficients.size(); i++)
        cout << setw(8) << MBPT_coefficients.at(i) << "\t";
    cout << endl;
}

double System::MBPT_getCoefficient(int i, int j, int a, int b)//get coefficient C_ij^ab
{
    double Ea = spOrbitals.SP_States.at(a).spEnergy;
    double Eb = spOrbitals.SP_States.at(b).spEnergy;
    double Ei = spOrbitals.SP_States.at(i).spEnergy;
    double Ej = spOrbitals.SP_States.at(j).spEnergy;
    return V2B(i,j,a,b) / (Ei+Ej-Ea-Eb);
}

void System::MBPT_calculateDeltaE()//calculate correlation energy for MBPT
{
    MBPT_deltaE = 0.0;
    for (int i = 0; i < A; i++)
    {
        for (int j = 0; j < A; j++)
        {
            for (int a = A; a < spOrbitals.SP_States.size(); a++)
            {
                for (int b = A; b < spOrbitals.SP_States.size(); b++)
                {
                    MBPT_deltaE += V2B(a,b,i,j) * MBPT_getCoefficient(i,j,a,b);
                }
            }
        }
    }
    MBPT_deltaE *= 0.25;
}

void System::generateConfigurations()
{
  int config = first(A);
  int configLast = last(spOrbitals.SP_States.size(),A);
  configurations.push_back(config);
  do
    {
      config = next(config);
      int NumUnpaired=0;
      for(int n=0;n<spOrbitals.SP_States.size()/2;n++)
	{
	  if(getBit(config,2*n)!=getBit(config,2*n+1))
	    ++NumUnpaired;
	}
      if(NumUnpaired==0)
	configurations.push_back(config);
    }
  while(config != configLast);
}

void System::printConfigurations()
{
    for (int i = 0; i < configurations.size(); i++)
        cout << bitset<8>(configurations.at(i)) << endl;
}

//////////////////////////////

int main()
{
    System system(4,4,0.4);//A, PMax, g

//    for (double g = -2.0; g <= 2.0; g += 0.1)
//    {
//        system.setG(g);
//        cout << g << "\t" << system.deltaE() << "\t";
//        system.diagonalization();
//    }

    system.generateConfigurations();
    system.printConfigurations();
    cout<<system.getH(0,2)<<endl;
    // system.MBPT();
    // system.MBPT_printCoefficients();
    // cout << system.MBPT_deltaE << endl;
    // system.diagonalization();
    // cout << system.diag_E_GS << endl;
    // system.diag_printCoefficients();
    return 0;
}
