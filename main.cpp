#include <iostream>
#include <vector>
#include <cmath>
#include <bitset>
//#include "help.h"
#include "Eigen/Dense"
#include <iomanip>

using namespace std;
using namespace Eigen;

//////////////////////////////







int first(int A)
{
    int result = 0;
    for (int i = 0; i < A; i++)
    {
        result += pow(2,i);
    }
    return result;
}

int last(int n, int A)
{
    int result = 0;
    for (int i = 0; i < A; i++)
    {
        result += pow(2,n-1-i);
    }
    return result;
}

int next(int x)
{
   int smallest, ripple, new_smallest, ones;

   if (x == 0) return 0;
   smallest     = (x & -x);
   ripple       = x + smallest;
   new_smallest = (ripple & -ripple);
   ones         = ((new_smallest/smallest) >> 1) - 1;
   return ripple | ones;
}

bool getBit(int x,int n)//get the nth bit of integer x
{
  return ((1<<n)&x)>>n;
}










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
    VectorXd diag_coefficients;

    //MBPT
    void MBPT();
    double MBPT_getCoefficient(int,int,int,int);
    void MBPT_printCoefficients();
    void MBPT_calculateDeltaE();
    vector<double> MBPT_coefficients;
    double MBPT_deltaE;

    //CCD
    void CCD_generateMatrices();
    void CCD_calculateTau();
    MatrixXd CCD_V_ph, CCD_V_pp, CCD_V_hh, CCD_e_ph, CCD_Tau, CCD_t;
    double CCD_deltaE;
};

System::System(int _A, int _numberSP) : A(_A), numberSP(_numberSP)//generate system with A particles
{
}

void System::generateH()//calculate the hamiltonian matrix
{
  generateConfigurations();
  int dim=configurations.size();
  H.resize(dim,dim);
  for(int i=0;i<dim;i++)
    for(int j=i;j<dim;j++)
      {
    H(i,j)=getH(i,j);
    H(j,i)=H(i,j);
      }
}

double System::getH(int a,int b)//return the matrix element between configuration a and b
{
  int bra=configurations[a];
  int ket=configurations[b];
  int diffCount=0;
  vector<int> OccupiedPos;
  vector<int> UnoccupiedPos;
  vector<int> NumOrbitalsKet;
  vector<int> NumOrbitalsBra;
  int numKet(0),numBra(0);
  for(int i=0;i<numberSP;i++)
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
      double temp=0;
      for(int j=0;j<numberSP;j++)
    {
      if(getBit(ket,j))
        {
          temp+=V2B(UnoccupiedPos[0],j,OccupiedPos[0],j);
        }
    }
      return phase*(V1B(UnoccupiedPos[0],OccupiedPos[0])+temp);
    }
  else if(diffCount==0)
    {
      double temp=0;
      for(int i=0;i<numberSP;i++)
    {
      if(getBit(ket,i))
        {
          temp+=V1B(i,i);
          for(int j=0;j<numberSP;j++)
        {
          if(getBit(ket,j))
            {
              temp+=0.5*V2B(i,j,i,j);
            }
        }
        }
    }
      return temp;
    }
  else
    return 0;
}

void System::diagonalization()//perform diagonalization
{
  generateH();
  SelfAdjointEigenSolver<MatrixXd> solver(H);
  if(solver.info()!=Success) abort();
  diag_E_GS=solver.eigenvalues()[0];
  diag_coefficients=solver.eigenvectors().col(0);
}

void System::diag_printCoefficients()//print ground state eigenvector with C_0 = 1
{
    for (int i = 0; i < diag_coefficients.size(); i++)
        cout << diag_coefficients[i]/diag_coefficients[0] << "\t";
    cout << endl;
}

void System::MBPT()
{
    MBPT_calculateDeltaE();
}

void System::MBPT_printCoefficients()//printing coefficients for MBPT
{
    for (int i = 0; i < MBPT_coefficients.size(); i++)
        cout << setw(8) << MBPT_coefficients.at(i) << "\t";
    cout << endl;
}

double System::MBPT_getCoefficient(int i, int j, int a, int b)//get coefficient C_ij^ab
{
    double Ea = SP_States.at(a)->spEnergy;
    double Eb = SP_States.at(b)->spEnergy;
    double Ei = SP_States.at(i)->spEnergy;
    double Ej = SP_States.at(j)->spEnergy;
    return V2B(i,j,a,b) / (Ei+Ej-Ea-Eb);
}

void System::MBPT_calculateDeltaE()//calculate correlation energy for MBPT
{
    MBPT_deltaE = 0.0;
    for (int i = 0; i < A; i++)
    {
        for (int j = 0; j < A; j++)
        {
            for (int a = A; a < SP_States.size(); a++)
            {
                for (int b = A; b < SP_States.size(); b++)
                {
                    MBPT_deltaE += V2B(a,b,i,j) * MBPT_getCoefficient(i,j,a,b);
                }
            }
        }
    }
    MBPT_deltaE *= 0.25;
}

void System::CCD_generateMatrices()
{
  int size_p = (numberSP-A)*(numberSP-A-1)/2;
  int size_h = A*(A-1)/2;;
  int index1;
  int index2;
  CCD_V_ph.resize(size_p,size_h);
  index2 = 0;
  for (int i = 0; i < A; i++)
    for (int j = i+1; j < A; j++)
      {
    index1 = 0;
    for (int a = A; a < numberSP; a++)
      for (int b = a+1; b < numberSP; b++)
        {
          CCD_V_ph(index1,index2) = V2B(a,b,i,j);
          index1++;
        }
    index2++;
      }

  CCD_V_pp.resize(size_p,size_p);
  index2 = 0;
  for (int c = A; c < numberSP; c++)
    for (int d = c+1; d < numberSP; d++)
      {
        index1 = 0;
    for (int a = A; a < numberSP; a++)
      for (int b = a+1; b < numberSP; b++)
        {
          CCD_V_pp(index1,index2) = V2B(a,b,c,d);
          index1++;
        }
    index2++;
      }

  CCD_V_hh.resize(size_h,size_h);
  index2 = 0;
  for (int i = 0; i < A; i++)
    for (int j = i+1; j < A; j++)
      {
        index1 = 0;
        for (int k = 0; k < A; k++)
          for (int l = k+1; l < A; l++)
        {
          CCD_V_hh(index1,index2) = V2B(k,l,i,j);
          index1++;
        }
    index2++;
      }

  CCD_e_ph.resize(size_p,size_h);
  index2 = 0;
  for (int i = 0; i < A; i++)
    for (int j = i+1; j < A; j++)
      {
        index1 = 0;
    for (int c = A; c < numberSP; c++)
      for (int d = c+1; d < numberSP; d++)
        {
          double Ec = SP_States.at(c)->spEnergy;
          double Ed = SP_States.at(d)->spEnergy;
          double Ei = SP_States.at(i)->spEnergy;
          double Ej = SP_States.at(j)->spEnergy;
          CCD_e_ph(index1,index2) = 1/(Ei+Ej-Ec-Ed);
          index1++;
        }
    index2++;
      }

}

void System::CCD_calculateTau()
{
    int size1 = (numberSP-A)*(numberSP-A-1)/2;
    int size2 = A*(A-1)/2;
    CCD_Tau.resize(size1,size2);
    CCD_Tau = CCD_V_ph;
    MatrixXd TauHelp;
    //double factor = 0.5;//TODO -0.0498237
    //double factor = 0.25;//TODO -0.0481901
    //double factor = 1;//TODO -0.0534759
    //double factor = 0;//TODO -0.0466667
    double factor = 1;//TODO
    TauHelp.resize(size1,size2);
    int i = 0;
    do
    {
        i++;
        TauHelp = CCD_Tau;
        MatrixXd help = CCD_e_ph.array() * TauHelp.array();
        CCD_Tau = CCD_V_ph + factor * CCD_V_pp * help;
        //cout << endl << i << endl << CCD_Tau << endl;
    }
    while(abs(CCD_Tau.norm() - TauHelp.norm()) > 1e-6);
    i = 0;
    do
    {
        i++;
        TauHelp = CCD_Tau;
        MatrixXd help = CCD_e_ph.array() * TauHelp.array();
        CCD_Tau = CCD_V_ph + factor * CCD_V_pp * help + factor * help * CCD_V_hh;
    }
    while(abs(CCD_Tau.norm() - TauHelp.norm()) > 1e-6);
    CCD_t = CCD_e_ph.array() * CCD_Tau.array();
    CCD_deltaE = (CCD_t * CCD_V_ph.transpose()).trace();
    //cout << "CCD  deltaE " << CCD_deltaE << endl;
}

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

Pairing::Pairing(int _A, int _numberSP, double _g) : System(_A, _numberSP)
{
    g =_g;
    generateSP_States(numberSP);
}

void Pairing::generateSP_States(int numberSP)//generate s.p. states (p orbits, spin up/down)
{
    for (int p = 1; p <= numberSP/2; p++)
    {
        SP_States.push_back(new SP_Pairing(p,0));//spin up
        SP_States.push_back(new SP_Pairing(p,1));//spin down
    }
}

void Pairing::generateConfigurations()
{
    int config = first(A);
    int configLast = last(numberSP,A);
    configurations.clear();
    configurations.push_back(config);
    do
    {
        config = next(config);
        int NumUnpaired=0;
        for(int n=0;n<numberSP/2;n++)
        {
            if(getBit(config,2*n)!=getBit(config,2*n+1))
            ++NumUnpaired;
        }
        if(NumUnpaired==0)
            configurations.push_back(config);
    }
    while(config != configLast);
}

void Pairing::printConfigurations()
{
    for (int i = 0; i < configurations.size(); i++)
        cout << bitset<8>(configurations.at(i)) << endl;//TODO 8
}

double Pairing::V1B(int a,int b)
{
    if((((SP_Pairing*)SP_States[a])->p == ((SP_Pairing*)SP_States[b])->p) &&
        (((SP_Pairing*)SP_States[a])->spin == ((SP_Pairing*)SP_States[b])->spin))
        return ((SP_Pairing*)SP_States[a])->p-1;
}

double Pairing::V2B(int i,int j,int a,int b)//calculate 2body matrix element <ij|V_pairing|ab>
{
    if ((((SP_Pairing*)SP_States[a])->p == ((SP_Pairing*)SP_States[b])->p) && //same 'p'
        (((SP_Pairing*)SP_States[a])->spin != ((SP_Pairing*)SP_States[b])->spin) &&//different spin
        (((SP_Pairing*)SP_States[i])->p == ((SP_Pairing*)SP_States[j])->p) &&
        (((SP_Pairing*)SP_States[i])->spin != ((SP_Pairing*)SP_States[j])->spin))
        return -0.5*g;
    else
        return 0;
}

int main()
{
  Pairing pairing(4,8,0.4);//A, numberSP, g

  cout << "g\t" << "MBPT\t" << "diag\t" << "CCD" << endl;

  for (double g = -1.0; g <= 1.0; g += 0.1)
    {
      pairing.g = g;
      pairing.CCD_generateMatrices();
      pairing.CCD_calculateTau();
      pairing.MBPT();

      pairing.diagonalization();

      cout << g << "\t" << pairing.MBPT_deltaE+2-pairing.g<< "\t";
      cout << pairing.diag_E_GS << "\t";
      cout << pairing.CCD_deltaE +2-pairing.g<< endl;
    }
  return 0;
}
