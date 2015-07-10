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
    PairingSP(int _p, bool _spin): p(_p),spin(_spin),spEnergy(_p-1){}
    int p; //from 1
    bool spin; //0 -- spin up, 1 -- spin down
    double spEnergy; //single particle energy
};

//////////////////////////////

class Pairing
{
public:
    vector<PairingSP> SP_States;
    void generateSP_States(int);
};

void Pairing::generateSP_States(int numberSP)//generate s.p. states (p orbits, spin up/down)
{
    for (int p = 1; p <= numberSP/2; p++)
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
    int numberSP; //number of single particle orbitals

    //pairing
    Pairing spOrbitals;
    void set_g(double);

    //configurations
    void generateConfigurations();
    void printConfigurations();
    vector<int> configurations;

    //system
    double g;
    MatrixXd H;
    double V2B(int,int,int,int);
    double V1B(int,int);
    double getH(int,int);
    void generateH();

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

    //CCD
    void CCD_generateMatrices();
    void CCD_calculateTau();
    MatrixXd CCD_V, CCD_V_tilde, CCD_e, CCD_Tau, CCD_t;
    double CCD_deltaE;

    //GF
    void GF_generateMatrices();
    void GF_diag();
    double GF_E_GS;
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

void System::set_g(double _g)//set pairing g value
{
    g = _g;
}

System::System(int _A, int _numberSP, double g) : A(_A),numberSP(_numberSP)//generate system with A particles
{
    set_g(g);
    spOrbitals.generateSP_States(numberSP);
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
        cout << setw(8) << diag_coefficients[i]/diag_coefficients[0] << "\t";
    cout << endl;
}

double System::V2B(int i,int j,int a,int b)//calculate 2body matrix element <ij|V_pairing|ab>
{
    if ((spOrbitals.SP_States[a].p == spOrbitals.SP_States[b].p) && //same 'p'
        (spOrbitals.SP_States[a].spin != spOrbitals.SP_States[b].spin) &&//different spin
        (spOrbitals.SP_States[i].p == spOrbitals.SP_States[j].p) &&
        (spOrbitals.SP_States[i].spin != spOrbitals.SP_States[j].spin))
        return -0.5*g;
    else
        return 0;
}
double System::V1B(int a,int b)
{
  if((spOrbitals.SP_States[a].p == spOrbitals.SP_States[b].p) && 
     (spOrbitals.SP_States[a].spin == spOrbitals.SP_States[b].spin))
    return spOrbitals.SP_States[a].p-1;
}

void System::MBPT()
{
    MBPT_calculateCoefficients();
    MBPT_calculateDeltaE();
}

void System::MBPT_calculateCoefficients()//calculate the coefficients for MBPT
{
    MBPT_coefficients.clear();
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
void System::generateConfigurations()
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

void System::printConfigurations()
{
    for (int i = 0; i < configurations.size(); i++)
        cout << bitset<8>(configurations.at(i)) << endl;
}

//////////////////////////////

void System::CCD_generateMatrices()
{
  int size1 = (numberSP-A)*(numberSP-A-1)/2;
  int size2 = A*(A-1)/2;;
  int index1;
  int index2;
  CCD_V.resize(size1,size2);
  index2 = 0;
  for (int i = 0; i < A; i++)
    for (int j = i+1; j < A; j++)
      {
	index1 = 0;
	for (int a = A; a < numberSP; a++)
	  for (int b = a+1; b < numberSP; b++)
	    {
	      CCD_V(index1,index2) = V2B(a,b,i,j);
	      index1++;
	    }
	index2++;
      }

  CCD_V_tilde.resize(size1,size1);
  index2 = 0;
  for (int c = A; c < numberSP; c++)
    for (int d = c+1; d < numberSP; d++)
      {
	index1 = 0;
	for (int a = A; a < numberSP; a++)
	  for (int b = a+1; b < numberSP; b++)
	    {
	      CCD_V_tilde(index1,index2) = V2B(a,b,c,d);
	      index1++;
	    }
	index2++;
      }

  CCD_e.resize(size1,size2);
  index2 = 0;
  for (int i = 0; i < A; i++)
    for (int j = i+1; j < A; j++)
      {
	index1 = 0;
	for (int c = A; c < numberSP; c++)
	  for (int d = c+1; d < numberSP; d++)
	    {
	      double Ec = spOrbitals.SP_States.at(c).spEnergy;
	      double Ed = spOrbitals.SP_States.at(d).spEnergy;
	      double Ei = spOrbitals.SP_States.at(i).spEnergy;
	      double Ej = spOrbitals.SP_States.at(j).spEnergy;
	      CCD_e(index1,index2) = 1/(Ei+Ej-Ec-Ed);
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
    CCD_Tau = CCD_V;
    MatrixXd TauHelp;
    //double factor = 0.5;//TODO -0.0498237
    //double factor = 0.25;//TODO -0.0481901
    //double factor = 1;//TODO -0.0534759
    //double factor = 0;//TODO -0.0466667
    double factor = 0.5;//TODO
    TauHelp.resize(size1,size2);
    int i = 0;
    do
    {
        i++;
        TauHelp = CCD_Tau;
        MatrixXd help = CCD_e.array() * TauHelp.array();
        CCD_Tau = CCD_V + factor * CCD_V_tilde * help;
    }
    while(abs(CCD_Tau.norm() - TauHelp.norm()) > 1e-6);
    CCD_t = CCD_e.array() * CCD_Tau.array();
    CCD_deltaE = (CCD_t * CCD_V.transpose()).trace();
    //cout << "CCD  deltaE " << CCD_deltaE << endl;
}


void System:: GF_generateMatrices()
{
  int num_p=numberSP-A;
  int num_h=A;
  int dim_2h1p=(num_h*(num_h-1))/2*num_p;
  int dim_2p1h=(num_p*(num_p-1))/2*num_h;
  GF_Sigma_static.resize(numberSP,numberSP);
  GF_Mstatic.resize(numberSP,numberSP);
  GF_Mr.resize(numberSP,dim_2h1p);
  GF_Mq.resize(numberSP,dim_2p1h);
  GF_Er.resize(dim_2h1p,dim_2h1p);
  GF_Eq.resize(dim_2p1h,dim_2p1h);
  GF_Mat.resize(numberSP+dim_2h1p+dim_2p1h,numberSP+dim_2h1p+dim_2p1h);

  for(int i=0;i<numberSP;i++)
    {
      for(int j=i;j<numberSP;j++)
	{
	  GF_Sigma_static(i,j)=0;
	  for(int k=0;k<A;k++)
	    {
	      GF_Sigma_static(i,j)+=V2B(i,k,j,k);
	    }
	  GF_Sigma_static(j,i)=GF_Sigma_static(i,j);
	}
    }
  
  for(int i=0;i<numberSP;i++)
    {
      for(int j=i;j<numberSP;j++)
	{
	  GF_Mstatic(i,j)=GF_Sigma_static(i,j);
	  if(i==j)
	    {
	      GF_Mstatic(i,i)+=spOrbitals.SP_States[i].spEnergy;
	    }
	  else
	    {
	      GF_Mstatic(j,i)=GF_Mstatic(i,j);
	    }
	}
    }
  GF_Er=MatrixXd::Zero(dim_2h1p,dim_2h1p);
  GF_Eq=MatrixXd::Zero(dim_2p1h,dim_2p1h);

  int j=0;
  for(int n=A;n<numberSP;n++)
    {
      for(int k1=0;k1<A;k1++)
	for(int k2=k1+1;k2<A;k2++)
	  {
	    GF_Er(j,j)=spOrbitals.SP_States[k1].spEnergy+spOrbitals.SP_States[k2].spEnergy-spOrbitals.SP_States[n].spEnergy;
	    for(int i=0;i<numberSP;i++)
	      {
		GF_Mr(i,j)=V2B(i,n,k1,k2);
	      }
	    ++j;
	  }
    }
  
  j=0;  
  for(int k=0;k<A;k++)
    {
      for(int n1=A;n1<numberSP;n1++)
	for(int n2=n1+1;n2<numberSP;n2++)
	  {
	    GF_Eq(j,j)=spOrbitals.SP_States[n1].spEnergy+spOrbitals.SP_States[n2].spEnergy-spOrbitals.SP_States[k].spEnergy;
	    for(int i=0;i<numberSP;i++)
	      {
		GF_Mq(i,j)=V2B(i,k,n1,n2);
	      }
	    ++j;
	  }
    }
  GF_Mat.topLeftCorner(numberSP,numberSP)=GF_Mstatic;
  GF_Mat.block(0,numberSP,numberSP,dim_2h1p)=GF_Mr;
  GF_Mat.topRightCorner(numberSP,dim_2p1h)=GF_Mq;
  GF_Mat.block(numberSP,0,dim_2h1p,numberSP)=GF_Mr.adjoint();
  GF_Mat.bottomLeftCorner(dim_2p1h,numberSP)=GF_Mq.adjoint();
  GF_Mat.block(numberSP,numberSP,dim_2h1p,dim_2h1p)=GF_Er;
  GF_Mat.block(numberSP+dim_2h1p,numberSP+dim_2h1p,dim_2p1h,dim_2p1h)=GF_Eq;

  // cout<<GF_Mstatic<<endl;
  // cout<<"=========================\n";
  // cout<<GF_Mr<<endl;
  // cout<<"=========================\n";
  // cout<<GF_Mq<<endl;
  // cout<<"=========================\n";
  // cout<<GF_Er<<endl;
  // cout<<"=========================\n";
  // cout<<GF_Eq<<endl;
  // cout<<"=========================\n";
  // cout<<GF_Mat<<endl;
}
void System::GF_diag()
{
  int num_p=numberSP-A;
  int num_h=A;
  int dim_2h1p=(num_h*(num_h-1))/2*num_p;
  int dim_2p1h=(num_p*(num_p-1))/2*num_h;
  SelfAdjointEigenSolver<MatrixXd> solver(GF_Mat);
  if(solver.info()!=Success) abort();
  GF_e=solver.eigenvalues();
  GF_Z=solver.eigenvectors().topLeftCorner(numberSP,numberSP+dim_2p1h+dim_2h1p);
  GF_W=solver.eigenvectors().bottomLeftCorner(dim_2p1h+dim_2h1p,numberSP+dim_2p1h+dim_2h1p);

  double sep_energy=(spOrbitals.SP_States[A-1].spEnergy+spOrbitals.SP_States[A].spEnergy)/2;
  GF_E_GS=0;
  
  for(int j=0;j<numberSP;j++)
    {
      for(int k=0;k<numberSP+dim_2p1h+dim_2h1p;k++)
  	{
  	  if(GF_e[k]<sep_energy)
  	    {
  	      GF_E_GS+=(spOrbitals.SP_States[j].spEnergy+GF_e[k])*GF_Z(j,k)*GF_Z(j,k);
  	    }
  	}
    }
  GF_E_GS*=0.5;
}

//////////////////////////////

int main()
{
  System system(4,8,0.4);//A, numberSP, g

  cout << "g\t" << "MBPT\t" << "diag\t" << "CCD" << "GF" << endl;

  for (double g = -2.0; g <= 2.0; g += 0.1)
    {
      system.set_g(g);
      system.CCD_generateMatrices();
      system.CCD_calculateTau();
      system.MBPT();
      
      system.diagonalization();

      system.GF_generateMatrices();
      system.GF_diag();

      cout << g << "\t" << system.MBPT_deltaE<< "\t";
      cout << system.diag_E_GS-2+system.g << "\t";
      cout << system.CCD_deltaE<< "\t";
      cout <<system.GF_E_GS-2+system.g<<endl;
    }
 
  return 0;
}
