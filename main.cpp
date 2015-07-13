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
    void GF_init();
    void GF_generateMatrices();
    void GF_diag();
    void GF_cal();
    void GF_iterate(int);
    double GF_E_GS;
    double GF_num;
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
    int Index_max_h;
    int Index_min_p;
    int dim_r,dim_q;
    double emin,emax;
    double sep_energy;
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
    {
      if((i-j)*(a-b)>0)
        return -0.5*g;
      else
	return 0.5*g;
    }
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
  //  MBPT_calculateCoefficients();
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
            for (int a = A; a < numberSP; a++)
            {
                for (int b = A; b < numberSP; b++)
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
    double factor = 1.0;//TODO
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
// void System::GF_init()
// {
//   GF_Z=MatrixXd::Zero(numberSP,numberSP);
//   GF_e.resize(numberSP,1);
//   for(int i=0;i<numberSP;i++)
//     {
//       GF_Z(i,i)=1;
//       GF_e(i)=spOrbitals.SP_States[i].spEnergy;
//     }
//   GF_Sigma_static.resize(numberSP,numberSP);
//   GF_Mstatic.resize(numberSP,numberSP);
//   Index_max_h=A-1;
//   Index_min_p=A;
//   emin=spOrbitals.SP_States[0].spEnergy+spOrbitals.SP_States[1].spEnergy-spOrbitals.SP_States[numberSP-1].spEnergy-2;
//   emax=spOrbitals.SP_States[numberSP-1].spEnergy+spOrbitals.SP_States[numberSP-2].spEnergy-spOrbitals.SP_States[0].spEnergy+2;
//   sep_energy=0.5*(spOrbitals.SP_States[A-1].spEnergy+spOrbitals.SP_States[A].spEnergy);

// }
// void System:: GF_generateMatrices()
// {
//   dim_r=0;
//   dim_q=0;
//   for(int k1=0;k1<=Index_max_h;k1++)
//     {
//       for(int k2=k1+1;k2<=Index_max_h;k2++)
// 	{
// 	  for(int n=Index_min_p;n<GF_e.rows();n++)
// 	    {
// 	      if((GF_e[k1]+GF_e[k2]-GF_e[n])>emin)
// 		++dim_r;
// 	    }
// 	}
//     }
//   for(int n1=Index_min_p;n1<GF_e.rows();n1++)
//     {
//       for(int n2=n1+1;n2<GF_e.rows();n2++)
// 	{
// 	  for(int k=0;k<=Index_max_h;k++)
// 	    {
// 	      if((GF_e[n1]+GF_e[n2]-GF_e[k])<emax)
// 		++dim_q;
// 	    }
// 	}
//     }
//   //  cout<<dim_r<<"\t"<<dim_q<<endl;
  
//   GF_Mr.resize(numberSP,dim_r);
//   GF_Mq.resize(numberSP,dim_q);
//   GF_Er.resize(dim_r,dim_r);
//   GF_Eq.resize(dim_q,dim_q);
//   GF_Mat.resize(numberSP+dim_r+dim_q,numberSP+dim_r+dim_q);

//   for(int i=0;i<numberSP;i++)
//     {
//       for(int j=i;j<numberSP;j++)
// 	{
// 	  GF_Sigma_static(i,j)=0;
// 	  for(int k=0;k<numberSP;k++)
// 	    {
// 	      for(int l=0;l<numberSP;l++)
// 		{
// 		  double rou_kl=0;
// 		  for(int index=0;index<=Index_max_h;index++)
// 		    rou_kl+=GF_Z(k,index)*GF_Z(l,index);
// 		  GF_Sigma_static(i,j)+=V2B(i,k,j,l)*rou_kl;
// 		}
// 	    }
// 	  GF_Sigma_static(j,i)=GF_Sigma_static(i,j);
// 	}
//     }
//   for(int i=0;i<numberSP;i++)
//     {
//       for(int j=i;j<numberSP;j++)
// 	{
// 	  GF_Mstatic(i,j)=GF_Sigma_static(i,j);
// 	  if(i==j)
// 	    {
// 		GF_Mstatic(i,i)+=spOrbitals.SP_States[i].spEnergy;
// 	    }
// 	  else
// 	    {
// 	      GF_Mstatic(j,i)=GF_Mstatic(i,j);
// 	    }
// 	}
//     }
//   GF_Er=MatrixXd::Zero(dim_r,dim_r);
//   GF_Eq=MatrixXd::Zero(dim_q,dim_q);

//   int index=0;
//   for(int k1=0;k1<=Index_max_h;k1++)
//     {
//       for(int k2=k1+1;k2<=Index_max_h;k2++)
//   	{
//   	  for(int n=Index_min_p;n<GF_e.rows();n++)
//   	    {
//   	      if((GF_e[k1]+GF_e[k2]-GF_e[n])>emin)
//   		{
//   		  GF_Er(index,index)=GF_e[k1]+GF_e[k2]-GF_e[n];
//   		  for(int a=0;a<numberSP;a++)
//   		    {
//   		      double Mr=0;
//   		      for(int b=0;b<numberSP;b++)
//   			{
//   			  for(int c=0;c<numberSP;c++)
//   			    {
//   			      for(int d=0;d<numberSP;d++)
//   				{
//   				  Mr+=V2B(a,b,c,d)*GF_Z(b,n)*GF_Z(c,k1)*GF_Z(d,k2);
//   				}
//   			    }
//   			}
// 		      GF_Mr(a,index)=Mr;
//   		    }
// 		  ++index;
//   		}
//   	    }
//   	}
//     }

//   index=0;
//   for(int n1=Index_min_p;n1<GF_e.size();n1++)
//     {
//       for(int n2=n1+1;n2<GF_e.rows();n2++)
//   	{
//   	  for(int k=0;k<=Index_max_h;k++)
//   	    {
//   	      if((GF_e[n1]+GF_e[n2]-GF_e[k])<emax)
//   		{
//   		  GF_Eq(index,index)=GF_e[n1]+GF_e[n2]-GF_e[k];
//   		  for(int a=0;a<numberSP;a++)
//   		    {
//   		      double Mq=0;
//   		      for(int b=0;b<numberSP;b++)
//   			{
//   			  for(int c=0;c<numberSP;c++)
//   			    {
//   			      for(int d=0;d<numberSP;d++)
//   				{
//   				  Mq+=V2B(a,b,c,d)*GF_Z(b,k)*GF_Z(c,n1)*GF_Z(d,n2);
//   				}
//   			    }
//   			}
//   		      GF_Mq(a,index)=Mq;
//   		    }
// 		  ++index;
//   		}
//   	    }
//   	}
//     }
  
//   GF_Mat.topLeftCorner(numberSP,numberSP)=GF_Mstatic;
//   GF_Mat.block(0,numberSP,numberSP,dim_r)=GF_Mr;
//   GF_Mat.topRightCorner(numberSP,dim_q)=GF_Mq;
//   GF_Mat.block(numberSP,0,dim_r,numberSP)=GF_Mr.adjoint();
//   GF_Mat.bottomLeftCorner(dim_q,numberSP)=GF_Mq.adjoint();
//   GF_Mat.block(numberSP,numberSP,dim_r,dim_r)=GF_Er;
//   GF_Mat.block(numberSP+dim_r,numberSP+dim_r,dim_q,dim_q)=GF_Eq;

//   // cout<<GF_Mstatic<<endl;
//   // cout<<"=========================\n";
//   // cout<<GF_Mr<<endl;
//   // cout<<"=========================\n";
//   // cout<<GF_Mq<<endl;
//   // cout<<"=========================\n";
//   // cout<<GF_Er<<endl;
//   // cout<<"=========================\n";
//   // cout<<GF_Eq<<endl;
//   // cout<<"=========================\n";
//   //  cout<<GF_Mat<<endl;
// }


void System::GF_init()
{
  GF_Z=MatrixXd::Zero(numberSP,numberSP);
  GF_e.resize(numberSP,1);
  for(int i=0;i<numberSP;i++)
    {
      GF_Z(i,i)=1;
      if(i<A)
	GF_e(i)=spOrbitals.SP_States[i].spEnergy-g/2;
      else
	GF_e(i)=spOrbitals.SP_States[i].spEnergy;
    }
  GF_Sigma_static.resize(numberSP,numberSP);
  GF_Mstatic.resize(numberSP,numberSP);
  Index_max_h=A-1;
  Index_min_p=A;
  emin=spOrbitals.SP_States[0].spEnergy+spOrbitals.SP_States[1].spEnergy-spOrbitals.SP_States[numberSP-1].spEnergy-1-g;
  emax=spOrbitals.SP_States[numberSP-1].spEnergy+spOrbitals.SP_States[numberSP-2].spEnergy-spOrbitals.SP_States[0].spEnergy+1+g/2;
  sep_energy=0.5*(spOrbitals.SP_States[A-1].spEnergy-g/2+spOrbitals.SP_States[A].spEnergy);

}
void System:: GF_generateMatrices()
{
  dim_r=0;
  dim_q=0;
  for(int k1=0;k1<=Index_max_h;k1++)
    {
      for(int k2=k1+1;k2<=Index_max_h;k2++)
	{
	  for(int n=Index_min_p;n<GF_e.rows();n++)
	    {
	      if((GF_e[k1]+GF_e[k2]-GF_e[n])>emin)
		++dim_r;
	    }
	}
    }
  for(int n1=Index_min_p;n1<GF_e.rows();n1++)
    {
      for(int n2=n1+1;n2<GF_e.rows();n2++)
	{
	  for(int k=0;k<=Index_max_h;k++)
	    {
	      if((GF_e[n1]+GF_e[n2]-GF_e[k])<emax)
		++dim_q;
	    }
	}
    }
  //  cout<<dim_r<<"\t"<<dim_q<<endl;
  
  GF_Mr.resize(numberSP,dim_r);
  GF_Mq.resize(numberSP,dim_q);
  GF_Er.resize(dim_r,dim_r);
  GF_Eq.resize(dim_q,dim_q);
  GF_Mat.resize(numberSP+dim_r+dim_q,numberSP+dim_r+dim_q);

  for(int i=0;i<numberSP;i++)
    {
      for(int j=i;j<numberSP;j++)
	{
	  GF_Sigma_static(i,j)=0;
	  for(int k=0;k<numberSP;k++)
	    {
	      for(int l=0;l<numberSP;l++)
		{
		  double rou_kl=0;
		  for(int index=0;index<=Index_max_h;index++)
		    rou_kl+=GF_Z(k,index)*GF_Z(l,index);
		  GF_Sigma_static(i,j)+=V2B(i,k,j,l)*rou_kl;
		}
	      if(k<A)
		GF_Sigma_static(i,j)-=V2B(i,k,j,k);
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
	      if(i<A)
		GF_Mstatic(i,i)+=spOrbitals.SP_States[i].spEnergy-g/2;
	      else
		GF_Mstatic(i,i)+=spOrbitals.SP_States[i].spEnergy;
	    }
	  else
	    {
	      GF_Mstatic(j,i)=GF_Mstatic(i,j);
	    }
	}
    }
  GF_Er=MatrixXd::Zero(dim_r,dim_r);
  GF_Eq=MatrixXd::Zero(dim_q,dim_q);

  int index=0;
  for(int k1=0;k1<=Index_max_h;k1++)
    {
      for(int k2=k1+1;k2<=Index_max_h;k2++)
  	{
  	  for(int n=Index_min_p;n<GF_e.rows();n++)
  	    {
  	      if((GF_e[k1]+GF_e[k2]-GF_e[n])>emin)
  		{
  		  GF_Er(index,index)=GF_e[k1]+GF_e[k2]-GF_e[n];
  		  for(int a=0;a<numberSP;a++)
  		    {
  		      double Mr=0;
  		      for(int b=0;b<numberSP;b++)
  			{
  			  for(int c=0;c<numberSP;c++)
  			    {
  			      for(int d=0;d<numberSP;d++)
  				{
  				  Mr+=V2B(a,b,c,d)*GF_Z(b,n)*GF_Z(c,k1)*GF_Z(d,k2);
  				}
  			    }
  			}
		      GF_Mr(a,index)=Mr;
  		    }
		  ++index;
  		}
  	    }
  	}
    }

  index=0;
  for(int n1=Index_min_p;n1<GF_e.size();n1++)
    {
      for(int n2=n1+1;n2<GF_e.rows();n2++)
  	{
  	  for(int k=0;k<=Index_max_h;k++)
  	    {
  	      if((GF_e[n1]+GF_e[n2]-GF_e[k])<emax)
  		{
  		  GF_Eq(index,index)=GF_e[n1]+GF_e[n2]-GF_e[k];
  		  for(int a=0;a<numberSP;a++)
  		    {
  		      double Mq=0;
  		      for(int b=0;b<numberSP;b++)
  			{
  			  for(int c=0;c<numberSP;c++)
  			    {
  			      for(int d=0;d<numberSP;d++)
  				{
  				  Mq+=V2B(a,b,c,d)*GF_Z(b,k)*GF_Z(c,n1)*GF_Z(d,n2);
  				}
  			    }
  			}
		      GF_Mq(a,index)=Mq;
  		    }
		  ++index;
  		}
  	    }
  	}
    }
  GF_Mat=MatrixXd::Zero(numberSP+dim_r+dim_q,numberSP+dim_r+dim_q);
  GF_Mat.topLeftCorner(numberSP,numberSP)=GF_Mstatic;
  GF_Mat.block(0,numberSP,numberSP,dim_r)=GF_Mr;
  GF_Mat.topRightCorner(numberSP,dim_q)=GF_Mq;
  GF_Mat.block(numberSP,0,dim_r,numberSP)=GF_Mr.adjoint();
  GF_Mat.bottomLeftCorner(dim_q,numberSP)=GF_Mq.adjoint();
  GF_Mat.block(numberSP,numberSP,dim_r,dim_r)=GF_Er;
  GF_Mat.block(numberSP+dim_r,numberSP+dim_r,dim_q,dim_q)=GF_Eq;

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
  //  cout<<GF_Mat<<endl;
}
void System::GF_diag()
{
  SelfAdjointEigenSolver<MatrixXd> solver(GF_Mat);
  if(solver.info()!=Success) abort();
  VectorXd GF_e_temp;
  MatrixXd GF_Z_temp;
  GF_e_temp=solver.eigenvalues();
  GF_Z_temp=solver.eigenvectors().topLeftCorner(numberSP,numberSP+dim_r+dim_q);
  GF_W=solver.eigenvectors().bottomLeftCorner(dim_r+dim_q,numberSP+dim_r+dim_q);
  int dim=0;
  for(int i=0;i<GF_Z_temp.cols();i++)
    {
      if(GF_Z_temp.col(i).squaredNorm()>1e-5)
	dim++;
    }
  GF_e.resize(dim,1);
  GF_Z.resize(numberSP,dim);
  int count=0;
  for(int i=0;i<GF_Z_temp.cols();i++)
    {
      if(GF_Z_temp.col(i).squaredNorm()>1e-5)
	{
	  GF_e(count)=GF_e_temp(i);
	  GF_Z.col(count)=GF_Z_temp.col(i);
	  count++;
	}
    }
  for(int i=0;i<GF_e.rows();i++)
    {
      if(GF_e[i]<sep_energy)
	{
	  Index_max_h=i;
	  Index_min_p=i+1;
	}
      else
	{
	  break;
	}
    }
}
void System::GF_cal()
{
  GF_E_GS=0;
  
  for(int j=0;j<numberSP;j++)
    {
      for(int k=0;k<=Index_max_h;k++)
  	{
	  GF_E_GS+=(spOrbitals.SP_States[j].spEnergy+GF_e[k])*GF_Z(j,k)*GF_Z(j,k);
  	}
    }
  GF_E_GS*=0.5;

  GF_num=0;
  for(int k=0;k<=Index_max_h;k++)
    {
      double Sh=0;
      for(int j=0;j<numberSP;j++)
	{
	  Sh+=GF_Z(j,k)*GF_Z(j,k);
	}
      //      cout<<GF_e[k]<<"\t"<<Sh<<endl;
      GF_num+=Sh;
    }
}
void System::GF_iterate(int max_num_iter)
{
  GF_init();
  double delta=0;
  int i=0;
  GF_E_GS=0;
  do
    {
      double E_GS_old=GF_E_GS;
      GF_generateMatrices();
      GF_diag();
      GF_cal();
      delta=abs(GF_E_GS-E_GS_old);
      //      cout<<GF_E_GS<<"\t"<<GF_num<<"\t"<<delta<<endl;
      ++i;
    }while(delta>1e-6&&i<max_num_iter);
}

//////////////////////////////

int main()
{
  System system(4,8,0.4);//A, numberSP, g
 //  cout << "g\t" << "MBPT\t" << "diag\t" << "CCD\t" << "GF" << endl;
  
  for (double g = -1.0; g <= 1.0; g += 0.1)
    {
      system.set_g(g);
      system.CCD_generateMatrices();
      system.CCD_calculateTau();
      system.MBPT();
      system.diagonalization();
      system.GF_iterate(1);

      cout << g << "\t" << system.MBPT_deltaE<< "\t";
      cout << system.diag_E_GS-2+system.g << "\t";
      cout << system.CCD_deltaE<< "\t";


      cout<<system.GF_E_GS-2+system.g<<"\t"<<system.GF_num<<endl;
    }

  // system.CCD_generateMatrices();
  // system.CCD_calculateTau();
  // cout << system.CCD_deltaE+2-system.g<<endl;
  
  // system.diagonalization();
  // cout << system.diag_E_GS<<endl;

  // system.MBPT();
  // cout<<system.MBPT_deltaE+2-system.g<<endl;
    
  // system.GF_iterate(10);
  // cout<<system.GF_E_GS<<"\t"<<system.GF_num<<endl;
  return 0;
}
