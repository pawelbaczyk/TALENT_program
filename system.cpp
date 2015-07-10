#include "system.h"

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
        cout << MBPT_coefficients.at(i) << "\t";
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
                GF_Mstatic(i,i)+=SP_States[i]->spEnergy;
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
                GF_Er(j,j)=SP_States[k1]->spEnergy+SP_States[k2]->spEnergy-SP_States[n]->spEnergy;
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
                GF_Eq(j,j)=SP_States[n1]->spEnergy+SP_States[n2]->spEnergy-SP_States[k]->spEnergy;
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

    double sep_energy=(SP_States[A-1]->spEnergy+SP_States[A]->spEnergy)/2;
    GF_E_GS=0;

    for(int j=0;j<numberSP;j++)
    {
        for(int k=0;k<numberSP+dim_2p1h+dim_2h1p;k++)
        {
            if(GF_e[k]<sep_energy)
            {
                GF_E_GS+=(SP_States[j]->spEnergy+GF_e[k])*GF_Z(j,k)*GF_Z(j,k);
            }
        }
    }
    GF_E_GS*=0.5;
}

//////////////////////////////
