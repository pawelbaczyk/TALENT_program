#include "system.h"

System::System(int _A) : A(_A)//generate system with A particles
{
}

System::~System()
{
    for (int i = SP_States.size()-1; i >= 0; i--)
    {
        delete SP_States[i];
    }
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
    double Ea = V1B(a,a);
    double Eb = V1B(b,b);
    double Ei = V1B(i,i);
    double Ej = V1B(j,j);
    return V2B(i,j,a,b) / (Ei+Ej-Ea-Eb);
}

void System::MBPT_calculateDeltaE()//calculate correlation energy for MBPT
{
    MBPT_deltaE = 0.0;
    for (int i = 0; i < A; i++)
    {
        for (int j = i+1; j < A; j++)
        {
            for (int a = A; a < SP_States.size(); a++)
            {
                for (int b = a+1; b < SP_States.size(); b++)
                {
                    MBPT_deltaE += V2B(a,b,i,j) * MBPT_getCoefficient(i,j,a,b);
                }
            }
        }
    }
}

void System::CCD_OnFlight()
{
    MatrixXd CCD_V_ph, CCD_V_pp, CCD_V_hh, CCD_e_ph, CCD_Tau;

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
                    double Ec = V1B(c,c);
                    double Ed = V1B(d,d);
                    double Ei = V1B(i,i);
                    double Ej = V1B(j,j);
                    CCD_e_ph(index1,index2) = 1/(Ei+Ej-Ec-Ed);
                    index1++;
                }
            index2++;
        }


    CCD_Tau.resize(size_p,size_h);
    CCD_Tau = CCD_V_ph;
    MatrixXd TauIter;
    TauIter.resize(size_p,size_h);
    CCD_OnFlight_t_m.resize(size_p,size_h);
    int iter = 0;
    double CCD_deltaE_iter;
    do
    {
        iter++;
        TauIter = CCD_Tau;
        CCD_OnFlight_t_m = CCD_e_ph.array() * TauIter.array();
        CCD_deltaE_iter = (CCD_OnFlight_t_m * CCD_V_ph.transpose()).trace();
//        CCD_Tau = CCD_V_ph
//                + CCD_V_pp * CCD_t_m;
//               // + CCD_t_m * CCD_V_hh;
        int index1;
        int index2;
        index2 = 0;
        for (int i = 0; i < A; i++)
            for (int j = i+1; j < A; j++)
            {
                index1 = 0;
                for (int a = A; a < numberSP; a++)
                    for (int b = a+1; b < numberSP; b++)
                    {
                        double help = 0.0;
                        CCD_Tau(index1,index2) = V2B(a,b,i,j);
                        for (int c = A; c < numberSP; c++)
                            for (int d = A; d < numberSP; d++)
                                help += 0.5 * V2B(a,b,c,d) * CCD_OnFlight_t(i,j,c,d);
                        for (int k = 0; k < A; k++)
                            for (int l = 0; l < A; l++)
                                help += 0.5 * V2B(k,l,i,j) * CCD_OnFlight_t(k,l,a,b);
                        for (int k = 0; k < A; k++)
                            for (int c = A; c < numberSP; c++)
                                help += V2B(k,b,c,j) * CCD_OnFlight_t(i,k,a,c)
                                      - V2B(k,b,c,i) * CCD_OnFlight_t(j,k,a,c)
                                      - V2B(k,a,c,j) * CCD_OnFlight_t(i,k,b,c)
                                      + V2B(k,a,c,i) * CCD_OnFlight_t(j,k,b,c);
                        for (int k = 0; k < A; k++)
                            for (int l = 0; l < A; l++)
                                for (int c = A; c < numberSP; c++)
                                    for (int d = A; d < numberSP; d++)
                                    {
                                        help += 0.25 * V2B(k,l,c,d) * CCD_OnFlight_t(i,j,c,d) * CCD_OnFlight_t(k,l,a,b);
                                        help += V2B(k,l,c,d) * CCD_OnFlight_t(i,k,a,c) * CCD_OnFlight_t(j,l,b,d)
                                              - V2B(k,l,c,d) * CCD_OnFlight_t(j,k,a,c) * CCD_OnFlight_t(i,l,b,d);
                                        help += -0.5 * V2B(k,l,c,d) * CCD_OnFlight_t(i,k,d,c) * CCD_OnFlight_t(l,j,a,b)
                                              +  0.5 * V2B(k,l,c,d) * CCD_OnFlight_t(j,k,d,c) * CCD_OnFlight_t(l,i,a,b);
                                        help += -0.5 * V2B(k,l,c,d) * CCD_OnFlight_t(l,k,a,c) * CCD_OnFlight_t(i,j,d,b)
                                              +  0.5 * V2B(k,l,c,d) * CCD_OnFlight_t(l,k,b,c) * CCD_OnFlight_t(i,j,d,a);
                                    }

                        CCD_Tau(index1,index2) += help;
                        index1++;
                    }
                index2++;
            }
        CCD_OnFlight_t_m = CCD_e_ph.array() * CCD_Tau.array();
        CCD_deltaE = (CCD_OnFlight_t_m * CCD_V_ph.transpose()).trace();
    }
    while (abs(CCD_deltaE_iter - CCD_deltaE) > 1e-8);
    CCD_OnFlight_t_m = CCD_e_ph.array() * CCD_Tau.array();
    CCD_deltaE = (CCD_OnFlight_t_m * CCD_V_ph.transpose()).trace();
}

double System::CCD_OnFlight_t(int i,int j,int a,int b)
{
    if ((i == j) || (a == b))
        return 0;
    int h1 = min(i,j);
    int h2 = max(i,j);
    int h_sign = 2*(int)(i<j)-1;
    int p1 = min(a,b) - A;
    int p2 = max(a,b) - A;
    int p_sign = 2*(int)(a<b)-1;
    int nh = (2*A-h1)*(h1+1)/2 - A + h2 - h1 - 1;
    int B = numberSP - A;
    int np = (2*B-p1)*(p1+1)/2 - B + p2 - p1 - 1;
    return h_sign * p_sign * CCD_OnFlight_t_m(np,nh);//TODO check indices
}

void System::CCD_SparseMatrices()
{
    SparseMatrix<double> CCD_V_ph, CCD_V_pp, CCD_V_hh, CCD_e_ph, CCD_Tau, CCD_t_m;

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(1e8);//TODO release!

    int size_p = (numberSP-A)*(numberSP-A-1)/2;
    int size_h = A*(A-1)/2;
    cerr << size_p*size_p << endl;
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
                    //CCD_V_ph(index1,index2) = V2B(a,b,i,j);
                    double V = V2B(a,b,i,j);
                    if (V != 0.0)//TODO
                        tripletList.push_back(T(index1,index2,V));
                    index1++;
                }
            index2++;
        }
    CCD_V_ph.setFromTriplets(tripletList.begin(), tripletList.end());

    tripletList.clear();
    CCD_V_pp.resize(size_p,size_p);
    index2 = 0;
    for (int c = A; c < numberSP; c++)
        for (int d = c+1; d < numberSP; d++)
        {
            index1 = 0;
            for (int a = A; a < numberSP; a++)
                for (int b = a+1; b < numberSP; b++)
                {
                    //CCD_V_pp(index1,index2) = V2B(a,b,c,d);
                    double V = V2B(a,b,c,d);
                    if (V != 0.0)//TODO
                        tripletList.push_back(T(index1,index2,V));
                    index1++;
                }
            index2++;
        }
    CCD_V_pp.setFromTriplets(tripletList.begin(), tripletList.end());

    tripletList.clear();
    CCD_V_hh.resize(size_h,size_h);
    index2 = 0;
    for (int i = 0; i < A; i++)
        for (int j = i+1; j < A; j++)
        {
            index1 = 0;
            for (int k = 0; k < A; k++)
                for (int l = k+1; l < A; l++)
                {
                    //CCD_V_hh(index1,index2) = V2B(k,l,i,j);
                    double V = V2B(k,l,i,j);
                    if (V != 0.0)//TODO
                        tripletList.push_back(T(index1,index2,V));
                    index1++;
                }
            index2++;
        }
    CCD_V_hh.setFromTriplets(tripletList.begin(), tripletList.end());

    tripletList.clear();
    CCD_e_ph.resize(size_p,size_h);
    index2 = 0;
    for (int i = 0; i < A; i++)
        for (int j = i+1; j < A; j++)
        {
            index1 = 0;
            for (int c = A; c < numberSP; c++)
                for (int d = c+1; d < numberSP; d++)
                {
                    double Ec = V1B(c,c);
                    double Ed = V1B(d,d);
                    double Ei = V1B(i,i);
                    double Ej = V1B(j,j);
                    //CCD_e_ph(index1,index2) = 1/(Ei+Ej-Ec-Ed);
                    tripletList.push_back(T(index1,index2,1/(Ei+Ej-Ec-Ed)));
                    index1++;
                }
            index2++;
        }
    CCD_e_ph.setFromTriplets(tripletList.begin(), tripletList.end());

    //calculating

    CCD_Tau.resize(size_p,size_h);
    CCD_Tau = CCD_V_ph;
    SparseMatrix<double> TauIter, help;
    TauIter.resize(size_p,size_h);
    help.resize(size_p,size_p);
    CCD_t_m.resize(size_p,size_h);
    int iter = 0;
    double CCD_deltaE_iter;
    do
    {
        iter++;
        TauIter = CCD_Tau;
        CCD_t_m = CCD_e_ph.cwiseProduct(TauIter);
        help = CCD_t_m * SparseMatrix<double>(CCD_V_ph.transpose());
        double sum =0;
        for (int k=0; k < help.outerSize(); ++k)
            sum += help.coeff(k,k);
        CCD_deltaE_iter = sum;
        CCD_Tau = CCD_V_ph
                + CCD_V_pp * CCD_t_m;
                //+ CCD_t_m * CCD_V_hh;
        cerr << iter << "\t" <<CCD_deltaE_iter<< endl;

        CCD_t_m = CCD_e_ph.cwiseProduct(CCD_Tau);
        help = CCD_t_m * SparseMatrix<double>(CCD_V_ph.transpose());
         sum =0;
        for (int k=0; k < help.outerSize(); ++k)
            sum += help.coeff(k,k);
        CCD_deltaE = sum;
    }
    while (abs(CCD_deltaE_iter - CCD_deltaE) > 1e-8);
    help = CCD_t_m * SparseMatrix<double>(CCD_V_ph.transpose());
    double sum =0;
    for (int k=0; k < help.outerSize(); ++k)
        sum += help.coeff(k,k);
    CCD_deltaE = sum;
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
                GF_Mstatic(i,i)+=V1B(i,i);
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
                GF_Er(j,j)=V1B(k1,k1)+V1B(k2,k2)-V1B(n,n);
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
                GF_Eq(j,j)=V1B(n1,n1)+V1B(n2,n2)-V1B(k,k);
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

    double sep_energy=(V1B(A-1,A-1)+V1B(A,A))/2;
    GF_E_GS=0;

    for(int j=0;j<numberSP;j++)
    {
        for(int k=0;k<numberSP+dim_2p1h+dim_2h1p;k++)
        {
            if(GF_e[k]<sep_energy)
            {
                GF_E_GS+=(V1B(j,j)+GF_e[k])*GF_Z(j,k)*GF_Z(j,k);
            }
        }
    }
    GF_E_GS*=0.5;
}

//////////////////////////////
