#include "infinite.h"


bool TwoBody_compare::operator() (TwoBody_State *P1, TwoBody_State *P2)
{
    TwoBody_Infinite* TwoBody1 = (TwoBody_Infinite*)P1;
    TwoBody_Infinite* TwoBody2 = (TwoBody_Infinite*)P2;
    if(TwoBody1->Sz != TwoBody2->Sz)
        return TwoBody1->Sz < TwoBody2->Sz;
    else if(TwoBody1->Tz != TwoBody2->Tz)
        return TwoBody1->Tz < TwoBody2->Tz;
    else if(TwoBody1->Nx != TwoBody2->Nx)
        return TwoBody1->Nx < TwoBody2->Nx;
    else if(TwoBody1->Ny != TwoBody2->Ny)
        return TwoBody1->Ny < TwoBody2->Ny;
    else if(TwoBody1->Nz != TwoBody2->Nz)
        return TwoBody1->Nz < TwoBody2->Nz;
    return false;
}

bool Channel::Include(TwoBody_Infinite *TwoBody)
{
    if((TwoBody->Nx == Nx) && (TwoBody->Ny == Ny) && (TwoBody->Nz == Nz) &&
            (TwoBody->Sz == Sz) && (TwoBody->Tz==Tz))
        return true;
    return false;
}

Infinite::Infinite(int _A, int _g_s, double _rho, int _nMax) : System(_A),gauss_x(100),gauss_w(100)
{
    g_s = _g_s;
    nMax = _nMax;
    setRho(_rho);
    gauleg(-1,1,gauss_x,gauss_w);
}

void Infinite::generateSP_States(int A)
{
    SP_States.clear();
    k_F = pow(6*M_PI*M_PI*rho/g_s,1./3);
    L = pow(A/rho,1./3);
    for (int N = 0; N <= 3*nMax*nMax; N++)
        for (int nX = -nMax; nX <= nMax; nX++)
            for (int nY = -nMax; nY <= nMax; nY++)
                for (int nZ = -nMax; nZ <= nMax; nZ++)
                    if (nX*nX + nY*nY + nZ*nZ == N)
                    {
                        double kx = 2 * M_PI * nX / L;
                        double ky = 2 * M_PI * nY / L;
                        double kz = 2 * M_PI * nZ / L;
                        if (g_s == 2 || g_s == 4)
                        {
                            SP_States.push_back(new SP_Infinite(nX,nY,nZ,kx,ky,kz,0,0));
                            SP_States.push_back(new SP_Infinite(nX,nY,nZ,kx,ky,kz,1,0));
                        }
                        if (g_s == 4)
                        {
                            SP_States.push_back(new SP_Infinite(nX,nY,nZ,kx,ky,kz,0,1));
                            SP_States.push_back(new SP_Infinite(nX,nY,nZ,kx,ky,kz,1,1));
                        }
                    }
    numberSP = SP_States.size();
    //HF_calculateE0;
}

void Infinite::setRho(double _rho)
{
    rho = _rho;
    generateSP_States(A);
    //    cout << "Generating matrices" << endl;
    //    map_generateV2B();
    //    cout << "Done" << endl;
}

void Infinite::printSP_States()
{
    for (int i = 0; i < SP_States.size(); i++)
    {
        SP_Infinite* state = (SP_Infinite*)SP_States[i];
        cout << std::fixed << i+1 << "\t" << state->spEnergy << "\t"
             << state->kx << "\t" << state->ky << "\t" << state->kz << "\t"
             << (state->spin == 0 ? "+" : "-") << "\t" << (state->isospin == 0 ? "n" : "p")
             << endl;
    }
}



void Infinite::generateTwoBody_States()
{
    for(int i=0;i<A;i++)
        for(int j=i+1;j<A;j++)
        {
            SP_Infinite* I=(SP_Infinite*)SP_States[i];
            SP_Infinite* J=(SP_Infinite*)SP_States[j];
            int Nx=I->nx + J->nx;
            int Ny=I->ny + J->ny;
            int Nz=I->nz + J->nz;
            int Sz=I->spin + J->spin;
            int Tz=I->isospin + J->isospin;
            TwoBody_States_hh.push_back(new TwoBody_Infinite(i,j,Nx,Ny,Nz,Sz,Tz));
        }
    for(int i=A;i<numberSP;i++)
        for(int j=0;j<A;j++)
        {
            SP_Infinite* I=(SP_Infinite*)SP_States[i];
            SP_Infinite* J=(SP_Infinite*)SP_States[j];
            int Nx=I->nx + J->nx;
            int Ny=I->ny + J->ny;
            int Nz=I->nz + J->nz;
            int Sz=I->spin + J->spin;
            int Tz=I->isospin + J->isospin;
            TwoBody_States_ph.push_back(new TwoBody_Infinite(i,j,Nx,Ny,Nz,Sz,Tz));
        }
    for(int i=A;i<numberSP;i++)
        for(int j=i+1;j<numberSP;j++)
        {
            SP_Infinite* I=(SP_Infinite*)SP_States[i];
            SP_Infinite* J=(SP_Infinite*)SP_States[j];
            int Nx=I->nx + J->nx;
            int Ny=I->ny + J->ny;
            int Nz=I->nz + J->nz;
            int Sz=I->spin + J->spin;
            int Tz=I->isospin + J->isospin;
            TwoBody_States_pp.push_back(new TwoBody_Infinite(i,j,Nx,Ny,Nz,Sz,Tz));
        }
    sort(TwoBody_States_hh.begin(),TwoBody_States_hh.end(),TwoBody_compare());
    sort(TwoBody_States_ph.begin(),TwoBody_States_ph.end(),TwoBody_compare());
    sort(TwoBody_States_pp.begin(),TwoBody_States_pp.end(),TwoBody_compare());
}

void Infinite::printTwoBody_States()//TODO
{
    for(int i=0;i<TwoBody_States_pp.size()-1;i++)
    {
        TwoBody_Infinite* P1=(TwoBody_Infinite*)TwoBody_States_pp[i];
        cout<<i+1<<"\t";
        cout<<P1->p<<"\t";
        cout<<P1->q<<"\t";
        cout<<P1->Sz<<"\t";
        cout<<P1->Tz<<"\t";
        cout<<P1->Nx<<"\t";
        cout<<P1->Ny<<"\t";
        cout<<P1->Nz<<endl;
    }
}

void Infinite::generateConfigurations()
{

}

void Infinite::printConfigurations()
{

}

double Infinite::V1B(int a,int b)
{
    if (a == b)
        return ((SP_Infinite*)SP_States[a])->HF_spEnergy;
    return 0.0;
}

double Infinite::V2B(int p,int q,int r,int s)
{

    return V2B_sym(p,q,r,s)-V2B_sym(p,q,s,r);
}

double Infinite::V2B_sym(int p, int q, int r, int s)
{
    double v = 0.0;
    SP_Infinite* P = (SP_Infinite*)SP_States[p];
    SP_Infinite* Q = (SP_Infinite*)SP_States[q];
    SP_Infinite* R = (SP_Infinite*)SP_States[r];
    SP_Infinite* S = (SP_Infinite*)SP_States[s];
    if ((P->nx + Q->nx == R->nx + S->nx) &&
            (P->ny + Q->ny == R->ny + S->ny) &&
            (P->nz + Q->nz == R->nz + S->nz) &&
            (P->spin + Q->spin == R->spin + S->spin) &&
            (P->isospin + Q->isospin == R->isospin + S->isospin))
    {
        double q2 = pow(P->kx - R->kx, 2) + pow(P->ky - R->ky, 2) + pow(P->kz - R->kz, 2);
        double vr = V0R / pow(L,3) * pow(M_PI/KR, 3./2) * exp(-q2/4/KR);
        double vt = V0T / pow(L,3) * pow(M_PI/KT, 3./2) * exp(-q2/4/KT);
        double vs = V0S / pow(L,3) * pow(M_PI/KS, 3./2) * exp(-q2/4/KS);
        v += 0.5*vr*(deltaSpinIsospin(P,Q,R,S)-deltaSpinIsospin(P,Q,S,R));
        v += 0.25*vt*(deltaSpinIsospin(P,Q,R,S)-deltaSpin(P,Q,R,S)*deltaIsospin(P,Q,S,R)
                      +deltaSpin(P,Q,S,R)*deltaIsospin(P,Q,R,S)-deltaSpinIsospin(P,Q,S,R));
        v += 0.25*vs*(deltaSpinIsospin(P,Q,R,S)+deltaSpin(P,Q,R,S)*deltaIsospin(P,Q,S,R)
                      -deltaSpin(P,Q,S,R)*deltaIsospin(P,Q,R,S)-deltaSpinIsospin(P,Q,S,R));
    }
    return v;
}


bool Infinite::deltaSpin(SP_Infinite* P, SP_Infinite* Q, SP_Infinite* R, SP_Infinite* S)
{
    return ((P->spin == R->spin) && (Q->spin == S->spin));
}

bool Infinite::deltaIsospin(SP_Infinite* P, SP_Infinite* Q, SP_Infinite* R, SP_Infinite* S)
{
    return ((P->isospin == R->isospin) && (Q->isospin == S->isospin));
}

bool Infinite::deltaSpinIsospin(SP_Infinite* P, SP_Infinite* Q, SP_Infinite* R, SP_Infinite* S)
{
    return ((P->spin == R->spin) && (P->isospin == R->isospin) &&
            (Q->spin == S->spin) && (Q->isospin == S->isospin));
}

void Infinite::HF_calculateE0()
{
    HF_E0 = 0.0;
    for (int i = 0; i < numberSP; i++)
    {
        double HF_T = SP_States[i]->spEnergy;
        double HF_V = 0.0;
        for (int j = 0; j < A; j++)
            HF_V += V2B(i,j,i,j);
        ((SP_Infinite*)SP_States[i])->HF_spEnergy = HF_T + HF_V;
        if (i < A)
            HF_E0 += HF_T + 0.5 * HF_V;
    }
}

double Infinite::HF_exact_f(double r)
{
    int vs;
    int vt;
    if(g_s==4)
    {
        vs=2;
        vt=2;
    }
    else
    {
        vs=2;
        vt=1;
    }
    double VR=V0R*exp(-KR*r*r);
    double VS=V0S*exp(-KS*r*r);
    double VT=V0T*exp(-KT*r*r);
    double V1=0.5*VR+0.25*VT+0.25*VS;
    double V2=0.25*VT-0.25*VS;
    double V3=0.25*VS-0.25*VT;
    double V4=-0.5*VR-0.25*VT-0.25*VS;
    return vs*vs*vt*vt*V1 + vs*vt*vt*V2 + vt*vs*vs*V3 + vs*vt*V4 -
            pow(3*J1(k_F*r)/(k_F*r),2)*
            (vs*vs*vt*vt*V4 + vs*vt*vt*V3 + vt*vs*vs*V2 + vs*vt*V1);
}
void Infinite::HF_cal_exact_E0()
{
    HF_exact_E0=0.3 * hbar*hbar /m * k_F*k_F;
    for(int i=0;i<gauss_x.size();i++)
    {
        double temp=M_PI*0.25*(gauss_x[i] + 1);
        double r=tan(temp);
        HF_exact_E0+=rho * 2*M_PI/(g_s*g_s) * r*r * HF_exact_f(r) * gauss_w[i]* M_PI*0.25/pow(cos(temp),2);
    }
}

void Infinite::CCD_generateBlockMatrices()
{
    for (int i = 0, j = 0; i < TwoBody_States_hh.size(); i++)
    {
        TwoBody_Infinite *State_hh = (TwoBody_Infinite *)TwoBody_States_hh[i];
        TwoBody_Infinite *State_pp = (TwoBody_Infinite *)TwoBody_States_pp[j];
        TwoBody_compare LessThan;
        while( (LessThan(State_pp,State_hh)) && (j < TwoBody_States_pp.size()))
        {
            j++;
            State_pp = (TwoBody_Infinite*)TwoBody_States_pp[j];
        }
        if( j >= TwoBody_States_pp.size() )
            break;
        if( !LessThan(State_hh,State_pp) )
        {
            int Nx = State_hh->Nx;
            int Ny = State_hh->Ny;
            int Nz = State_hh->Nz;
            int Sz = State_hh->Sz;
            int Tz = State_hh->Tz;
            CCD_V_hhpp.push_back(Channel(Nx,Ny,Nz,Sz,Tz));
            CCD_V_hhpp[CCD_V_hhpp.size()-1].bra.push_back(i);
            CCD_V_hhpp[CCD_V_hhpp.size()-1].ket.push_back(j);
            if( i < TwoBody_States_hh.size()-1 )
            {
                TwoBody_Infinite *Ph = (TwoBody_Infinite*)TwoBody_States_hh[i+1];
                while( (!LessThan(Ph,State_hh)) && (!LessThan(State_hh,Ph)) )
                {
                    i++;
                    CCD_V_hhpp[CCD_V_hhpp.size()-1].bra.push_back(i);
                    if( i >= TwoBody_States_hh.size()-1) break;
                    State_hh = (TwoBody_Infinite*)TwoBody_States_hh[i];
                    Ph = (TwoBody_Infinite*)TwoBody_States_hh[i+1];
                }
            }
            if(j < TwoBody_States_pp.size()-1 )
            {
                TwoBody_Infinite *Pp = (TwoBody_Infinite*)TwoBody_States_pp[j+1];
                while( (!LessThan(Pp,State_pp)) && (!LessThan(State_pp,Pp)) )
                {
                    j++;
                    CCD_V_hhpp[CCD_V_hhpp.size()-1].ket.push_back(j);
                    if(j >= TwoBody_States_pp.size()-1) break;
                    State_pp = (TwoBody_Infinite*)TwoBody_States_pp[j];
                    Pp = (TwoBody_Infinite*)TwoBody_States_pp[j+1];
                }
            }
            j++;
        }
    }

    for(int i = 0; i < CCD_V_hhpp.size(); i++)
    {
        int Nx = CCD_V_hhpp[i].Nx;
        int Ny = CCD_V_hhpp[i].Ny;
        int Nz = CCD_V_hhpp[i].Nz;
        int Sz = CCD_V_hhpp[i].Sz;
        int Tz = CCD_V_hhpp[i].Tz;
        CCD_V_hhhh.push_back(Channel(Nx,Ny,Nz,Sz,Tz));
        int size = CCD_V_hhpp[i].bra.size();
        CCD_V_hhhh[CCD_V_hhhh.size()-1].bra.resize(size);
        CCD_V_hhhh[CCD_V_hhhh.size()-1].ket.resize(size);
        for(int j = 0; j < size; j++)
        {
            CCD_V_hhhh[CCD_V_hhhh.size()-1].bra[j] = CCD_V_hhpp[i].bra[j];
            CCD_V_hhhh[CCD_V_hhhh.size()-1].ket[j] = CCD_V_hhpp[i].bra[j];
        }
    }


    for(int i = 0; i < CCD_V_hhpp.size(); i++)
    {
        int Nx = CCD_V_hhpp[i].Nx;
        int Ny = CCD_V_hhpp[i].Ny;
        int Nz = CCD_V_hhpp[i].Nz;
        int Sz = CCD_V_hhpp[i].Sz;
        int Tz = CCD_V_hhpp[i].Tz;
        CCD_V_pppp.push_back(Channel(Nx,Ny,Nz,Sz,Tz));
        int size = CCD_V_hhpp[i].ket.size();
        CCD_V_pppp[CCD_V_pppp.size()-1].bra.resize(size);
        CCD_V_pppp[CCD_V_pppp.size()-1].ket.resize(size);
        for(int j = 0; j < size; j++)
        {
            CCD_V_pppp[CCD_V_pppp.size()-1].bra[j] = CCD_V_hhpp[i].ket[j];
            CCD_V_pppp[CCD_V_pppp.size()-1].ket[j] = CCD_V_hhpp[i].ket[j];
        }
    }

    CCD_T_hhpp = CCD_V_hhpp;
    CCD_e_hhpp = CCD_V_hhpp;
    for(int i = 0; i < CCD_V_hhpp.size(); i++)
    {
        int bra_dim = CCD_V_hhpp[i].bra.size();
        int ket_dim = CCD_V_hhpp[i].ket.size();
        CCD_V_hhpp[i].mat.resize(bra_dim,ket_dim);
        CCD_T_hhpp[i].mat.resize(bra_dim,ket_dim);
        CCD_e_hhpp[i].mat.resize(bra_dim,ket_dim);
        for(int Ibra=0; Ibra < bra_dim; Ibra++)
        {
            for(int Iket = 0; Iket < ket_dim; Iket++)
            {
                TwoBody_Infinite *State_bra = (TwoBody_Infinite*)TwoBody_States_hh[ CCD_V_hhpp[i].bra[Ibra] ];
                TwoBody_Infinite *State_ket = (TwoBody_Infinite*)TwoBody_States_pp[ CCD_V_hhpp[i].ket[Iket] ];
                int bra_p = State_bra->p;
                int bra_q = State_bra->q;
                int ket_p = State_ket->p;
                int ket_q = State_ket->q;
                double V = V2B(bra_p,bra_q,ket_p,ket_q);
                double de = 1.0/(V1B(bra_p,bra_p)+V1B(bra_q,bra_q)-V1B(ket_p,ket_p)-V1B(ket_q,ket_q));
                CCD_V_hhpp[i].mat(Ibra,Iket) = V;
                CCD_e_hhpp[i].mat(Ibra,Iket) = de;
                CCD_T_hhpp[i].mat(Ibra,Iket) = V*de;
            }
        }
    }

    for(int i = 0; i < CCD_V_hhhh.size(); i++)
    {
        int dim = CCD_V_hhhh[i].bra.size();
        CCD_V_hhhh[i].mat.resize(dim,dim);
        for(int Ibra=0; Ibra < dim; Ibra++)
        {
            for(int Iket = Ibra; Iket < dim; Iket++)
            {
                TwoBody_Infinite *State_bra = (TwoBody_Infinite*)TwoBody_States_hh[ CCD_V_hhhh[i].bra[Ibra] ];
                TwoBody_Infinite *State_ket = (TwoBody_Infinite*)TwoBody_States_hh[ CCD_V_hhhh[i].ket[Iket] ] ;
                int bra_p = State_bra->p;
                int bra_q = State_bra->q;
                int ket_p = State_ket->p;
                int ket_q = State_ket->q;
                CCD_V_hhhh[i].mat(Ibra,Iket) = V2B(bra_p,bra_q,ket_p,ket_q);
                CCD_V_hhhh[i].mat(Iket,Ibra) = CCD_V_hhhh[i].mat(Ibra,Iket);
            }
        }
    }

    for(int i = 0; i < CCD_V_pppp.size(); i++)
    {
        int dim = CCD_V_pppp[i].bra.size();
        CCD_V_pppp[i].mat.resize(dim,dim);
        for(int Ibra = 0; Ibra < dim; Ibra++)
        {
            for(int Iket = Ibra; Iket < dim; Iket++)
            {
                TwoBody_Infinite *State_bra = (TwoBody_Infinite*)TwoBody_States_pp[ CCD_V_pppp[i].bra[Ibra] ];
                TwoBody_Infinite *State_ket = (TwoBody_Infinite*)TwoBody_States_pp[ CCD_V_pppp[i].ket[Iket] ];
                int bra_p = State_bra->p;
                int bra_q = State_bra->q;
                int ket_p = State_ket->p;
                int ket_q = State_ket->q;
                CCD_V_pppp[i].mat(Ibra,Iket) = V2B(bra_p,bra_q,ket_p,ket_q);
                CCD_V_pppp[i].mat(Iket,Ibra) = CCD_V_pppp[i].mat(Ibra,Iket);
            }
        }
    }
}

void Infinite::CCD_BlockMatrices()
{
    generateTwoBody_States();
    CCD_generateBlockMatrices();
    vector<Channel> CCD_T_old;
    int iter = 0;
    CCD_deltaE = 0.0;
    double CCD_deltaE_old;
    do
    {
        iter++;
        CCD_deltaE_old = CCD_deltaE;
        CCD_T_old = CCD_T_hhpp;
        for(int i = 0; i < CCD_V_hhpp.size(); i++)
        {
            CCD_T_hhpp[i].mat = CCD_V_hhpp[i].mat
                    + CCD_T_old[i].mat * CCD_V_pppp[i].mat
                          + CCD_V_hhhh[i].mat * CCD_T_old[i].mat;
            CCD_T_hhpp[i].mat = CCD_T_hhpp[i].mat.cwiseProduct(CCD_e_hhpp[i].mat);
        }
        double sum = 0;
        for(int i = 0; i < CCD_V_hhpp.size(); i++)
        {
            sum += (CCD_V_hhpp[i].mat*(CCD_T_hhpp[i].mat.transpose())).trace();
        }
        CCD_deltaE = sum;
        cout << iter << "\t" << setprecision(10) << sum << endl;
    }
    while (abs(CCD_deltaE_old - CCD_deltaE) > 1e-8);
}

void Infinite::CCD_BlockMatricesIntermediates()
{
    generateTwoBody_States();
    CCD_generateBlockMatrices();
    vector<Channel> CCD_T_old;
    int iter = 0;
    CCD_deltaE = 0.0;

    vector<Channel> CCD_intermediate_pp, CCD_intermediate_hh,
            CCD_intermediate_hhhh, CCD_intermediate_hpph, CCD_intermediate_pppp;
    CCD_intermediate_pppp = CCD_V_hhpp;
    CCD_intermediate_hhhh = CCD_V_hhpp;

    double CCD_deltaE_old;
    do
    {
        iter++;
        CCD_deltaE_old = CCD_deltaE;
        CCD_T_old = CCD_T_hhpp;

        for(int i = 0; i < CCD_V_hhpp.size(); i++)
        {
            CCD_intermediate_pppp[i].mat = CCD_V_pppp[i].mat;
            CCD_intermediate_hhhh[i].mat = CCD_V_hhhh[i].mat
                    + CCD_V_hhpp[i].mat * CCD_T_hhpp[i].mat.transpose();
        }

        for(int i = 0; i < CCD_V_hhpp.size(); i++)
        {
            CCD_T_hhpp[i].mat = CCD_V_hhpp[i].mat
                    + CCD_T_old[i].mat * CCD_intermediate_pppp[i].mat
                          + CCD_intermediate_hhhh[i].mat * CCD_T_old[i].mat;
            CCD_T_hhpp[i].mat = CCD_T_hhpp[i].mat.cwiseProduct(CCD_e_hhpp[i].mat);
        }
        double sum = 0;
        for(int i = 0; i < CCD_V_hhpp.size(); i++)
        {
            sum += (CCD_V_hhpp[i].mat*(CCD_T_hhpp[i].mat.transpose())).trace();
        }
        CCD_deltaE = sum;
        cout << iter << "\t" << setprecision(10) << sum << endl;
    }
    while (abs(CCD_deltaE_old - CCD_deltaE) > 1e-8);
}

//void Infinite::MatrixElement(vector<Channel> &channels, int p, int q, int r, int s)
//{
//    SP_Infinite* P = (SP_Infinite*)SP_States[p];
//    SP_Infinite* Q = (SP_Infinite*)SP_States[q];
//    int Nx_bra = P->nx + Q->nx;
//    int Ny_bra = P->ny + Q->ny;
//    int Nz_bra = P->nz + Q->nz;
//    int Sz_bra = P->spin + Q->spin;
//    int Tz_bra = P->isospin + Q->isospin;
//    TwoBody_Infinite bra = TwoBody_Infinite(p,q,Nx_bra,Ny_bra,Nz_bra,Sz_bra,Tz_bra);

//    SP_Infinite* R = (SP_Infinite*)SP_States[r];
//    SP_Infinite* S = (SP_Infinite*)SP_States[s];
//    int Nx_ket = R->nx + S->nx;
//    int Ny_ket = R->ny + S->ny;
//    int Nz_ket = R->nz + S->nz;
//    int Sz_ket = R->spin + S->spin;
//    int Tz_ket = R->isospin + S->isospin;
//    TwoBody_Infinite ket = TwoBody_Infinite(r,s,Nx_ket,Ny_ket,Nz_ket,Sz_ket,Tz_ket);

//    int i = -1, j = -1;
//    while ( !(channels[i++ +1].Include(&bra)) );
//    while ( !(channels[j++ +1].Include(&ket)) );

//}
