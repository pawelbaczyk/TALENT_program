#include "infinite.h"

bool SP_Infinite_compare::operator() (SP_State *P1, SP_State *P2)
{
    if (P1->spEnergy != P2->spEnergy)
        return P1->spEnergy < P2->spEnergy;
    else
    {
        SP_Infinite* S1 = (SP_Infinite*)P1;
        SP_Infinite* S2 = (SP_Infinite*)P2;
        if (S1->kx != S2->kx)
            return S1->kx < S2->kx;
        if (S1->ky != S2->ky)
            return S1->ky < S2->ky;
        if (S1->kz != S2->kz)
            return S1->kz < S2->kz;
        if (S1->spin != S2->spin)
            return S1->spin > S2->spin;
        if (S1->isospin != S2->isospin)
            return S1->isospin > S2->isospin;
        else return 0;
    }

}

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

bool Channel::Include(TwoBody_Infinite &TwoBody)
{
    if((TwoBody.Nx == Nx) && (TwoBody.Ny == Ny) && (TwoBody.Nz == Nz) &&
            (TwoBody.Sz == Sz) && (TwoBody.Tz==Tz))
        return true;
    return false;
}

Infinite::Infinite(int _A, int _g_s, double _rho, int _nMax)
    : System(_A), g_s(_g_s), nMax(_nMax), gauss_x(100), gauss_w(100)
{
    thetaX = 0.0;
    thetaY = 0.0;
    thetaZ = 0.0;
    setRho(_rho);
    gauleg(-1,1,gauss_x,gauss_w);
}

Infinite::~Infinite()
{
    for (int i = 0; i < TwoBody_States_hmp.size(); i++)
    {
        delete TwoBody_States_hmp[i];
    }
    for (int i = 0; i < TwoBody_States_hpm.size(); i++)
    {
        delete TwoBody_States_hpm[i];
    }
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
                        double kx = (2 * M_PI * nX + thetaX) / L;
                        double ky = (2 * M_PI * nY + thetaY) / L;
                        double kz = (2 * M_PI * nZ + thetaZ) / L;
                        if (g_s == 2 || g_s == 4)
                        {
                            SP_States.push_back(new SP_Infinite(nX,nY,nZ,kx,ky,kz,1,1));
                            SP_States.push_back(new SP_Infinite(nX,nY,nZ,kx,ky,kz,-1,1));
                        }
                        if (g_s == 4)
                        {
                            SP_States.push_back(new SP_Infinite(nX,nY,nZ,kx,ky,kz,1,-1));
                            SP_States.push_back(new SP_Infinite(nX,nY,nZ,kx,ky,kz,-1,-1));
                        }
                    }
    sort(SP_States.begin(),SP_States.end(),SP_Infinite_compare());
    numberSP = SP_States.size();

    HF_calculateE0();//to set HF s.p. energies
}

void Infinite::generateMagicNumbers()
{
    magicNumbers.clear();
    double energy = SP_States[0]->spEnergy;
    for (int i = 0; i < numberSP; i++)
    {
        if (abs(SP_States[i]->spEnergy - energy) > 1e-6)
        {
            energy = SP_States[i]->spEnergy;
            magicNumbers.push_back(i);
        }
    }
    cout << "MAGIC NUMBERS" << endl;
    cout << "Shell\tCapacity\n";
    cout << magicNumbers[0] << "\t" << magicNumbers[0] << "\n";
    for (int i = 1; i < magicNumbers.size(); i++)
        cout << magicNumbers[i] << "\t" << magicNumbers[i]-magicNumbers[i-1] << "\n";
}

void Infinite::printSP_States()
{
    cout << "***************************************************************\n";
    cout << "*                   Single particles states                   *\n";
    cout << "***************************************************************\n";
    cout << "Number\t" << "Energy\t\t" << "nx\tny\tnz\t" << "Spin\t" << "Isospin\n";
    cout << "***************************************************************\n";
    for (int i = 0; i < SP_States.size(); i++)
    {
        SP_Infinite* state = (SP_Infinite*)SP_States[i];
        cout << std::fixed << i+1 << "\t" << state->spEnergy << "\t"
             << state->nx << "\t" << state->ny << "\t" << state->nz << "\t"
             << (state->spin == 1 ? "+" : "-") << "\t" << (state->isospin == 1 ? "n" : "p")
             << endl;
        if (i+1 == A)
            cout << "************************* Fermi level *************************\n";
    }
    cout << "***************************************************************\n";
}

void Infinite::generateTwoBody_States()
{
    TwoBody_States_hh.clear();
    for (int i = 0; i < A; i++)
        for (int j = i+1; j < A; j++)
        {
            SP_Infinite* I = (SP_Infinite*)SP_States[i];
            SP_Infinite* J = (SP_Infinite*)SP_States[j];
            int Nx = I->nx + J->nx;
            int Ny = I->ny + J->ny;
            int Nz = I->nz + J->nz;
            int Sz = (I->spin + J->spin)/2;
            int Tz = (I->isospin + J->isospin)/2;
            TwoBody_States_hh.push_back(new TwoBody_Infinite(i,j,Nx,Ny,Nz,Sz,Tz));
        }
    sort(TwoBody_States_hh.begin(),TwoBody_States_hh.end(),TwoBody_compare());

    TwoBody_States_hp.clear();
    for (int i = 0; i < A; i++)
        for (int j = A; j < numberSP; j++)
        {
            SP_Infinite* I = (SP_Infinite*)SP_States[i];
            SP_Infinite* J = (SP_Infinite*)SP_States[j];
            int Nx = I->nx + J->nx;
            int Ny = I->ny + J->ny;
            int Nz = I->nz + J->nz;
            int Sz = (I->spin + J->spin)/2;
            int Tz = (I->isospin + J->isospin)/2;
            TwoBody_States_hp.push_back(new TwoBody_Infinite(i,j,Nx,Ny,Nz,Sz,Tz));
        }
    sort(TwoBody_States_hp.begin(),TwoBody_States_hp.end(),TwoBody_compare());

    TwoBody_States_pp.clear();
    for (int i = A; i < numberSP; i++)
        for (int j = i+1; j < numberSP; j++)
        {
            SP_Infinite* I = (SP_Infinite*)SP_States[i];
            SP_Infinite* J = (SP_Infinite*)SP_States[j];
            int Nx = I->nx + J->nx;
            int Ny = I->ny + J->ny;
            int Nz = I->nz + J->nz;
            int Sz = (I->spin + J->spin)/2;
            int Tz = (I->isospin + J->isospin)/2;
            TwoBody_States_pp.push_back(new TwoBody_Infinite(i,j,Nx,Ny,Nz,Sz,Tz));
        }
    sort(TwoBody_States_pp.begin(),TwoBody_States_pp.end(),TwoBody_compare());

    TwoBody_States_hmp.clear();
    for (int i = 0; i < A; i++)
        for (int j = A; j <numberSP; j++)
        {
            SP_Infinite* I = (SP_Infinite*)SP_States[i];
            SP_Infinite* J = (SP_Infinite*)SP_States[j];
            int Nx = -I->nx + J->nx;
            int Ny = -I->ny + J->ny;
            int Nz = -I->nz + J->nz;
            int Sz = (-I->spin + J->spin)/2;
            int Tz = (-I->isospin + J->isospin)/2;
            TwoBody_States_hmp.push_back(new TwoBody_Infinite(i,j,Nx,Ny,Nz,Sz,Tz));
        }
    sort(TwoBody_States_hmp.begin(),TwoBody_States_hmp.end(),TwoBody_compare());

    TwoBody_States_hpm.clear();
    for (int i = 0; i < A; i++)
        for (int j = A; j <numberSP; j++)
        {
            SP_Infinite* I = (SP_Infinite*)SP_States[i];
            SP_Infinite* J = (SP_Infinite*)SP_States[j];
            int Nx = I->nx - J->nx;
            int Ny = I->ny - J->ny;
            int Nz = I->nz - J->nz;
            int Sz = (I->spin - J->spin)/2;
            int Tz = (I->isospin - J->isospin)/2;
            TwoBody_States_hpm.push_back(new TwoBody_Infinite(i,j,Nx,Ny,Nz,Sz,Tz));
        }
    sort(TwoBody_States_hpm.begin(),TwoBody_States_hpm.end(),TwoBody_compare());
}

void Infinite::printTwoBody_States()
{
    cout << "***************************************************************\n";
    cout << "*                    Two particle states hh                   *\n";
    cout << "***************************************************************\n";
    cout << "Number\t" << "p\tq\t" << "nx\tny\tnz\t" << "Spin\t" << "Isospin\n";
    cout << "***************************************************************\n";
    for (int i = 0; i < TwoBody_States_hh.size(); i++)
    {
        TwoBody_Infinite* state = (TwoBody_Infinite*)TwoBody_States_hh[i];
        cout << i+1 << "\t" << state->p+1 << "\t" << state->q+1 << "\t"
             << state->Nx<<"\t" << state->Ny << "\t" << state->Nz << "\t"
             << state->Sz << "\t" << state->Tz << "\n";
    }
    cout << "***************************************************************\n";

    cout << "***************************************************************\n";
    cout << "*                    Two particle states hp                   *\n";
    cout << "***************************************************************\n";
    cout << "Number\t" << "p\tq\t" << "nx\tny\tnz\t" << "Spin\t" << "Isospin\n";
    cout << "***************************************************************\n";
    for (int i = 0; i < TwoBody_States_hp.size(); i++)
    {
        TwoBody_Infinite* state = (TwoBody_Infinite*)TwoBody_States_hp[i];
        cout << i+1 << "\t" << state->p+1 << "\t" << state->q+1 << "\t"
             << state->Nx<<"\t" << state->Ny << "\t" << state->Nz << "\t"
             << state->Sz << "\t" << state->Tz << "\n";
    }
    cout << "***************************************************************\n";

    cout << "***************************************************************\n";
    cout << "*                    Two particle states pp                   *\n";
    cout << "***************************************************************\n";
    cout << "Number\t" << "p\tq\t" << "nx\tny\tnz\t" << "Spin\t" << "Isospin\n";
    cout << "***************************************************************\n";
    for (int i = 0; i < TwoBody_States_pp.size(); i++)
    {
        TwoBody_Infinite* state = (TwoBody_Infinite*)TwoBody_States_pp[i];
        cout << i+1 << "\t" << state->p+1 << "\t" << state->q+1 << "\t"
             << state->Nx<<"\t" << state->Ny << "\t" << state->Nz << "\t"
             << state->Sz << "\t" << state->Tz << "\n";
    }
    cout << "***************************************************************\n";

    cout << "***************************************************************\n";
    cout << "*                   Two particle states h-1p                  *\n";
    cout << "***************************************************************\n";
    cout << "Number\t" << "p\tq\t" << "nx\tny\tnz\t" << "Spin\t" << "Isospin\n";
    cout << "***************************************************************\n";
    for (int i = 0; i < TwoBody_States_hmp.size(); i++)
    {
        TwoBody_Infinite* state = (TwoBody_Infinite*)TwoBody_States_hmp[i];
        cout << i+1 << "\t" << state->p+1 << "\t" << state->q+1 << "\t"
             << state->Nx<<"\t" << state->Ny << "\t" << state->Nz << "\t"
             << state->Sz << "\t" << state->Tz << "\n";
    }
    cout << "***************************************************************\n";

    cout << "***************************************************************\n";
    cout << "*                   Two particle states hp-1                  *\n";
    cout << "***************************************************************\n";
    cout << "Number\t" << "p\tq\t" << "nx\tny\tnz\t" << "Spin\t" << "Isospin\n";
    cout << "***************************************************************\n";
    for (int i = 0; i < TwoBody_States_hpm.size(); i++)
    {
        TwoBody_Infinite* state = (TwoBody_Infinite*)TwoBody_States_hpm[i];
        cout << i+1 << "\t" << state->p+1 << "\t" << state->q+1 << "\t"
             << state->Nx<<"\t" << state->Ny << "\t" << state->Nz << "\t"
             << state->Sz << "\t" << state->Tz << "\n";
    }
    cout << "***************************************************************\n";
}

void Infinite::generateConfigurations()
{

}

void Infinite::printConfigurations()
{

}

void Infinite::setRho(double _rho)
{
    rho = _rho;
    generateSP_States(A);
}

void Infinite::setA(double _A)
{
    A = _A;
    generateSP_States(A);
}

void Infinite::setTheta(double _thetaX, double _thetaY, double _thetaZ)
{
    thetaX = _thetaX;
    thetaY = _thetaY;
    thetaZ = _thetaZ;
    generateSP_States(A);
}

void Infinite::TA_calculateE0(int nA, double (Infinite::*func)(), double (Infinite::*ref)())
{
    double refValue = (this->*ref)();
    double SP_diff = abs(refValue);
    SP_deltaE = 0.0;
    TA_deltaE = 0.0;
    twisted_x.resize(nA);
    twisted_w.resize(nA);
    gauleg(0,M_PI,twisted_x,twisted_w);
    for (int x = 0; x < twisted_x.size(); x++)
        for (int y = 0; y < twisted_x.size(); y++)
            for (int z = 0; z < twisted_x.size(); z++)
            {
                setTheta(twisted_x[x],twisted_x[y],twisted_x[z]);
                double value = (this->*func)();
                TA_deltaE += twisted_w[x] * twisted_w[y] * twisted_w[z] * value;
                if (abs(value/A - refValue) < SP_diff)
                {
                    SP_diff = abs(value/A - refValue);
                    SP_x = x;
                    SP_y = y;
                    SP_z = z;
                    SP_deltaE = value;
                }
            }
    TA_deltaE /= pow(M_PI,3);
    setTheta(0.0,0.0,0.0);
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

double Infinite::HF_calculateE0()
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
    return HF_E0;
}

double Infinite::HF_calculateT0()
{
    double HF_T0 = 0.0;
    for (int j = 0; j < A; j++)
    {
        double HF_T = SP_States[j]->spEnergy;
        HF_T0 += HF_T;
    }
    return HF_T0;
}

double Infinite::HF_calculateV0()
{
    double HF_V0 = 0.0;
    for (int i = 0; i < A; i++)
        for (int j = 0; j < A; j++)
            HF_V0 += 0.5 * V2B(i,j,i,j);
    return HF_V0;
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

double Infinite::HF_exactE0()
{
    HF_exact_E0 = 0.3 * hbar*hbar / m * k_F * k_F;
    for(int i = 0; i < gauss_x.size(); i++)
    {
        double temp = M_PI * 0.25 * (gauss_x[i] + 1);
        double r = tan(temp);
        HF_exact_E0 += rho * 2*M_PI/(g_s*g_s)
                * r * r * HF_exact_f(r) * gauss_w[i] * M_PI * 0.25 / pow(cos(temp),2);
    }
    return HF_exact_E0;
}

double Infinite::HF_exactT0()
{
    return 0.3 * hbar*hbar / m * k_F * k_F;
}

double Infinite::HF_exactV0()
{
    double HF_exact_V0 = 0.0;
    for(int i = 0; i < gauss_x.size(); i++)
    {
        double temp = M_PI * 0.25 * (gauss_x[i] + 1);
        double r = tan(temp);
        HF_exact_V0 += rho * 2*M_PI/(g_s*g_s)
                * r * r * HF_exact_f(r) * gauss_w[i] * M_PI * 0.25 / pow(cos(temp),2);
    }
    return HF_exact_V0;
}

void Infinite::CCD_generateBlockMatrices()
{
    CCD_V_hhhh.clear();
    CCD_V_hhpp.clear();
    CCD_V_pppp.clear();
    CCD_V_hpmhpm.clear();
    CCD_V_hmphpm.clear();
    CCD_T_hhpp.clear();
    CCD_position = Matrix<Position,Dynamic,Dynamic>();
    
    //hhpp
    CCD_position.resize(numberSP,numberSP);
    TwoBody_compare LessThan;
    for (int i = 0, j = 0; i < TwoBody_States_hh.size(); i++)
    {
        TwoBody_Infinite *State_hh = (TwoBody_Infinite *)TwoBody_States_hh[i];
        TwoBody_Infinite *State_pp = (TwoBody_Infinite *)TwoBody_States_pp[j];
        while( (LessThan(State_pp,State_hh)) && (j < TwoBody_States_pp.size()-1))
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
            int channelNumber = CCD_V_hhpp.size()-1;

            CCD_V_hhpp[channelNumber].bra.push_back(i);
            CCD_position(State_hh->p,State_hh->q) =
                    Position(channelNumber, CCD_V_hhpp[channelNumber].bra.size()-1);
            CCD_position(State_hh->q,State_hh->p) =
                    Position(channelNumber, CCD_V_hhpp[channelNumber].bra.size()-1);

            CCD_V_hhpp[channelNumber].ket.push_back(j);
            CCD_position(State_pp->p,State_pp->q) =
                    Position(channelNumber, CCD_V_hhpp[channelNumber].ket.size()-1);
            CCD_position(State_pp->q,State_pp->p) =
                    Position(channelNumber, CCD_V_hhpp[channelNumber].ket.size()-1);

            if( i < TwoBody_States_hh.size()-1 )
            {
                TwoBody_Infinite *Ph = (TwoBody_Infinite*)TwoBody_States_hh[i+1];
                while( (!LessThan(Ph,State_hh)) && (!LessThan(State_hh,Ph)) )
                {
                    i++;
                    CCD_V_hhpp[channelNumber].bra.push_back(i);
                    CCD_position(Ph->p,Ph->q) =
                            Position(channelNumber, CCD_V_hhpp[channelNumber].bra.size()-1);
                    CCD_position(Ph->q,Ph->p) =
                            Position(channelNumber, CCD_V_hhpp[channelNumber].bra.size()-1);
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
                    CCD_V_hhpp[channelNumber].ket.push_back(j);
                    CCD_position(Pp->p,Pp->q) =
                            Position(channelNumber, CCD_V_hhpp[channelNumber].ket.size()-1);
                    CCD_position(Pp->q,Pp->p) =
                            Position(channelNumber, CCD_V_hhpp[channelNumber].ket.size()-1);
                    if(j >= TwoBody_States_pp.size()-1) break;
                    State_pp = (TwoBody_Infinite*)TwoBody_States_pp[j];
                    Pp = (TwoBody_Infinite*)TwoBody_States_pp[j+1];
                }
            }
            j++;
        }
    }

    //hhhh
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

    //pppp
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

    //interaction for hhpp
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

    //interaction for hhhh
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

    //interaction for pppp
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

    //hmphpm
    for (int i = 0, j = 0; (i < TwoBody_States_hmp.size()) && (j < TwoBody_States_hpm.size()); i++)
    {
        TwoBody_Infinite *State_hmp = (TwoBody_Infinite *)TwoBody_States_hmp[i];
        TwoBody_Infinite *State_hpm = (TwoBody_Infinite *)TwoBody_States_hpm[j];
        while( (LessThan(State_hpm,State_hmp)) && (j < TwoBody_States_hpm.size()-1))
        {
            j++;
            State_hpm = (TwoBody_Infinite*)TwoBody_States_hpm[j];
        }
        if( j >= TwoBody_States_hpm.size() )
            break;
        if( !LessThan(State_hmp,State_hpm) )
        {
            int Nx = State_hmp->Nx;
            int Ny = State_hmp->Ny;
            int Nz = State_hmp->Nz;
            int Sz = State_hmp->Sz;
            int Tz = State_hmp->Tz;

            CCD_V_hmphpm.push_back(Channel(Nx,Ny,Nz,Sz,Tz));
            int channelNumber = CCD_V_hmphpm.size()-1;

            CCD_V_hmphpm[channelNumber].bra.push_back(i);
            //CCD_position(State_hmp->p,State_hmp->q) =
            //        Position(channelNumber, CCD_V_hmphpm[channelNumber].bra.size()-1);
            //CCD_position(State_hmp->q,State_hmp->p) =
            //        Position(channelNumber, CCD_V_hmphpm[channelNumber].bra.size()-1);

            CCD_V_hmphpm[channelNumber].ket.push_back(j);
            CCD_position(State_hpm->p,State_hpm->q) =
                    Position(channelNumber, CCD_V_hmphpm[channelNumber].ket.size()-1);
            CCD_position(State_hpm->q,State_hpm->p) =
                    Position(channelNumber, CCD_V_hmphpm[channelNumber].ket.size()-1);

            if( i < TwoBody_States_hmp.size()-1 )
            {
                TwoBody_Infinite *Ph = (TwoBody_Infinite*)TwoBody_States_hmp[i+1];
                while( (!LessThan(Ph,State_hmp)) && (!LessThan(State_hmp,Ph)) )
                {
                    i++;
                    CCD_V_hmphpm[channelNumber].bra.push_back(i);
                    //CCD_position(Ph->p,Ph->q) =
                    //        Position(channelNumber, CCD_V_hmphpm[channelNumber].bra.size()-1);
                    //CCD_position(Ph->q,Ph->p) =
                    //        Position(channelNumber, CCD_V_hmphpm[channelNumber].bra.size()-1);
                    if( i >= TwoBody_States_hmp.size()-1) break;
                    State_hmp = (TwoBody_Infinite*)TwoBody_States_hmp[i];
                    Ph = (TwoBody_Infinite*)TwoBody_States_hmp[i+1];
                }
            }
            if(j < TwoBody_States_hpm.size()-1 )
            {
                TwoBody_Infinite *Pp = (TwoBody_Infinite*)TwoBody_States_hpm[j+1];
                while( (!LessThan(Pp,State_hpm)) && (!LessThan(State_hpm,Pp)) )
                {
                    j++;
                    CCD_V_hmphpm[channelNumber].ket.push_back(j);
                    CCD_position(Pp->p,Pp->q) =
                            Position(channelNumber, CCD_V_hmphpm[channelNumber].ket.size()-1);
                    CCD_position(Pp->q,Pp->p) =
                            Position(channelNumber, CCD_V_hmphpm[channelNumber].ket.size()-1);
                    if(j >= TwoBody_States_hpm.size()-1) break;
                    State_hpm = (TwoBody_Infinite*)TwoBody_States_hpm[j];
                    Pp = (TwoBody_Infinite*)TwoBody_States_hpm[j+1];
                }
            }
            j++;
        }
    }

    //hpmhpm
    for(int i = 0; i < CCD_V_hmphpm.size(); i++)
    {
        int Nx = CCD_V_hmphpm[i].Nx;
        int Ny = CCD_V_hmphpm[i].Ny;
        int Nz = CCD_V_hmphpm[i].Nz;
        int Sz = CCD_V_hmphpm[i].Sz;
        int Tz = CCD_V_hmphpm[i].Tz;
        CCD_V_hpmhpm.push_back(Channel(Nx,Ny,Nz,Sz,Tz));
        int size = CCD_V_hmphpm[i].ket.size();
        CCD_V_hpmhpm[CCD_V_hpmhpm.size()-1].bra.resize(size);
        CCD_V_hpmhpm[CCD_V_hpmhpm.size()-1].ket.resize(size);
        for(int j = 0; j < size; j++)
        {
            CCD_V_hpmhpm[CCD_V_hpmhpm.size()-1].bra[j] = CCD_V_hmphpm[i].ket[j];
            CCD_V_hpmhpm[CCD_V_hpmhpm.size()-1].ket[j] = CCD_V_hmphpm[i].ket[j];
        }
    }

    //interaction for hmphpm
    for(int i = 0; i < CCD_V_hmphpm.size(); i++)
    {
        int bra_dim = CCD_V_hmphpm[i].bra.size();
        int ket_dim = CCD_V_hmphpm[i].ket.size();
        CCD_V_hmphpm[i].mat.resize(bra_dim,ket_dim);
        for(int Ibra=0; Ibra < bra_dim; Ibra++)
        {
            for(int Iket = 0; Iket < ket_dim; Iket++)
            {
                TwoBody_Infinite *State_bra = (TwoBody_Infinite*)TwoBody_States_hmp[ CCD_V_hmphpm[i].bra[Ibra] ];
                TwoBody_Infinite *State_ket = (TwoBody_Infinite*)TwoBody_States_hpm[ CCD_V_hmphpm[i].ket[Iket] ];
                int bra_p = State_bra->p;
                int bra_q = State_bra->q;
                int ket_p = State_ket->p;
                int ket_q = State_ket->q;
                CCD_V_hmphpm[i].mat(Ibra,Iket) = V2B(bra_p,ket_p,bra_q,ket_q);
            }
        }
    }

    //interaction for hpmhpm
    for(int i = 0; i < CCD_V_hpmhpm.size(); i++)
    {
        int dim = CCD_V_hpmhpm[i].bra.size();
        CCD_V_hpmhpm[i].mat.resize(dim,dim);
        for(int Ibra = 0; Ibra < dim; Ibra++)
        {
            for(int Iket = Ibra; Iket < dim; Iket++)
            {
                TwoBody_Infinite *State_bra = (TwoBody_Infinite*)TwoBody_States_hpm[ CCD_V_hpmhpm[i].bra[Ibra] ];
                TwoBody_Infinite *State_ket = (TwoBody_Infinite*)TwoBody_States_hpm[ CCD_V_hpmhpm[i].ket[Iket] ];
                int bra_p = State_bra->p;
                int bra_q = State_bra->q;
                int ket_p = State_ket->p;
                int ket_q = State_ket->q;
                CCD_V_hpmhpm[i].mat(Ibra,Iket) = V2B(bra_p,ket_q,bra_q,ket_p);
                CCD_V_hpmhpm[i].mat(Iket,Ibra) = CCD_V_hpmhpm[i].mat(Ibra,Iket);
            }
        }
    }
    

    /////
//    for (int channel = 0; channel < CCD_V_hmphpm.size(); channel++)
//    {
//        for (int Ibra = 0; Ibra < CCD_V_hmphpm[channel].bra.size(); Ibra++)
//            for (int Iket = 0; Iket < CCD_V_hmphpm[channel].ket.size(); Iket++)
//            {
//                TwoBody_Infinite *State_bra = (TwoBody_Infinite*)TwoBody_States_hmp[ CCD_V_hmphpm[channel].bra[Ibra] ];
//                TwoBody_Infinite *State_ket = (TwoBody_Infinite*)TwoBody_States_hpm[ CCD_V_hmphpm[channel].ket[Iket] ];
//                int sign = +1;
//                int i = State_bra->p;
//                int b = State_bra->q;
//                int j = State_ket->p;
//                int a = State_ket->q;
//                if (b > a)
//                    sign *= -1;
//                if (i > j)
//                    sign *= -1;
//                Position pos_pp = CCD_position(a,b);
//                Position pos_hh = CCD_position(i,j);
//                if (pos_pp.channel != -1)
//                   if (abs(CCD_V_hmphpm[channel].mat(Ibra,Iket)
//                           - sign * CCD_V_hhpp[pos_pp.channel].mat(pos_hh.index,pos_pp.index)) > 1e-6)
//                       cout << "NIE\t" << CCD_V_hmphpm[channel].mat(Ibra,Iket)
//                            << "\t"<< sign * CCD_V_hhpp[pos_pp.channel].mat(pos_hh.index,pos_pp.index)
//                               << endl;
//            }
//    }



//        for (int i = 0; i < CCD_V_hpmhpm.size(); i++)
//        {
//            cout << i << endl << "\t";
//            for (int j = 0; j < CCD_V_hpmhpm[i].bra.size(); j++)
//                cout << CCD_V_hpmhpm[i].bra[j] << "\t";
//            cout << endl << "\t";
//            for (int j = 0; j < CCD_V_hpmhpm[i].ket.size(); j++)
//                cout << CCD_V_hpmhpm[i].ket[j] << "\t";
//            cout << endl << endl;
//        }

//        printTwoBody_States();
}

void Infinite::CCD_BlockMatricesLadders()
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
        //cout << iter << "\t" << setprecision(10) << sum << endl;
    }
    while (abs(CCD_deltaE_old - CCD_deltaE) > 1e-8);
}


void Infinite::CCD_BlockMatricesIntermediates()
{
    generateTwoBody_States();
    CCD_generateBlockMatrices();
    vector<Channel> CCD_T_old = CCD_T_hhpp;
    int iter = 0;
    CCD_deltaE = 0.0;

    vector<double> CCD_intermediate_1B(numberSP,0.0);
    vector<Channel> helpMatrix = CCD_V_pppp;

    vector<Channel> CCD_intermediate_hhhh, CCD_intermediate_pppp;
    CCD_intermediate_pppp = CCD_V_pppp;
    CCD_intermediate_hhhh = CCD_V_hhhh;

    vector<Channel> CCD_intermediate_hmphpm, CCD_intermediate_hpmhpm, CCD_T_hmphpm;
    CCD_intermediate_hmphpm = CCD_V_hmphpm;
    CCD_intermediate_hpmhpm = CCD_V_hpmhpm;
    CCD_T_hmphpm = CCD_V_hmphpm;

    double CCD_deltaE_old;
    do
    {
        iter++;
        CCD_deltaE_old = CCD_deltaE;
        //CCD_T_old = CCD_T_hhpp;
        double alpha = 0.0;
        for(int i = 0; i < CCD_V_hhpp.size(); i++)
            CCD_T_old[i].mat = alpha * CCD_T_old[i].mat + (1-alpha) * CCD_T_hhpp[i].mat;

        //one body intermediate
        //particle part
        for(int i = 0; i < CCD_V_hhpp.size(); i++)
            helpMatrix[i].mat = CCD_T_old[i].mat.transpose() * CCD_V_hhpp[i].mat;
        for (int b = A; b < numberSP; b++)
        {
            CCD_intermediate_1B[b] = V1B(b,b);
            for (int d = A; d < numberSP; d++)
            {
                Position pos = CCD_position(b,d);
                if (pos.channel != -1)
                    CCD_intermediate_1B[b] +=
                            -1.0 * helpMatrix[pos.channel].mat(pos.index,pos.index);
            }
        }
        //hole part
        helpMatrix = CCD_V_hhhh;
        for(int i = 0; i < CCD_V_hhpp.size(); i++)
            helpMatrix[i].mat = CCD_T_old[i].mat * CCD_V_hhpp[i].mat.transpose();
        for (int k = 0; k < A; k++)
        {
            CCD_intermediate_1B[k] = V1B(k,k);
            for (int j = 0; j < A; j++)
            {
                Position pos = CCD_position(k,j);
                if (pos.channel != -1)
                    CCD_intermediate_1B[k] +=
                            +1.0 * helpMatrix[pos.channel].mat(pos.index,pos.index);
            }
        }
        //pppp part
        CCD_intermediate_pppp = CCD_V_pppp;
        //hhhh part
        for(int i = 0; i < CCD_V_hhpp.size(); i++)
            CCD_intermediate_hhhh[i].mat = CCD_V_hhhh[i].mat
	      //+ CCD_V_hhpp[i].mat * CCD_T_old[i].mat.transpose();
	        + CCD_T_old[i].mat * CCD_V_hhpp[i].mat.transpose();
        //difficult part
        for (int channel = 0; channel < CCD_V_hmphpm.size(); channel++)
        {
            for (int Ibra = 0; Ibra < CCD_V_hmphpm[channel].bra.size(); Ibra++)
                for (int Iket = 0; Iket < CCD_V_hmphpm[channel].ket.size(); Iket++)
                {
                    TwoBody_Infinite *State_bra = (TwoBody_Infinite*)TwoBody_States_hmp[ CCD_V_hmphpm[channel].bra[Ibra] ];
                    TwoBody_Infinite *State_ket = (TwoBody_Infinite*)TwoBody_States_hpm[ CCD_V_hmphpm[channel].ket[Iket] ];

                    int i = State_bra->p;
        int a = State_bra->q;
                    int j = State_ket->p;
        int b = State_ket->q;
                    Position pos_pp = CCD_position(a,b);
                    Position pos_hh = CCD_position(i,j);
        int sign=(i-j)*(a-b) > 0 ? 1 : -1;
        if (pos_pp.channel != -1 && pos_hh.channel != -1 && pos_pp.channel==pos_hh.channel)
                       CCD_T_hmphpm[channel].mat(Ibra,Iket)
                               = sign * CCD_T_old[pos_pp.channel].mat(pos_hh.index,pos_pp.index);
                    else
                        CCD_T_hmphpm[channel].mat(Ibra,Iket) = 0.0;

                       //to remove
                    if (pos_pp.channel != -1 && pos_hh.channel != -1)
                   if (abs(CCD_T_hmphpm[channel].mat(Ibra,Iket)
                           - sign * CCD_T_old[pos_pp.channel].mat(pos_hh.index,pos_pp.index)) > 1e-6)
                       cout << "NO\t" << CCD_T_hmphpm[channel].mat(Ibra,Iket)
                            << "\t"<< sign * CCD_T_old[pos_pp.channel].mat(pos_hh.index,pos_pp.index)
                               << endl;
                }
        }
        for (int i = 0; i < CCD_V_hmphpm.size(); i++)
        {
            CCD_intermediate_hpmhpm[i].mat = CCD_V_hpmhpm[i].mat +
                    0.5 * CCD_V_hmphpm[i].mat.transpose() * CCD_T_hmphpm[i].mat;
        }
        helpMatrix = CCD_V_hmphpm;
        for (int i = 0; i < CCD_V_hmphpm.size(); i++)
        {
            helpMatrix[i].mat = CCD_T_hmphpm[i].mat * CCD_intermediate_hpmhpm[i].mat;
        }

        //calculations
        for (int i = 0; i < CCD_V_hhpp.size(); i++)
        {
            CCD_T_hhpp[i].mat = CCD_V_hhpp[i].mat
                    + CCD_T_old[i].mat * CCD_intermediate_pppp[i].mat
                    + CCD_intermediate_hhhh[i].mat * CCD_T_old[i].mat;
        }
        for (int channel = 0; channel < helpMatrix.size(); channel++)
        {
            for (int Ibra = 0; Ibra < helpMatrix[channel].bra.size(); Ibra++)
                for (int Iket = 0; Iket < helpMatrix[channel].ket.size(); Iket++)
                {
                    TwoBody_Infinite *State_bra = (TwoBody_Infinite*)TwoBody_States_hmp[ helpMatrix[channel].bra[Ibra] ];
                    TwoBody_Infinite *State_ket = (TwoBody_Infinite*)TwoBody_States_hpm[ helpMatrix[channel].ket[Iket] ];

                    int i = State_bra->p;
                    int a = State_bra->q;
                    int j = State_ket->p;
                    int b = State_ket->q;
            int sign=(i-j)*(a-b) > 0 ? 1 : -1;
                    Position pos_pp = CCD_position(a,b);
                    Position pos_hh = CCD_position(i,j);//check
            if ((pos_pp.channel != -1) && (pos_hh.channel != - 1)
                && (pos_pp.channel == pos_hh.channel))
                    CCD_T_hhpp[pos_pp.channel].mat(pos_hh.index,pos_pp.index)
                            += sign * helpMatrix[channel].mat(Ibra,Iket);
                }
        }
        //dividing by 1B terms
        for (int i = 0; i < CCD_V_hhpp.size(); i++)
        {
            for (int Ibra = 0; Ibra < CCD_T_old[i].bra.size(); Ibra++)
                for (int Iket = 0; Iket < CCD_T_old[i].ket.size(); Iket++)
                    CCD_T_hhpp[i].mat(Ibra,Iket) /=
                            + CCD_intermediate_1B[TwoBody_States_hh[ CCD_T_old[i].bra[Ibra] ]->p]
                            + CCD_intermediate_1B[TwoBody_States_hh[ CCD_T_old[i].bra[Ibra] ]->q]
                            - CCD_intermediate_1B[TwoBody_States_pp[ CCD_T_old[i].ket[Iket] ]->p]
                            - CCD_intermediate_1B[TwoBody_States_pp[ CCD_T_old[i].ket[Iket] ]->q];
        }

        double sum = 0;
        for(int i = 0; i < CCD_V_hhpp.size(); i++)
        {
            sum += (CCD_V_hhpp[i].mat*(CCD_T_hhpp[i].mat.transpose())).trace();
        }
        CCD_deltaE = sum;
        //cout << iter << "\t" << setprecision(10) << sum << endl;
    }
    while (abs(CCD_deltaE_old - CCD_deltaE) > 1e-9);
}

//void Infinite::CCD_BlockMatricesIntermediates()
//{
//    generateTwoBody_States();
//    CCD_generateBlockMatrices();
//    vector<Channel> CCD_T_old = CCD_T_hhpp;
//    int iter = 0;
//    CCD_deltaE = 0.0;

//    vector<double> CCD_intermediate_1B(numberSP,0.0);
//    vector<Channel> helpMatrix = CCD_V_pppp;

//    vector<Channel> CCD_intermediate_hhhh, CCD_intermediate_pppp;
//    CCD_intermediate_pppp = CCD_V_pppp;
//    CCD_intermediate_hhhh = CCD_V_hhhh;

//    vector<Channel> CCD_intermediate_hmphpm, CCD_intermediate_hpmhpm, CCD_T_hmphpm;
//    CCD_intermediate_hmphpm = CCD_V_hmphpm;
//    CCD_intermediate_hpmhpm = CCD_V_hpmhpm;
//    CCD_T_hmphpm = CCD_V_hmphpm;

//    double CCD_deltaE_old;
//    do
//    {
//        iter++;
//        CCD_deltaE_old = CCD_deltaE;
//        //CCD_T_old = CCD_T_hhpp;
//        double alpha = 0.0;
//        for(int i = 0; i < CCD_V_hhpp.size(); i++)
//            CCD_T_old[i].mat = alpha * CCD_T_old[i].mat + (1-alpha) * CCD_T_hhpp[i].mat;

//        //one body intermediate
//        //particle part
//        for(int i = 0; i < CCD_V_hhpp.size(); i++)
//            helpMatrix[i].mat = CCD_T_old[i].mat.transpose() * CCD_V_hhpp[i].mat;
//        for (int b = A; b < numberSP; b++)
//        {
//            CCD_intermediate_1B[b] = V1B(b,b);
//            for (int d = A; d < numberSP; d++)
//            {
//                Position pos = CCD_position(b,d);
//                if (pos.channel != -1)
//                    CCD_intermediate_1B[b] +=
//                            -1.0 * helpMatrix[pos.channel].mat(pos.index,pos.index);
//            }
//        }
//        //hole part
//        helpMatrix = CCD_V_hhhh;
//        for(int i = 0; i < CCD_V_hhpp.size(); i++)
//            helpMatrix[i].mat = CCD_T_old[i].mat * CCD_V_hhpp[i].mat.transpose();
//        for (int k = 0; k < A; k++)
//        {
//            CCD_intermediate_1B[k] = V1B(k,k);
//            for (int j = 0; j < A; j++)
//            {
//                Position pos = CCD_position(k,j);
//                if (pos.channel != -1)
//                    CCD_intermediate_1B[k] +=
//                            +1.0 * helpMatrix[pos.channel].mat(pos.index,pos.index);
//            }
//        }
//        //pppp part
//        CCD_intermediate_pppp = CCD_V_pppp;
//        //hhhh part
//        for(int i = 0; i < CCD_V_hhpp.size(); i++)
//            CCD_intermediate_hhhh[i].mat = CCD_V_hhhh[i].mat
//                    + CCD_V_hhpp[i].mat * CCD_T_old[i].mat.transpose();
//        //difficult part
//        for (int channel = 0; channel < CCD_V_hmphpm.size(); channel++)
//        {
//            for (int Ibra = 0; Ibra < CCD_V_hmphpm[channel].bra.size(); Ibra++)
//                for (int Iket = 0; Iket < CCD_V_hmphpm[channel].ket.size(); Iket++)
//                {
//                    TwoBody_Infinite *State_bra = (TwoBody_Infinite*)TwoBody_States_hmp[ CCD_V_hmphpm[channel].bra[Ibra] ];
//                    TwoBody_Infinite *State_ket = (TwoBody_Infinite*)TwoBody_States_hpm[ CCD_V_hmphpm[channel].ket[Iket] ];
//                    int sign = +1;
//                    int i = State_bra->p;
//                    int a = State_bra->q;
//                    int j = State_ket->p;
//                    int b = State_ket->q;
//                    if (a > b)
//                        sign *= -1;
//                    if (i > j)
//                        sign *= -1;
//                    Position pos_pp = CCD_position(a,b);
//                    Position pos_hh = CCD_position(i,j);
//                    if (pos_pp.channel != -1)
//                       CCD_T_hmphpm[channel].mat(Ibra,Iket)
//                               = sign * CCD_T_old[pos_pp.channel].mat(pos_hh.index,pos_pp.index);
//                    else
//                        CCD_T_hmphpm[channel].mat(Ibra,Iket) = 0.0;

//                       //to remove
//                    if (pos_pp.channel != -1)
//                   if (abs(CCD_T_hmphpm[channel].mat(Ibra,Iket)
//                           - sign * CCD_T_old[pos_pp.channel].mat(pos_hh.index,pos_pp.index)) > 1e-6)
//                       cout << "NO\t" << CCD_T_hmphpm[channel].mat(Ibra,Iket)
//                            << "\t"<< sign * CCD_T_old[pos_pp.channel].mat(pos_hh.index,pos_pp.index)
//                               << endl;
//                }
//        }
//        for (int i = 0; i < CCD_V_hmphpm.size(); i++)
//        {
//            CCD_intermediate_hpmhpm[i].mat = CCD_V_hpmhpm[i].mat +
//                    0.5 * CCD_V_hmphpm[i].mat.transpose() * CCD_T_hmphpm[i].mat;
//        }
//        helpMatrix = CCD_V_hmphpm;
//        for (int i = 0; i < CCD_V_hmphpm.size(); i++)
//        {
//            helpMatrix[i].mat = CCD_T_hmphpm[i].mat * CCD_intermediate_hpmhpm[i].mat;
//        }

//        //calculations
//        for (int i = 0; i < CCD_V_hhpp.size(); i++)
//        {
//            CCD_T_hhpp[i].mat = CCD_V_hhpp[i].mat
//                    + CCD_T_old[i].mat * CCD_intermediate_pppp[i].mat
//                    + CCD_intermediate_hhhh[i].mat * CCD_T_old[i].mat;
//        }
//        for (int channel = 0; channel < helpMatrix.size(); channel++)
//        {
//            for (int Ibra = 0; Ibra < helpMatrix[channel].bra.size(); Ibra++)
//                for (int Iket = 0; Iket < helpMatrix[channel].ket.size(); Iket++)
//                {
//                    TwoBody_Infinite *State_bra = (TwoBody_Infinite*)TwoBody_States_hmp[ helpMatrix[channel].bra[Ibra] ];
//                    TwoBody_Infinite *State_ket = (TwoBody_Infinite*)TwoBody_States_hpm[ helpMatrix[channel].ket[Iket] ];
//                    int sign = +1;
//                    int i = State_bra->p;
//                    int a = State_bra->q;
//                    int j = State_ket->p;
//                    int b = State_ket->q;
//                    if (i > j)
//                        sign *= -1;
//                    if (a > b)
//                        sign *= -1;
//                    Position pos_pp = CCD_position(a,b);
//                    Position pos_hh = CCD_position(i,j);//check
//                    if ((pos_pp.channel == -1) || (pos_hh.channel == -1)
//                         || (pos_pp.channel != pos_hh.channel))
//                            continue;
//                    CCD_T_hhpp[pos_pp.channel].mat(pos_hh.index,pos_pp.index)
//                            += sign * helpMatrix[channel].mat(Ibra,Iket);
//                }
//        }
//        //dividing by 1B terms
//        for (int i = 0; i < CCD_V_hhpp.size(); i++)
//        {
//            for (int Ibra = 0; Ibra < CCD_T_old[i].bra.size(); Ibra++)
//                for (int Iket = 0; Iket < CCD_T_old[i].ket.size(); Iket++)
//                    CCD_T_hhpp[i].mat(Ibra,Iket) /=
//                            + CCD_intermediate_1B[TwoBody_States_hh[ CCD_T_old[i].bra[Ibra] ]->p]
//                            + CCD_intermediate_1B[TwoBody_States_hh[ CCD_T_old[i].bra[Ibra] ]->q]
//                            - CCD_intermediate_1B[TwoBody_States_pp[ CCD_T_old[i].ket[Iket] ]->p]
//                            - CCD_intermediate_1B[TwoBody_States_pp[ CCD_T_old[i].ket[Iket] ]->q];
//        }

//        double sum = 0;
//        for(int i = 0; i < CCD_V_hhpp.size(); i++)
//        {
//            sum += (CCD_V_hhpp[i].mat*(CCD_T_hhpp[i].mat.transpose())).trace();
//        }
//        CCD_deltaE = sum;
//        cout << iter << "\t" << setprecision(10) << sum << endl;
//    }
//    while (abs(CCD_deltaE_old - CCD_deltaE) > 1e-9);
//}

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
