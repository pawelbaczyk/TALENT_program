#include "infinite.h"

Infinite::Infinite(int _A, int _g_s, double _rho, int _nMax) : System(_A)
{
    g_s = _g_s;
    rho = _rho;
    nMax = _nMax;
    generateSP_States(_A);
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
                            SP_States.push_back(new SP_Infinite(kx,ky,kz,0,0));
                            SP_States.push_back(new SP_Infinite(kx,ky,kz,1,0));
                        }
                        if (g_s == 4)
                        {
                            SP_States.push_back(new SP_Infinite(kx,ky,kz,0,1));
                            SP_States.push_back(new SP_Infinite(kx,ky,kz,1,1));
                        }
                    }
}

void Infinite::setRho(double _rho)
{
    rho = _rho;
    generateSP_States(A);
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

void Infinite::generateConfigurations()
{

}

void Infinite::printConfigurations()
{

}

double Infinite::V1B(int a,int b)
{
    if (a == b)
        return SP_States[a]->spEnergy;
    else
        return 0.0;
}

double Infinite::V2B(int p,int q,int r,int s)
{
    return V2B_sym(p, q, r, s) - V2B_sym(p, q, s, r);
}

double Infinite::V2B_sym(int p, int q, int r, int s)
{
    double v = 0.0;
    SP_Infinite* P = (SP_Infinite*)SP_States[p];
    SP_Infinite* Q = (SP_Infinite*)SP_States[q];
    SP_Infinite* R = (SP_Infinite*)SP_States[r];
    SP_Infinite* S = (SP_Infinite*)SP_States[s];
    if ((P->kx + Q->kx == R->kx + S->kx) &&
        (P->ky + Q->ky == R->ky + S->ky) &&
        (P->kz + Q->kz == R->kz + S->kz) &&
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
    for (int i = 0; i < A; i++)
    {
        HF_E0 += V1B(i,i);
        for (int j = 0; j < A; j++)
            HF_E0 += 0.5 * V2B(i,j,i,j);
    }
}

