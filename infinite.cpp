#include "infinite.h"

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
                            SP_States.push_back(new SP_Infinite(kx,ky,kz,0,0));
                            SP_States.push_back(new SP_Infinite(kx,ky,kz,1,0));
                        }
                        if (g_s == 4)
                        {
                            SP_States.push_back(new SP_Infinite(kx,ky,kz,0,1));
                            SP_States.push_back(new SP_Infinite(kx,ky,kz,1,1));
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

void Infinite::CCD_BlockMatrices()
{
    vector<Channel> CCD_V_hhhh, CCD_V_hhpp;

    double nx, ny, nz;
    int S, T = -2;
    for (int i = 0; i < states2B_hh.size(); i++)
    {
        State2B currentState = states2B_hh[i];
        if ((currentState.nx != nx) || (currentState.ny != ny) || (currentState.nz != nz) ||
                (currentState.S != S) || (currentState.T != T))
        {
            nx = currentState.nx;
            ny = currentState.ny;
            nz = currentState.nz;
            S = currentState.S;
            T = currentState.T;
            CCD_V_hhhh.push_back(Channel(nx,ny,nz,S,T));
        }
        CCD_V_hhhh[CCD_V_hhhh.size()-1].bra.push_back(i);
        CCD_V_hhhh[CCD_V_hhhh.size()-1].ket.push_back(i);
    }


    for (int i = 0; i < states2B_hh.size(); i++)
    {
        bool good = 0;
        State2B currentState = states2B_hh[i];
        if ((currentState.nx != nx) || (currentState.ny != ny) || (currentState.nz != nz) ||
                (currentState.S != S) || (currentState.T != T))
        {
            nx = currentState.nx;
            ny = currentState.ny;
            nz = currentState.nz;
            S = currentState.S;
            T = currentState.T;
            for (int j = 0; j < states2B_pp.size(); j++)
            {
                State state_pp = states2B_pp[j];
                if ((state_pp.nx == nx) || (state_pp.ny == ny) || (state_pp.nz == nz) ||
                        (state_pp.S == S) || (state_pp.T == T))
                {
                    if (good == 0)
                    {
                        CCD_V_hhpp.push_back(Channel(nx,ny,nz,S,T));
                        good = 1;
                    }
                    CCD_V_hhpp[CCD_V_hhpp.size()-1].ket.push_back(j);
                }
            }
        }
        CCD_V_hhpp[CCD_V_hhhh.size()-1].bra.push_back(i);
    }





//    for (int i = 0; i < CCD_V_hhhh.size(); i++)
//    {
//        Channel currentChannel = CCD_V_hhhh[i];
//        for(int s1 = 0; s1 < currentChannel.states.size(); s1++)
//            for(int s2 = 0; s2 < currentChannel.states.size(); s2++)
//                currentChannel.m(s1,s2) = V2B(s1,s1,s2,s2);
//    }

//    for (int i = 0; i < CCD_V_hhhh.size(); i++)
//        for (int j = 0; j < CCD_t_hhpp.size(); j++)
//        {

//        }
}

