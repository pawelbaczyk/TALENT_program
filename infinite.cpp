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
        return SP_States[a]->spEnergy;
    else
        return 0.0;
}

double Infinite::V2B(int p,int q,int r,int s)
{
//    char key[100];
//    int a, b, c, d;
//    int sign = +1;
//    if (p == q)
//        return 0.0;
//    else if (p < q)
//    {
//        a = p;
//        b = q;
//    }
//    else
//    {
//        a = q;
//        b = p;
//        sign *= -1;
//    }
//    if (r == s)
//        return 0.0;
//    else if (r < s)
//    {
//        c = r;
//        d = s;
//    }
//    else
//    {
//        c = s;
//        d = r;
//        sign *= -1;
//    }
//    if (p+q <= r+s)
//        sprintf(key,"%.4d%.4d%.4d%.4d",a,b,c,d);
//    else
//        sprintf(key,"%.4d%.4d%.4d%.4d",c,d,a,b);
//    map<string,double>::iterator it;
//    it = map_V2B.find(string(key));
//    if (it != map_V2B.end())
//        return (it->second)*sign;
//    else
//        return 0.0;
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
    for (int i = 0; i < A; i++)
    {
        HF_E0 += V1B(i,i);
        for (int j = 0; j < A; j++)
            HF_E0 += 0.5 * V2B(i,j,i,j);
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

void Infinite::map_generateV2B()
{
    if (numberSP > 10000)
    {
        cout << "Too many single particle states! Change map_generateV2B function" << endl;
        exit(1);
    }
    map_V2B.clear();
    for (int a = 0; a < numberSP; a++)
        for (int b = a+1; b < numberSP; b++)
            for (int c = 0; c < numberSP; c++)
                for (int d = max(a+b-c,c+1); d < numberSP; d++)
                {
                    char key[100];
                    sprintf(key,"%.4d%.4d%.4d%.4d",a,b,c,d);
                    SP_Infinite* A = (SP_Infinite*)SP_States[a];
                    SP_Infinite* B = (SP_Infinite*)SP_States[b];
                    SP_Infinite* C = (SP_Infinite*)SP_States[c];
                    SP_Infinite* D = (SP_Infinite*)SP_States[d];
                    if ((A->kx + B->kx == C->kx + D->kx) &&
                        (A->ky + B->ky == C->ky + D->ky) &&
                        (A->kz + B->kz == C->kz + D->kz) &&
                        (A->spin + B->spin == C->spin + D->spin) &&
                        (A->isospin + B->isospin == C->isospin + D->isospin))
                        map_V2B.insert(pair<string,double>(string(key),V2B_sym(a,b,c,d)-V2B_sym(a,b,d,c)));
                }
}


void Infinite::CCD_calculateDeltaE()
{
    //CCD_Tau = CCD_V_ph;//first guess
    map_t.clear();
    map<string,double> map_t_prev;
    for (int i = 0; i < A; i++)
        for (int j = i+1; j < A; j++)
            for (int a = A; a < numberSP; a++)
                for (int b = a+1; b < numberSP; b++)
                {
                    char key[100];
                    sprintf(key,"%.4d%.4d%.4d%.4d",i,j,a,b);
                    SP_Infinite* A = (SP_Infinite*)SP_States[a];
                    SP_Infinite* B = (SP_Infinite*)SP_States[b];
                    SP_Infinite* I = (SP_Infinite*)SP_States[i];
                    SP_Infinite* J = (SP_Infinite*)SP_States[j];
                    if ((A->kx + B->kx == I->kx + J->kx) &&
                        (A->ky + B->ky == I->ky + J->ky) &&
                        (A->kz + B->kz == I->kz + J->kz) &&
                        (A->spin + B->spin == I->spin + J->spin) &&
                        (A->isospin + B->isospin == I->isospin + J->isospin))
                    {
                        double epsilon = I->spEnergy + J->spEnergy - A->spEnergy - B->spEnergy;
                        map_t.insert(pair<string,double>(string(key),V2B(i,j,a,b)/epsilon));
                    }
                }
    int iter = 0;
    double CCD_deltaE_prev = 0.0;
    do
    {
        iter++;
        map_t_prev = map_t;
        CCD_deltaE_prev = CCD_deltaE;
        for (map<string,double>::iterator it_t = map_t.begin(); it_t != map_t.end(); ++it_t)
        {
            int a, b, i, j;
            sscanf((it_t->first).c_str(),"%4d%4d%4d%4d",&i,&j,&a,&b);
            it_t->second = V2B(a,b,i,j);//t is tau

            double help = 0.0;
            for (int c = A; c < numberSP; c++)
                for (int d = c+1; d < numberSP; d++)
                    help += V2B(a,b,c,d) * CCD_t(map_t_prev,i,j,c,d);
            //for (int k = 0; k < A; k++)
                //for (int l = k+1; l < A; l++)
                    //help += V2B(k,l,i,j) * CCD_t(map_t_prev,k,l,a,b);
//                        for (int k = 0; k < A; k++)
//                            for (int c = A; c < numberSP; c++)
//                                help += V2B(k,b,c,j) * CCD_t(i,k,a,c)
//                                      - V2B(k,b,c,i) * CCD_t(j,k,a,c)
//                                      - V2B(k,a,c,j) * CCD_t(i,k,b,c)
//                                      + V2B(k,a,c,i) * CCD_t(j,k,b,c);
//                        for (int k = 0; k < A; k++)
//                            for (int l = 0; l < A; l++)
//                                for (int c = A; c < numberSP; c++)
//                                    for (int d = A; c < numberSP; c++)
//                                    {
//                                        help += 0.25 * V2B(k,l,c,d) * CCD_t(i,j,c,d) * CCD_t(k,l,a,b);
//                                        help += V2B(k,l,c,d) * CCD_t(i,k,a,c) * CCD_t(j,l,b,d)
//                                              - V2B(k,l,c,d) * CCD_t(j,k,a,c) * CCD_t(i,l,b,d);
//                                        help += -0.5 * V2B(k,l,c,d) * CCD_t(i,k,d,c) * CCD_t(l,j,a,b)
//                                              +  0.5 * V2B(k,l,c,d) * CCD_t(j,k,d,c) * CCD_t(l,i,a,b);
//                                        help += -0.5 * V2B(k,l,c,d) * CCD_t(l,k,a,c) * CCD_t(i,j,d,b)
//                                              +  0.5 * V2B(k,l,c,d) * CCD_t(l,k,b,c) * CCD_t(i,j,d,a);
//                                    }
            it_t->second += help;
        }
        CCD_deltaE = 0.0;
        for (map<string,double>::iterator it_t = map_t.begin(); it_t != map_t.end(); ++it_t)
        {
            int a, b, i, j;
            sscanf((it_t->first).c_str(),"%4d%4d%4d%4d",&i,&j,&a,&b);
            SP_Infinite* A = (SP_Infinite*)SP_States[a];
            SP_Infinite* B = (SP_Infinite*)SP_States[b];
            SP_Infinite* I = (SP_Infinite*)SP_States[i];
            SP_Infinite* J = (SP_Infinite*)SP_States[j];
            double epsilon = I->spEnergy + J->spEnergy - A->spEnergy - B->spEnergy;
            it_t->second *= 1/epsilon;//t is real t
            CCD_deltaE += it_t->second * V2B(a,b,i,j);

        }
        cout << CCD_deltaE << endl;
    }
    while (abs(CCD_deltaE_prev - CCD_deltaE) > 1e-8);
}

double Infinite::CCD_t(map<string,double> &_map, int p, int q, int r, int s)
{
    char key[100];
    int a, b, c, d;
    int sign = +1;
    if (p == q)
        return 0.0;
    else if (p < q)
    {
        a = p;
        b = q;
    }
    else
    {
        a = q;
        b = p;
        sign *= -1;
    }
    if (r == s)
        return 0.0;
    else if (r < s)
    {
        c = r;
        d = s;
    }
    else
    {
        c = s;
        d = r;
        sign *= -1;
    }
    if (p+q <= r+s)
        sprintf(key,"%.4d%.4d%.4d%.4d",a,b,c,d);
    else
        sprintf(key,"%.4d%.4d%.4d%.4d",c,d,a,b);
    map<string,double>::iterator it;
    it = _map.find(string(key));
    if (it != _map.end())
        return (it->second)*sign;
    else
        return 0.0;
}
