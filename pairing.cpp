#include "pairing.h"

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
