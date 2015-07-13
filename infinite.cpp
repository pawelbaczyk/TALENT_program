#include "infinite.h"

Infinite::Infinite(int _A, int _numberSP) : System(_A, _numberSP)
{
    generateSP_States(_A);
}

void Infinite::generateSP_States(int A)
{
    int nMax = 2;
    double L = 1;
    for (int N = 0; N*N <= 3*nMax*nMax; N++)
        for (int i = -nMax; i <= nMax; i++)
            for (int j = -nMax; j <= nMax; j++)
                for (int k = -nMax; k <= nMax; k++)
                    if (i*i + j*j + k*k == N*N)
                    {
                        double kx = 2 * M_PI * i / L;
                        double ky = 2 * M_PI * j / L;
                        double kz = 2 * M_PI * k / L;
                        cout << N << "\t" << i << "\t" << j << "\t" << k << endl;
                        SP_States.push_back(new SP_Infinite(kx,ky,kz,0,0));
                        SP_States.push_back(new SP_Infinite(kx,ky,kz,0,1));
                        SP_States.push_back(new SP_Infinite(kx,ky,kz,1,0));
                        SP_States.push_back(new SP_Infinite(kx,ky,kz,1,1));
                    }
}

void Infinite::generateConfigurations()
{

}

void Infinite::printConfigurations()
{

}

double Infinite::V1B(int,int)
{
    return 0.0;
}

double Infinite::V2B(int,int,int,int)
{
    return 0.0;
}
