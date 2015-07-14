#include "help.h"
#include "system.h"
#include "pairing.h"
#include "infinite.h"

int main()
{
    Infinite infinite(14,2,0.10,3);//A,g_s,rho,nMax
    //infinite.printSP_States();

    for (double rho = 0.002; rho <= 0.16; rho += 0.002)
    {
        infinite.setRho(rho);
        infinite.HF_calculateE0();
	infinite.HF_cal_exact_E0();
        cout << infinite.rho << "\t" << infinite.HF_E0 / infinite.A << "\t" <<infinite.HF_exact_E0<< endl;
    }

    return 0;
}
