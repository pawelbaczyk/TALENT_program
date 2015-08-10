#include "help.h"
#include "system.h"
#include "pairing.h"
#include "infinite.h"
#include <iomanip>
#include <time.h>

void TA_vs_HF(int nA, int nMax, double (Infinite::*func)(), double (Infinite::*ref)() )
{
    Infinite infinite(2,2,0.16,nMax);//A,g_s,rho,nMax

    cout << "PARAMETERS" << endl;
    cout << "g_s\tk_F\trho\tnMax\tnA\n";
    cout << infinite.g_s << "\t" << infinite.k_F << "\t" << infinite.rho << "\t"
         << nMax << "\t" << nA << endl;
    cout << "\n\n";

    infinite.generateMagicNumbers();
    cout << "\n\n";

    cout << "TWISTED AVERAGE" << endl;
    cout << "Shell\t1-HF/INF\t1-TA/INF\t1-SP/INF\tSP_x\tSP_y\tSP_z\tTime" << endl;

    for (int i = 1; i < infinite.magicNumbers.size(); i++)
    {
        clock_t tStart = clock();
        infinite.setA(infinite.magicNumbers[i]);
        double value = (infinite.*func)();
        double refValue = (infinite.*ref)();
        infinite.TA_calculateE0(nA,func,ref);
        cout << infinite.magicNumbers[i] << "\t"
             << abs(1.0 - value/infinite.A / refValue) << "\t"
             << abs(1.0 - (infinite.TA_deltaE)/infinite.A / refValue) << "\t"
             << abs(1.0 - (infinite.SP_deltaE)/infinite.A / refValue) << "\t"
             << infinite.SP_x << "\t" << infinite.SP_y << "\t" << infinite.SP_z << "\t"
             << (double)(clock() - tStart)/CLOCKS_PER_SEC << "\n";
    }
}

int main()
{
    int nA, nMax;
    char mode;
    cin >> nA >> nMax >> mode;
    if (mode == 'E')
    {
        cout << "Calculating total energy E0" << endl;
        TA_vs_HF(nA, nMax, &Infinite::HF_calculateE0, &Infinite::HF_exactE0);
    }
    else if (mode == 'T')
    {
        cout << "Calculating kinetic energy T0" << endl;
        TA_vs_HF(nA, nMax, &Infinite::HF_calculateT0, &Infinite::HF_exactT0);
    }
    else if (mode == 'V')
    {
        cout << "Calculating potential energy V0" << endl;
        TA_vs_HF(nA, nMax, &Infinite::HF_calculateV0, &Infinite::HF_exactV0);
    }
    else
        cout << "Wrong input!" << endl;

//    infinite.HF_cal_exact_E0();
//    infinite.HF_calculateE0();
//    //infinite.TA_calculateE0(10);

//    cout << (infinite.TA_deltaE)/infinite.A << "\t"
//         << infinite.HF_E0/infinite.A << "\t"
//         << infinite.HF_exact_E0 << endl
//         << 1.0 - (infinite.TA_deltaE)/infinite.A / infinite.HF_exact_E0 << "\t"
//         << 1.0 - infinite.HF_E0/infinite.A / infinite.HF_exact_E0 << endl;
    //infinite.CCD_BlockMatricesIntermediates();


    //infinite.MBPT();
    //  infinite.HF_cal_exact_E0();
    //  cout << "MBPT\t" << "CCD\t" << "TA\t" << "HF\t" << "HFinf" << endl;
    //  cout << (infinite.HF_E0+infinite.MBPT_deltaE)/infinite.A << "\t"
    //       << (infinite.HF_E0+infinite.CCD_deltaE)/infinite.A << "\t"
    //       << (infinite.TA_deltaE)/infinite.A << "\t"
    //       << setprecision(10) << infinite.HF_E0/infinite.A << "\t"
    //       << infinite.HF_exact_E0 << endl;

    //Pairing pairing(4,8,0.5);
    // pairing.CCD_OnFlight();
    //pairing.MBPT();
    // cout << pairing.MBPT_deltaE << "\t" << pairing.CCD_deltaE << endl;

    return 0;
}
