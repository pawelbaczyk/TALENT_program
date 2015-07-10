#include "help.h"
#include "system.h"
#include "pairing.h"

int main()
{
    Pairing pairing(4,8,0.4);//A, numberSP, g

    cout << "g\t" << "MBPT\t" << "diag\t" << "CCD\t" << "GF" << endl;

    for (double g = -1.0; g <= 1.0; g += 0.1)
    {
        pairing.g = g;
        pairing.CCD_generateMatrices();
        pairing.CCD_calculateTau();
        pairing.MBPT();

        pairing.diagonalization();

        pairing.GF_generateMatrices();
        pairing.GF_diag();

        cout << g << "\t" << pairing.MBPT_deltaE + 2 - pairing.g << "\t";
        cout << pairing.diag_E_GS << "\t";
        cout << pairing.CCD_deltaE + 2 - pairing.g<< "\t";
        cout << pairing.GF_E_GS <<endl;
    }
    return 0;
}
