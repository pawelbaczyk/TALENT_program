#include "help.h"
#include "system.h"
#include "pairing.h"
#include "infinite.h"

int main()
{
  Infinite infinite(2,2,0.10,1);//A,g_s,rho,nMax
 // infinite.printSP_States();



  for (double rho = 0.002; rho <= 0.16; rho += 0.002)
    {
      infinite.setRho(rho);
      infinite.HF_calculateE0();
      infinite.HF_cal_exact_E0();
      cout << infinite.rho << "\t" << infinite.HF_E0 / infinite.A << "\t" << infinite.HF_exact_E0 << endl;
    }



  //infinite.CCD_generateMatrices();
 // infinite.CCD_calculateTau();
 // cout << infinite.CCD_deltaE << endl;
  //cout << infinite.CCD_V_ph << endl;
 // infinite.MBPT();
  //cout << infinite.MBPT_deltaE<< endl;

//  Pairing pairing(4,8,1);
//  pairing.CCD_generateMatrices();
//  pairing.CCD_calculateTau();
//  cout << pairing.CCD_V_ph << endl;

  //    Pairing pairing(4,8,0.4);//A, numberSP, g

  //    cout << "g\t" << "diag\t" << "f_CCD\t" << "2_CCD" << endl;

  //    for (double g = -1.0; g <= 1.0; g += 0.1)
  //    {
  //        pairing.g = g;
  //        pairing.CCD_generateMatrices();
  //        pairing.CCD_calculateTau();
  //        pairing.MBPT();

  //        pairing.diagonalization();

  //        pairing.GF_generateMatrices();
  //        pairing.GF_diag();

  //        cout << g << "\t" ;//<< pairing.MBPT_deltaE + 2 - pairing.g << "\t";
  //        //cout << pairing.diag_E_GS << "\t";
  //        cout << pairing.CCD_deltaE + pairing.H(0,0) - pairing.diag_E_GS<< "\t";

  //        //cout << pairing.GF_E_GS <<endl;
  //    }
  return 0;
}
