#include "help.h"
#include "system.h"
#include "pairing.h"
#include "infinite.h"
#include <iomanip>

int main()
{
  Infinite infinite(14,2,0.16,2);//A,g_s,rho,nMax

  //  infinite.printSP_States();
  //  infinite.generateTwoBody_States();
  //   infinite.CCD_BlockMatrices();
  //  infinite.printTwoBody_States();
   infinite.HF_calculateE0();//IMPORTANT;
   infinite.CCD_BlockMatricesIntermediates();cout << endl;
   infinite.CCD_BlockMatricesLadders();

  // infinite.CCD_SparseMatrices();
//   cout<<infinite.CCD_deltaE<<endl;
  // infinite.HF_calculateE0();
  // infinite.CCD_SparseMatrices();
   //  infinite.MBPT();
   //  cout<<infinite.MBPT_deltaE<<endl;
  // cout << (infinite.HF_E0+infinite.MBPT_deltaE)/infinite.A << "\t" << (infinite.HF_E0+infinite.CCD_deltaE)/infinite.A << "\t"
  //      << setprecision(10) << infinite.HF_E0/infinite.A << "\t" << infinite.HF_exact_E0 << endl;

  // Pairing pairing(4,8,0.5);
  // pairing.CCD_OnFlight();
  // pairing.MBPT();
  // cout << pairing.MBPT_deltaE << "\t" << pairing.CCD_deltaE << endl;

  return 0;
}
