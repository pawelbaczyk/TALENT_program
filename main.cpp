#include "help.h"
#include "system.h"
#include "pairing.h"
#include "infinite.h"
#include <iomanip>

int main()
{
  Infinite infinite(2,2,0.16,1);//A,g_s,rho,nMax

<<<<<<< HEAD

=======
  infinite.HF_calculateE0();//IMPORTANT;

  
  infinite.CCD_BlockMatricesIntermediates();
  //  infinite.CCD_BlockMatricesLadders();

  //  infinite.CCD_SparseMatrices();
  //  cout<<infinite.CCD_deltaE<<endl;
  // infinite.HF_calculateE0();
  // infinite.CCD_SparseMatrices();
  //  infinite.MBPT();
  //  cout<<infinite.MBPT_deltaE<<endl;
  // cout << (infinite.HF_E0+infinite.MBPT_deltaE)/infinite.A << "\t" << (infinite.HF_E0+infinite.CCD_deltaE)/infinite.A << "\t"
  //      << setprecision(10) << infinite.HF_E0/infinite.A << "\t" << infinite.HF_exact_E0 << endl;

  Pairing pairing(4,8,0.5);
  // pairing.CCD_OnFlight();
  //pairing.MBPT();
  // cout << pairing.MBPT_deltaE << "\t" << pairing.CCD_deltaE << endl;
>>>>>>> f08fd8fbe494c8727c7a5f4caa7170e57e376196

  return 0;
}
