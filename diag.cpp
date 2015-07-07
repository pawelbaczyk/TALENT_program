#include"diag.h"
void generate_H(MatrixXd&H,double g)
{
  H.resize(6,6);
  for(int i=0;i<6;i++)
    for(int j=0;j<6;j++)
      {
	if(i!=j)
	  H(i,j)=-0.5*g;
	if(i+j==6)
	  H(i,j)=0;
	if(i==j)
	  H(i,j)=2*(i+1)-g;
      }
}
void diag(MatrixXd &H,double&E_gs)
{
  SelfAdjointEigenSolver<MatrixXd> solver(H);
  if(solver.info()!=Success) abort();
  E_gs=solver.eigenvalues()[0];
}
