#ifndef DIAG_H
#define DIAG_H
#include"Eigen/Dense"

using namespace std;
using namespace Eigen;

void generate_H(MatrixXd&H,double g);
void diag(MatrixXd &H,double&E_gs,VectorXd& coeff);
#endif
