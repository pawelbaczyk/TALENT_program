#ifndef HELP_H
#define HELP_H
#include<vector>
using namespace std;

int first(int A);

int last(int n, int A);

int next(int x);

bool getBit(int x,int n);

double J1(double x);//bessel_J1


void gauleg(const double x1, const double x2, vector<double> &x, vector<double> &w);//gauss-legendre integration

#endif
