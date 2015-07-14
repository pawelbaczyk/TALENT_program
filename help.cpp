#include "help.h"
#include <cmath>

int first(int A)
{
    int result = 0;
    for (int i = 0; i < A; i++)
    {
        result += pow(2,i);
    }
    return result;
}

int last(int n, int A)
{
    int result = 0;
    for (int i = 0; i < A; i++)
    {
        result += pow(2,n-1-i);
    }
    return result;
}

int next(int x)
{
   int smallest, ripple, new_smallest, ones;

   if (x == 0) return 0;
   smallest     = (x & -x);
   ripple       = x + smallest;
   new_smallest = (ripple & -ripple);
   ones         = ((new_smallest/smallest) >> 1) - 1;
   return ripple | ones;
}

bool getBit(int x,int n)//get the nth bit of integer x
{
  return ((1<<n)&x)>>n;
}



void gauleg(const double x1, const double x2, vector<double> &x, vector<double> &w)//gauss-legendre integration
{
  const double EPS=1.0e-14;
  double z1,z,xm,xl,pp,p3,p2,p1;
  int n=x.size();
  int m=(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  for (int i=0;i<m;i++)
    {
      z=cos(3.141592654*(i+0.75)/(n+0.5));
      do {
	p1=1.0;
	p2=0.0;
	for (int j=0;j<n;j++)
	  {
	    p3=p2;
	    p2=p1;
	    p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
	  }
	pp=n*(z*p1-p2)/(z*z-1.0);
	z1=z;
	z=z1-p1/pp;
      } while (abs(z-z1) > EPS);
      x[i]=xm-xl*z;
      x[n-1-i]=xm+xl*z;
      w[i]=2.0*xl/((1.0-z*z)*pp*pp);
      w[n-1-i]=w[i];
    }
}

double J1(double x)//bessel_J1
{
  return sin(x)/(x*x)-cos(x)/x;
}
