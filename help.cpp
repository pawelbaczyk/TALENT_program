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
