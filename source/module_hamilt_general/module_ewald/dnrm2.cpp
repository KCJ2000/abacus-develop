#include "dnrm2.h"

#include<math.h>
#include <iostream>

double dnrm2(const int n, const double *x, const int incx)
{
    // compute Euclidean length (12 norm) of std::vector x,
    // Actually the input is an double array type instead of std::vector type
    // The result is L2 norm
    // int n is the maximum of the index
    // int incx is the interval of two index
    if (n < 0 || incx <= 0)
    {
        std::cerr << "\n error in dnrm2, n < 0 or incx <= 0, ";
        return 0;
    }
    int len=sizeof(x)/sizeof(x[0]);
    if (n-incx > len)
    {
        std::cerr << "\n error in dnrm2, the index is out of the array";
        return 0;
    }
    if (n == 0)
    {
        return 0;
    }
    // all of judgement statements above are insurance code making sure that
    // the index must be up to the mustard.
    // Some C++ interpreters are not sensitive to index bug.
    // So this part is necessary.

    double norm2=0.0;
    for (int ix=0; ix<n; ix+=incx)
    {
        norm2 += x[ix] * x[ix];
    }

    return sqrt(norm2);
}