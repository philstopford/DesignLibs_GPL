﻿using System;

namespace Burkardt.Quadrature;

public static class CEIQF
{
    public static double ceiqf(int nt, double[] t, int[] mlt, int kind, double alpha,
            double beta, double a, double b, Func <double, int, double > f )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CEIQF constructs and applies a quadrature formula based on user knots.
        //
        //  Discussion:
        //
        //    The knots may have multiplicity.  The quadrature interval is over
        //    any valid A, B.  A classical weight function is selected by the user.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 January 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Sylvan Elhay, Jaroslav Kautsky,
        //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
        //    Interpolatory Quadrature,
        //    ACM Transactions on Mathematical Software,
        //    Volume 13, Number 4, December 1987, pages 399-415.
        //
        //  Parameters:
        //
        //    Input, int NT, the number of knots.
        //
        //    Input, double T[NT], the knots.
        //
        //    Input, int MLT[NT], the multiplicity of the knots.
        //
        //    Input, int KIND, the rule.
        //    1, Legendre,             (a,b)       1.0
        //    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
        //    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
        //    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
        //    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
        //    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
        //    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
        //    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
        //    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
        //
        //    Input, double ALPHA, the value of Alpha, if needed.
        //
        //    Input, double BETA, the value of Beta, if needed.
        //
        //    Input, double A, B, the interval endpoints.
        //
        //    Input, double F ( double X, int I ), the name of a routine which
        //    evaluates the function and some of its derivatives.  The routine
        //    must return in F the value of the I-th derivative of the function
        //    at X.  The highest value of I will be the maximum value in MLT minus
        //    one.  The value X will always be a knot.
        //
        //    Output, double CEIQF, the value of the quadrature formula 
        //    applied to F.
        //
    {
        int i;

        const int lu = 0;
        int n = 0;
        for (i = 0; i < nt; i++)
        {
            n += mlt[i];
        }

        const int key = 1;
        int[] ndx = new int[nt];

        double[] wts = CIQF.ciqf(nt, t, mlt, n, ref ndx, key, kind, alpha, beta, a, b, lu);

        double qfsum = EIQF.eiqf(nt, t, mlt, wts, n, ndx, key, f);

        return qfsum;
    }
}