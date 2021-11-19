using System;

namespace Burkardt.Quadrature;

public static class CEGQF
{
    public static double cegqf(int nt, int kind, double alpha, double beta, double a, double b,
            Func < double, int, double > f )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CEGQF computes a quadrature formula and applies it to a function.
        //
        //  Discussion:
        //
        //    The user chooses the quadrature formula to be used, as well as the
        //    interval (A,B) in which it is applied.
        //
        //    Note that the knots and weights of the quadrature formula are not
        //    returned to the user.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 February 2010
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
        //    Input, int KIND, the rule.
        //    1, Legendre,             (a,b)       1.0
        //    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
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
        //    at X.  The value I will always be 0.  The value X will always be a knot.
        //
        //    Output, double CEGQF, the value of the quadrature formula 
        //    applied to F.
        //
    {
        const int lo = 0;
        double[] t = new double[nt];
        double[] wts = new double[nt];

        CGQF.cgqf(nt, kind, alpha, beta, a, b, lo, ref t, ref wts);
        //
        //  Evaluate the quadrature sum.
        //
        double qfsum = EIQFS.eiqfs(nt, t, wts, f);

        return qfsum;
    }
}