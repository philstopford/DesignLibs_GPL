namespace Burkardt.Quadrature;

public static class CLIQF
{
    public static double[] cliqf(int nt, double[] t, int kind, double alpha, double beta,
            double a, double b, int lo )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CLIQF computes a classical quadrature formula, with optional printing.
        //
        //  Discussion:
        //
        //    This routine computes all the weights of an interpolatory
        //    quadrature formula with
        //    1. only simple knots and
        //    2. a classical weight function with any valid A and B, and
        //    3. optionally prints the knots and weights and a check of the moments.
        //
        //    To evaluate this quadrature formula for a given function F,
        //    call routine EIQFS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 January 2010
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
        //    Input, int LO, indicates what is to be done.
        //    > 0, compute and print weights and moments check.
        //    = 0, compute weights.
        //    < 0, compute and print weights.
        //
        //    Output, double CLIQF[NT], the weights.
        //
    {
        int i;

        const int key = 1;
        int[] mlt = new int[nt];
        for (i = 0; i < nt; i++)
        {
            mlt[i] = 1;
        }

        int[] ndx = new int[nt];

        double[] wts = CIQF.ciqf(nt, t, mlt, nt, ref ndx, key, kind, alpha, beta, a, b, lo);

        return wts;
    }
}