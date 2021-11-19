namespace Burkardt.Quadrature;

public static class CLIQFS
{
    public static double[] cliqfs(int nt, double[] t, int kind, double alpha, double beta,
            int lo )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CLIQFS computes the weights of a quadrature formula in the default interval.
        //
        //  Discussion:
        //
        //    This routine computes the weights of an interpolatory quadrature formula
        //    with a classical weight function, in the default interval A, B,
        //    using only simple knots.
        //
        //    It can optionally print knots and weights and a check of the moments.
        //
        //    To evaluate a quadrature computed by CLIQFS, call EIQFS.
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
        //    Input, int LO, chooses the printing option.
        //     > 0, compute weights, print them, print the moment check results.
        //     0, compute weights.
        //     < 0, compute weights and print them.
        //
        //    Output, double CLIQFS[NT], the weights.
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

        double[] wts = CIQFS.ciqfs(nt, t, mlt, nt, ref ndx, key, kind, alpha, beta, lo);

        return wts;
    }
}