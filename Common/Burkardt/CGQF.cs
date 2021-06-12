namespace Burkardt
{
    public static class CGQF
    {
        public static void cgqf(int nt, int kind, double alpha, double beta, double a, double b,
            double[] t, double[] wts )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CGQF computes knots and weights of a Gauss quadrature formula.
        //
        //  Discussion:
        //
        //    The user may specify the interval (A,B).
        //
        //    Only simple knots are produced.
        //
        //    Use routine EIQFS to evaluate this quadrature formula.
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
        //    Input, double A, B, the interval endpoints, or
        //    other parameters.
        //
        //    Output, double T[NT], the knots.
        //
        //    Output, double WTS[NT], the weights.
        //
        {
            int i;
            int[] mlt;
            int[] ndx;
            //
            //  Compute the Gauss quadrature formula for default values of A and B.
            //
            CDGQF.cdgqf(nt, kind, alpha, beta, t, wts);
            //
            //  Prepare to scale the quadrature formula to other weight function with 
            //  valid A and B.
            //
            mlt = new int[nt];
            for (i = 0; i < nt; i++)
            {
                mlt[i] = 1;
            }

            ndx = new int[nt];
            for (i = 0; i < nt; i++)
            {
                ndx[i] = i + 1;
            }

            SCQF.scqf(nt, t, mlt, wts, nt, ndx, ref wts, ref t, kind, alpha, beta, a, b);

            return;
        }
    }
}