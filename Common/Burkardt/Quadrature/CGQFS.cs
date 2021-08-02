using System;

namespace Burkardt
{
    public static class CGQFS
    {
        public static void cgqfs(int nt, int kind, double alpha, double beta, int lo, ref double[] t,
        ref double[] wts )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CGQFS computes knots and weights of a Gauss quadrature formula.
        //
        //  Discussion:
        //
        //    This routine computes the knots and weights of a Gauss quadrature
        //    formula with:
        //
        //    * a classical weight function with default values for A and B;
        //    * only simple knots
        //    * optionally print knots and weights and a check of the moments
        //
        //    Use routine EIQFS to evaluate a quadrature formula computed by 
        //    this routine.
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
        //    Input, int LO, selects the action.
        //    > 0, compute and print knots and weights.  Print moments check.
        //    = 0, compute knots and weights.
        //    < 0, compute and print knots and weights.
        //
        //    Output, double T[NT], the knots.
        //
        //    Output, double WTS[NT], the weights.
        //
        {
            int i;
            int key;
            int m;
            int mex;
            int[] mlt;
            int mmex;
            int mop;
            int[] ndx;
            double[] w;
            //
            //  Check there is enough workfield and assign workfield
            //
            key = 1;
            mop = 2 * nt;
            m = mop + 1;
            mex = m + 2;
            mmex = Math.Max(mex, 1);
            //
            //  Compute the Gauss quadrature formula for default values of A and B.
            //
            CDGQF.cdgqf(nt, kind, alpha, beta, ref t, ref wts);
            //
            //  Exit if no print required.
            //
            if (lo != 0)
            {
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

                w = new double[mmex];

                CHKQFS.chkqfs(t, wts, mlt, nt, nt, ndx, key, ref w, mop, mmex, kind, alpha,
                    beta, lo);
            }
        }
    }
}