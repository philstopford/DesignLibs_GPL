using System;
using Burkardt.MatrixNS;

namespace Burkardt.Quadrature
{
    public static class CIQFS
    {
        public static double[] ciqfs(int nt, double[] t, int[] mlt, int nwts, ref int[] ndx, int key,
        int kind, double alpha, double beta, int lo )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIQFS computes some weights of a quadrature formula in the default interval.
        //
        //  Discussion:
        //
        //    This routine computes some or all the weights of a quadrature formula 
        //    for a classical weight function with default values of A and B,
        //    and a given set of knots and multiplicities. 
        //
        //    The weights may be packed into the output array WTS according to a 
        //    user-defined pattern or sequentially. 
        //
        //    The routine will also optionally print knots and weights and a check of 
        //    the moments.
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
        //    Input, int MLT[NT], the multiplicity of the knots.
        //
        //    Input, int NWTS, the number of weights.
        //
        //    Input/output, int NDX[NT],  used to index the output
        //    array WTS.  If KEY = 1, then NDX need not be preset.  For more
        //    details see the comments in CAWIQ.
        //
        //    Input, int KEY, indicates the structure of the WTS
        //    array.  It will normally be set to 1.  For more details see
        //    the comments in CAWIQ.
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
        //    Input, int LO, selects the actions to perform.
        //     > 0, compute and print weights.  Print moments check.
        //     = 0, compute weights.
        //     < 0, compute and print weights.
        //
        //    Output, double CIQFS[NWTS], the weights.
        //
        {
            double[] aj;
            double[] bj;
            int j;
            int jdf;
            int l;
            int m;
            int mex;
            int mop;
            int n;
            int nst;
            double[] w;
            double[] wts;
            double zemu;

            jdf = 0;
            n = 0;
            l = Math.Abs(key);

            for (j = 0; j < nt; j++)
            {
                if (l == 1 || Math.Abs(ndx[j]) != 0)
                {
                    n = n + mlt[j];
                }
            }

            //
            //  N knots when counted according to multiplicity.
            //
            if (nwts < n)
            {
                Console.WriteLine("");
                Console.WriteLine("CIQFS - Fatal error!");
                Console.WriteLine("  NWTS < N.");
                return null;
            }

            m = n + 1;
            mex = 2 + m;
            nst = m / 2;
            //
            //  Get the Jacobi matrix.
            //
            aj = new double[nst];
            bj = new double[nst];

            zemu = Matrix.class_matrix(kind, nst, alpha, beta, ref aj, ref bj);
            //
            //  Call weights routine.
            //
            wts = CAWIQ.cawiq(nt, t, mlt, n, ref ndx, key, nst, ref aj, ref bj, ref jdf, zemu);

            //
            //
            //  Call checking routine.
            //
            if (lo != 0)
            {
                mop = m - 1;

                w = new double[mex];

                CHKQFS.chkqfs(t, wts, mlt, nt, n, ndx, key, ref w, mop, mex, kind,
                    alpha, beta, lo);

            }

            return wts;
        }
    }
}