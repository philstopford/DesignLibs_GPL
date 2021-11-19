using System;
using Burkardt.Interpolation;

namespace Burkardt.Quadrature;

public static class CIQF
{
    public static double[] ciqf(int nt, double[] t, int[] mlt, int nwts, ref int[] ndx, int key,
            int kind, double alpha, double beta, double a, double b, int lo )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIQF computes weights for a classical weight function and any interval.
        //
        //  Discussion:
        //
        //    This routine compute somes or all the weights of a quadrature formula
        //    for a classical weight function with any valid A, B and a given set of 
        //    knots and multiplicities.  
        //
        //    The weights may be packed into the output array WTS according to a 
        //    user-defined pattern or sequentially. 
        //
        //    The routine will also optionally print knots and weights and a check 
        //    of the moments.
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
        //    Input/output, int NDX[NT], used to index the output
        //    array WTS.  If KEY = 1, then NDX need not be preset.  For more
        //    details see the comments in CAWIQ.
        //
        //    Input, int KEY, indicates the structure of the WTS
        //    array.  It will normally be set to 1.  This will cause the weights to be 
        //    packed sequentially in array WTS.  For more details see the comments 
        //    in CAWIQ.
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
        //    Input, int LO, selects the actions to perform.
        //     > 0, compute and print weights.  Print moments check.
        //     = 0, compute weights.
        //     < 0, compute and print weights.
        //
        //    Output, double CIQF[NWTS], the weights.
        //
    {
        int j;

        int m = 1;
        int l = Math.Abs(key);

        for (j = 0; j < nt; j++)
        {
            if (l == 1 || Math.Abs(ndx[j]) != 0)
            {
                m += mlt[j];
            }
        }

        if (nwts + 1 < m)
        {
            Console.WriteLine("");
            Console.WriteLine("CIQF - Fatal error!");
            Console.WriteLine("  NWTS + 1 < M.");
            return null;
        }

        int mex = 2 + m;
        //
        //  Scale the knots to default A, B.
        //
        double[] st = SCT.sct(nt, t, kind, a, b);
        //
        //  Compute the weights.
        //
        int lu = 0;

        double[] wts = CIQFS.ciqfs(nt, st, mlt, nwts, ref ndx, key, kind, alpha, beta, lu);
        //
        //  Don't scale user's knots - only scale weights.
        //
        SCQF.scqf(nt, st, mlt, wts, nwts, ndx, ref wts, ref st, kind, alpha, beta, a, b);

        if (lo != 0)
        {
            int mop = m - 1;

            CHKQF.chkqf(t, wts, mlt, nt, nwts, ndx, key, mop, mex, kind,
                alpha, beta, lo, a, b);
        }

        return wts;
    }
}