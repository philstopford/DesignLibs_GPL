using System;
using Burkardt.Weight;

namespace Burkardt.Quadrature;

public static class CHKQF
{
    public static void chkqf(double[] t, double[] wts, int[] mlt, int nt, int nwts, int[] ndx,
            int key, int mop, int mex, int kind, double alpha, double beta, int lo,
            double a, double b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHKQF computes and prints the moments of a quadrature formula.
        //
        //  Discussion:
        //
        //    The quadrature formula is based on a clasical weight function with 
        //    any valid A, B.
        //
        //    No check can be made for non-classical weight functions.
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
        //    Input, double T[NT], the knots.
        //
        //    Input, double WTS[NWTS], the weights.
        //
        //    Input, int MLT[NT], the multiplicity of the knots.
        //
        //    Input, int NT, the number of knots.
        //
        //    Input, int NWTS, the number of weights.
        //
        //    Input, int NDX[NT], used to index the array WTS.  
        //    If KEY = 1, then NDX need not be preset.  For more details see the 
        //    comments in CAWIQ.
        //
        //    Input, int KEY, indicates the structure of the WTS
        //    array.  It will normally be set to 1.  This will cause the weights to be 
        //    packed sequentially in array WTS.  For more details see the comments 
        //    in CAWIQ.
        //
        //    Input, int MOP, the expected order of precision of the 
        //    quadrature formula.
        //
        //    Input, int MEX, the number of moments required to be 
        //    tested.  Set MEX = 1 and LO < 0 for no moments check.
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
        //    Input, int LO, selects the action to carry out.
        //     > 0, print weights and moment tests.
        //     = 0, print nothing. compute moment test.
        //     < 0, print weights only. don't compute moment tests.
        //
        //    Input, double A, B, the interval endpoints.
        //
    {
        double tmp = 0;

        double[] w = new double[mex];

        PARCHK.parchk(kind, mex, alpha, beta);

        if (lo != 0)
        {
            const int izero = 0;

            Console.WriteLine("");
            Console.WriteLine("  Interpolatory quadrature formula");
            Console.WriteLine("");
            Console.WriteLine("  Type  Interval       Weight function               Name");
            Console.WriteLine("");
            switch (kind)
            {
                case 1:
                    Console.WriteLine("    1    (a,b)              1.0                    Legendre");
                    break;
                case 2:
                    Console.WriteLine("    2    (a,b)      ((b-x)*(x-a))^(-0.5)          Chebyshev Type 1");
                    break;
                case 3:
                    Console.WriteLine("    3    (a,b)      ((b-x)*(x-a))^alpha           Gegenbauer");
                    break;
                case 4:
                    Console.WriteLine("    4    (a,b)    (b-x)^alpha*(x-a)^beta          Jacobi");
                    break;
                case 5:
                    Console.WriteLine("    5   (a,+oo)  (x-a)^alpha*exp(-b*(x-a))      Gen Laguerre");
                    break;
                case 6:
                    Console.WriteLine("    6  (-oo,+oo) |x-a|^alpha*exp(-b*(x-a)^2)  Gen Hermite");
                    break;
                case 7:
                    Console.WriteLine("    7    (a,b)      |x-(a+b)/2.0|^alpha        Exponential");
                    break;
                case 8:
                    Console.WriteLine("    8   (a,+oo)    (x-a)^alpha*(x+b)^beta         Rational");
                    break;
                case 9:
                    Console.WriteLine("    9   (a,b)     (b-x)*(x-a)^(+0.5)         Chebyshev Type 2");
                    break;
            }

            Console.WriteLine("");
            Console.WriteLine("     Parameters   A          " + a + "");
            Console.WriteLine("                  B          " + b + "");
            switch (kind)
            {
                case >= 3 and <= 8:
                    Console.WriteLine("                  alpha      " + alpha + "");
                    break;
            }

            switch (kind)
            {
                case 4:
                case 8:
                    Console.WriteLine("                  beta       " + beta + "");
                    break;
            }

            CHKQFS.chkqfs(t, wts, mlt, nt, nwts, ndx, key, ref w, mop, mex, izero,
                alpha, beta, -Math.Abs(lo));
        }

        switch (lo)
        {
            case >= 0:
            {
                //
                //  Compute the moments in W.
                //
                w = SCMM.scmm(mex, kind, alpha, beta, a, b);

                switch (kind)
                {
                    case 1:
                    case 2:
                    case 3:
                    case 4:
                    case 7:
                    case 9:
                        tmp = (b + a) / 2.0;
                        break;
                    case 5:
                    case 6:
                    case 8:
                        tmp = a;
                        break;
                }

                double[] t2 = new double[nt];

                int i;
                for (i = 0; i < nt; i++)
                {
                    t2[i] = t[i] - tmp;
                }

                const int neg = -1;
                //
                //  Check moments.
                //
                CHKQFS.chkqfs(t2, wts, mlt, nt, nwts, ndx, key, ref w, mop, mex, neg, alpha, beta,
                    lo);
                break;
            }
        }
    }
}