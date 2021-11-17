using System;
using Burkardt.Types;

namespace Burkardt.Weight;

public static class SCMM
{
    public static double[] scmm(int m, int kind, double alpha, double beta, double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SCMM computes moments of a classical weight function scaled to [A,B].
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
        //    Input, int M, the number of moments.
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
        //    Output, double W(M), the scaled moments.
        //
    {
        double al;
        double be;
        int i;
        double p = 0;
        double q = 0;

        double temp = typeMethods.r8_epsilon();

        switch (kind)
        {
            case 1:
            {
                al = 0.0;
                be = 0.0;
                if (Math.Abs(b - a) <= temp)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SCMM - Fatal error!");
                    Console.WriteLine("  B - A too small!");
                    return null;
                }

                q = (b - a) / 2.0;
                p = Math.Pow(q, al + be + 1.0);
                break;
            }
            case 2:
            {
                al = -0.5;
                be = -0.5;
                if (Math.Abs(b - a) <= temp)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SCMM - Fatal error!");
                    Console.WriteLine("  B - A too small!");
                    return null;
                }

                q = (b - a) / 2.0;
                p = Math.Pow(q, al + be + 1.0);
                break;
            }
            case 3:
            {
                al = alpha;
                be = alpha;
                if (Math.Abs(b - a) <= temp)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SCMM - Fatal error!");
                    Console.WriteLine("  B - A too small!");
                    return null;
                }

                q = (b - a) / 2.0;
                p = Math.Pow(q, al + be + 1.0);
                break;
            }
            case 4:
            {
                al = alpha;
                be = beta;
                if (Math.Abs(b - a) <= temp)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SCMM - Fatal error!");
                    Console.WriteLine("  B - A too small!");
                    return null;
                }

                q = (b - a) / 2.0;
                p = Math.Pow(q, al + be + 1.0);
                break;
            }
            case 5 when b <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("SCMM - Fatal error!");
                Console.WriteLine("  B <= 0!");
                return null;
            case 5:
                q = 1.0 / b;
                p = Math.Pow(q, alpha + 1.0);
                break;
            case 6 when b <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("SCMM - Fatal error!");
                Console.WriteLine("  B <= 0!");
                return null;
            case 6:
                q = 1.0 / Math.Sqrt(b);
                p = Math.Pow(q, alpha + 1.0);
                break;
            case 7:
            {
                al = alpha;
                be = 0.0;
                if (Math.Abs(b - a) <= temp)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SCMM - Fatal error!");
                    Console.WriteLine("  B - A too small!");
                    return null;
                }

                q = (b - a) / 2.0;
                p = Math.Pow(q, al + be + 1.0);
                break;
            }
            case 8 when a + b <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("SCMM - Fatal error!");
                Console.WriteLine("  A + B <= 0");
                return null;
            case 8:
                q = a + b;
                p = Math.Pow(q, alpha + beta + 1.0);
                break;
            case 9 when Math.Abs(b - a) <= temp:
                Console.WriteLine("");
                Console.WriteLine("SCMM - Fatal error!");
                Console.WriteLine("  B - A too small!");
                return null;
            case 9:
                q = (b - a) / 2.0;
                p = q * q;
                break;
        }

        //
        //  Compute the moments in W.
        //
        double[] w = WM.wm(m, kind, alpha, beta);

        double tmp = p;

        for (i = 0; i < m; i++)
        {
            w[i] *= tmp;
            tmp *= q;
        }

        return w;
    }
}