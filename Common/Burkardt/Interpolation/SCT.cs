using System;
using Burkardt.Types;

namespace Burkardt.Interpolation
{
    public static class SCT
    {
        public static double[] sct(int nt, double[] t, int kind, double a, double b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SCT rescales distinct knots to an interval [A,B].
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
        //    Input, double T[NT], the original knots.
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
        //    Input, double A, B, the interval endpoints for which the
        //    knots ST should be scaled.
        //
        //    Output, double SCT[NT], the scaled knots.
        //
        {
            double bma;
            int i;
            double shft = 0;
            double slp = 0;
            double[] st;
            double tmp;

            if (kind < 1 || 9 < kind)
            {
                Console.WriteLine("");
                Console.WriteLine("SCT - Fatal error!");
                Console.WriteLine("  KIND falls outside range of 1 to 8.");
                return null;
            }

            if (kind == 1 || kind == 2 || kind == 3 || kind == 4 || kind == 7 || kind == 9)
            {
                tmp = typeMethods.r8_epsilon();
                bma = b - a;

                if (bma <= tmp)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SCT - Fatal error!");
                    Console.WriteLine("  B - A too small.");
                    return null;
                }

                slp = 2.0 / bma;
                shft = -(a + b) / bma;
            }
            else if (kind == 5)
            {
                if (b < 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SCT - Fatal error!");
                    Console.WriteLine("  B < 0.");
                    return null;
                }

                slp = b;
                shft = -a * b;
            }
            else if (kind == 6)
            {
                if (b < 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SCT - Fatal error!");
                    Console.WriteLine("  B < 0.");
                    return null;
                }

                slp = Math.Sqrt(b);
                shft = -a * slp;
            }
            else if (kind == 8)
            {
                slp = 1.0 / (a + b);

                if (slp <= 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SCT - Fatal error.");
                    Console.WriteLine("  1 / ( A + B ) <= 0.");
                    return null;
                }

                shft = -a * slp;
            }

            st = new double[nt];

            for (i = 0; i < nt; i++)
            {
                st[i] = shft + slp * t[i];
            }

            return st;
        }
    }
}