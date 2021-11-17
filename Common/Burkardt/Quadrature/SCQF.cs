using System;
using Burkardt.Types;
using Burkardt.Weight;

namespace Burkardt.Quadrature;

public static class SCQF
{
    public static void scqf(int nt, double[] t, int[] mlt, double[] wts, int nwts, int[] ndx,
            ref double[] swts, ref double[] st, int kind, double alpha, double beta, double a,
            double b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SCQF scales a quadrature formula to a nonstandard interval.
        //
        //  Discussion:
        //
        //    The arrays WTS and SWTS may coincide.
        //
        //    The arrays T and ST may coincide.
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
        //    Input, int MLT[NT], the multiplicity of the knots.
        //
        //    Input, double WTS[NWTS], the weights.
        //
        //    Input, int NWTS, the number of weights.
        //
        //    Input, int NDX[NT], used to index the array WTS.  
        //    For more details see the comments in CAWIQ.
        //
        //    Output, double SWTS[NWTS], the scaled weights.
        //
        //    Output, double ST[NT], the scaled knots.
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
    {
        double al = 0;
        double be = 0;
        double shft = 0;
        double slp = 0;

        double temp = typeMethods.r8_epsilon();

        PARCHK.parchk(kind, 1, alpha, beta);

        switch (kind)
        {
            case 1:
            {
                al = 0.0;
                be = 0.0;
                if (Math.Abs(b - a) <= temp)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SCQF - Fatal error!");
                    Console.WriteLine("  |B - A| too small.");
                    return;
                }

                shft = (a + b) / 2.0;
                slp = (b - a) / 2.0;
                break;
            }
            case 2:
            {
                al = -0.5;
                be = -0.5;
                if (Math.Abs(b - a) <= temp)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SCQF - Fatal error!");
                    Console.WriteLine("  |B - A| too small.");
                    return;
                }

                shft = (a + b) / 2.0;
                slp = (b - a) / 2.0;
                break;
            }
            case 3:
            {
                al = alpha;
                be = alpha;
                if (Math.Abs(b - a) <= temp)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SCQF - Fatal error!");
                    Console.WriteLine("  |B - A| too small.");
                    return;
                }

                shft = (a + b) / 2.0;
                slp = (b - a) / 2.0;
                break;
            }
            case 4:
            {
                al = alpha;
                be = beta;

                if (Math.Abs(b - a) <= temp)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SCQF - Fatal error!");
                    Console.WriteLine("  |B - A| too small.");
                    return;
                }

                shft = (a + b) / 2.0;
                slp = (b - a) / 2.0;
                break;
            }
            case 5 when b <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("SCQF - Fatal error!");
                Console.WriteLine("  B <= 0");
                return;
            case 5:
                shft = a;
                slp = 1.0 / b;
                al = alpha;
                be = 0.0;
                break;
            case 6 when b <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("SCQF - Fatal error!");
                Console.WriteLine("  B <= 0.");
                return;
            case 6:
                shft = a;
                slp = 1.0 / Math.Sqrt(b);
                al = alpha;
                be = 0.0;
                break;
            case 7:
            {
                al = alpha;
                be = 0.0;
                if (Math.Abs(b - a) <= temp)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SCQF - Fatal error!");
                    Console.WriteLine("  |B - A| too small.");
                    return;
                }

                shft = (a + b) / 2.0;
                slp = (b - a) / 2.0;
                break;
            }
            case 8 when a + b <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("SCQF - Fatal error!");
                Console.WriteLine("  A + B <= 0.");
                return;
            case 8:
                shft = a;
                slp = a + b;
                al = alpha;
                be = beta;
                break;
            case 9:
            {
                al = 0.5;
                be = 0.5;
                if (Math.Abs(b - a) <= temp)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SCQF - Fatal error!");
                    Console.WriteLine("  |B - A| too small.");
                    return;
                }

                shft = (a + b) / 2.0;
                slp = (b - a) / 2.0;
                break;
            }
        }

        double p = Math.Pow(slp, al + be + 1.0);

        for (int k = 0; k < nt; k++)
        {
            st[k] = shft + slp * t[k];
            int l = Math.Abs(ndx[k]);

            if (l != 0)
            {
                double tmp = p;
                for (int i = l - 1; i <= l - 1 + mlt[k] - 1; i++)
                {
                    swts[i] = wts[i] * tmp;
                    tmp *= slp;
                }
            }
        }
    }
}