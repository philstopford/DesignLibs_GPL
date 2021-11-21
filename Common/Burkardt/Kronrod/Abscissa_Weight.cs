﻿using System;
using Burkardt.Types;

namespace Burkardt.Kronrod;

public static class Abscissa_Weight
{
    public static void abwe1(int n, int m, double eps, double coef2, bool even, double[] b,
            ref double x, ref double w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ABWE1 calculates a Kronrod abscissa and weight.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 August 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Robert Piessens, Maria Branders.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Robert Piessens, Maria Branders,
        //    A Note on the Optimal Addition of Abscissas to Quadrature Formulas
        //    of Gauss and Lobatto,
        //    Mathematics of Computation,
        //    Volume 28, Number 125, January 1974, pages 135-139.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the Gauss rule.
        //
        //    Input, int M, the value of ( N + 1 ) / 2.
        //
        //    Input, double EPS, the requested absolute accuracy of the
        //    abscissas.
        //
        //    Input, double COEF2, a value needed to compute weights.
        //
        //    Input, bool EVEN, is TRUE if N is even.
        //
        //    Input, double B[M+1], the Chebyshev coefficients.
        //
        //    Input/output, double *X; on input, an estimate for
        //    the abscissa, and on output, the computed abscissa.
        //
        //    Output, double *W, the weight.
        //
    {
        double ai;
        double b0 = 0;
        double d0;
        double d1;
        double d2 = 0;
        double delta = 0;
        double fd = 0;
        int iter;
        int k;

        int ka = x switch
        {
            0.0 => 1,
            _ => 0
        };

        //
        //  Iterative process for the computation of a Kronrod abscissa.
        //
        for (iter = 1; iter <= 50; iter++)
        {
            double b1 = 0.0;
            double b2 = b[m];
            double yy = 4.0 * x * x - 2.0;
            d1 = 0.0;

            double dif;
            switch (even)
            {
                case true:
                    ai = m + m + 1;
                    d2 = ai * b[m];
                    dif = 2.0;
                    break;
                default:
                    ai = m + 1;
                    d2 = 0.0;
                    dif = 1.0;
                    break;
            }

            for (k = 1; k <= m; k++)
            {
                ai -= dif;
                int i = m - k + 1;
                b0 = b1;
                b1 = b2;
                d0 = d1;
                d1 = d2;
                b2 = yy * b1 - b0 + b[i - 1];
                switch (even)
                {
                    case false:
                        i += 1;
                        break;
                }

                d2 = yy * d1 - d0 + ai * b[i - 1];
            }

            double f;
            switch (even)
            {
                case true:
                    f = x * (b2 - b1);
                    fd = d2 + d1;
                    break;
                default:
                    f = 0.5 * (b2 - b0);
                    fd = 4.0 * x * d2;
                    break;
            }

            //
            //  Newton correction.
            //
            delta = f / fd;
            x -= delta;

            if (ka == 1)
            {
                break;
            }

            if (Math.Abs(delta) <= eps)
            {
                ka = 1;
            }
        }

        //
        //  Catch non-convergence.
        //
        if (ka != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("ABWE1 - Fatal error!");
            Console.WriteLine("  Iteration limit reached.");
            Console.WriteLine("  EPS is " + eps + "");
            Console.WriteLine("  Last DELTA was " + delta + "");
            return;
        }

        //
        //  Computation of the weight.
        //
        d0 = 1.0;
        d1 = x;
        ai = 0.0;
        for (k = 2; k <= n; k++)
        {
            ai += 1.0;
            d2 = ((ai + ai + 1.0) * x * d1 - ai * d0) / (ai + 1.0);
            d0 = d1;
            d1 = d2;
        }

        w = coef2 / (fd * d2);
    }

    public static void abwe2(int n, int m, double eps, double coef2, bool even, double[] b,
            ref double x, ref double w1, ref double w2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ABWE2 calculates a Gaussian abscissa and two weights.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 April 2013
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Robert Piessens, Maria Branders.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Robert Piessens, Maria Branders,
        //    A Note on the Optimal Addition of Abscissas to Quadrature Formulas
        //    of Gauss and Lobatto,
        //    Mathematics of Computation,
        //    Volume 28, Number 125, January 1974, pages 135-139.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the Gauss rule.
        //
        //    Input, int M, the value of ( N + 1 ) / 2.
        //
        //    Input, double EPS, the requested absolute accuracy of the
        //    abscissas.
        //
        //    Input, double COEF2, a value needed to compute weights.
        //
        //    Input, bool EVEN, is TRUE if N is even.
        //
        //    Input, double B[M+1], the Chebyshev coefficients.
        //
        //    Input/output, double *X; on input, an estimate for
        //    the abscissa, and on output, the computed abscissa.
        //
        //    Output, double *W1, the Gauss-Kronrod weight.
        //
        //    Output, double *W2, the Gauss weight.
        //
    {
        double delta = 0;
        int iter;
        int k;
        double p0 = 0;
        double p1;
        double p2 = 0;
        double pd2 = 0;

        int ka = x switch
        {
            0.0 => 1,
            _ => 0
        };

        //
        //  Iterative process for the computation of a Gaussian abscissa.
        //
        for (iter = 1; iter <= 50; iter++)
        {
            p0 = 1.0;
            p1 = x;
            double pd0 = 0.0;
            double pd1 = 1.0;
            switch (n)
            {
                //
                //  When N is 1, we need to initialize P2 and PD2 to avoid problems with DELTA.
                //
                case <= 1 when typeMethods.r8_epsilon() < Math.Abs(x):
                    p2 = (3.0 * x * x - 1.0) / 2.0;
                    pd2 = 3.0 * x;
                    break;
                case <= 1:
                    p2 = 3.0 * x;
                    pd2 = 3.0;
                    break;
            }

            double ai = 0.0;
            for (k = 2; k <= n; k++)
            {
                ai += 1.0;
                p2 = ((ai + ai + 1.0) * x * p1 - ai * p0) / (ai + 1.0);
                pd2 = ((ai + ai + 1.0) * (p1 + x * pd1) - ai * pd0)
                      / (ai + 1.0);
                p0 = p1;
                p1 = p2;
                pd0 = pd1;
                pd1 = pd2;
            }

            //
            //  Newton correction.
            //
            delta = p2 / pd2;
            x -= delta;

            if (ka == 1)
            {
                break;
            }

            if (Math.Abs(delta) <= eps)
            {
                ka = 1;
            }
        }

        //
        //  Catch non-convergence.
        //
        if (ka != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("ABWE2 - Fatal error!");
            Console.WriteLine("  Iteration limit reached.");
            Console.WriteLine("  EPS is " + eps + "");
            Console.WriteLine("  Last DELTA was " + delta + "");
            return;
        }

        //
        //  Computation of the weight.
        //
        double an = n;

        w2 = 2.0 / (an * pd2 * p0);

        p1 = 0.0;
        p2 = b[m];
        double yy = 4.0 * x * x - 2.0;
        for (k = 1; k <= m; k++)
        {
            int i = m - k + 1;
            p0 = p1;
            p1 = p2;
            p2 = yy * p1 - p0 + b[i - 1];
        }

        w1 = even switch
        {
            true => w2 + coef2 / (pd2 * x * (p2 - p1)),
            _ => w2 + 2.0 * coef2 / (pd2 * (p2 - p0))
        };
    }
}