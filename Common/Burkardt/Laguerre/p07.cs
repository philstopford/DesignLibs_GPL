﻿using System;

namespace Burkardt.Laguerre
{
    public static partial class Integrands
    {
        public static double p07_alpha()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P07_ALPHA returns ALPHA for problem 7.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double P07_ALPHA, the value of ALPHA.
            //
        {
            double alpha;

            alpha = 2.0;

            return alpha;
        }

        public static double p07_exact()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P07_EXACT returns the exact integral for problem 7.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double P07_EXACT, the value of the integral.
            //
        {
            double exact;

            exact = 0.16266891;

            return exact;
        }

        public static double[] p07_fun(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P07_FUN evaluates the integrand for problem 7.
            //
            //  Discussion:
            //
            //    D&R gives "exact" value as 0.16266891...
            //    Mathematica does not return a value.
            //    D&R gives Laguerre(16) as  0.097083064...
            //
            //  Integral:
            //
            //    exp ( -2 ) Integral ( 2 <= x < +oo ) sin ( x - 1 ) 
            //      / sqrt ( x * ( x - 2 ) ) dx
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Philip Davis, Philip Rabinowitz,
            //    Methods of Numerical Integration,
            //    Second Edition,
            //    Dover, 2007,
            //    ISBN: 0486453391,
            //    LC: QA299.3.D28.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double X[N], the evaluation points.
            //
            //    Output, double P07_FUN[N], the function values.
            //
        {
            double[] f;
            int i;

            f = new double[n];

            for (i = 0; i < n; i++)
            {
                if (x[i] == 2.0)
                {
                    f[i] = 0.0;
                }
                else
                {
                    f[i] = Math.Exp(-2.0)
                        * Math.Sin(x[i] - 1.0) / Math.Sqrt(x[i] * (x[i] - 2.0));
                }
            }

            return f;
        }

        public static string p07_title()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P07_TITLE returns the title for problem 7.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, string P07_TITLE, the title of the problem.
            //
        {
            string title;

            title = "Bessel function";

            return title;
        }

    }
}