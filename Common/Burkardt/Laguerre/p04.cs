﻿using System;

namespace Burkardt.Laguerre
{
    public static partial class Integrands
    {
        public static double p04_alpha()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P04_ALPHA returns ALPHA for problem 4.
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
            //    Output, double P04_ALPHA, the value of ALPHA.
            //
        {
            double alpha;

            alpha = 2.0;

            return alpha;
        }

        public static double p04_exact()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P04_EXACT returns the estimated integral for problem 4.
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
            //    Output, double P04_EXACT, the estimated value of the integral.
            //
        {
            double exact;

            exact = -0.0046848541335080643181;

            return exact;
        }

        public static double[] p04_fun(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P04_FUN evaluates the integrand for problem 4.
            //
            //  Discussion:
            //
            //    D&R gives "exact" value as -0.0046984...
            //    Mathematica returns        -0.0046848541335080643181...
            //    D&R gives Laguerre(16) as  -0.039258696...
            //
            //  Integral:
            //
            //    exp ( -2 ) Integral ( 2 <= x < +oo ) sin ( x ) / x dx
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
            //    Output, double P04_FUN[N], the function values.
            //
        {
            double[] f;
            int i;

            f = new double[n];

            for (i = 0; i < n; i++)
            {
                if (x[i] == 0.0)
                {
                    f[i] = Math.Exp(-2.0);
                }
                else
                {
                    f[i] = Math.Exp(-2.0) * Math.Sin(x[i]) / x[i];
                }
            }

            return f;
        }

        public static string p04_title()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P04_TITLE returns the title for problem 4.
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
            //    Output, string P04_TITLE, the title of the problem.
            //
        {
            string title;

            title = "Sine integral";

            return title;
        }

    }
}