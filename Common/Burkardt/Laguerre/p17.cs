﻿using System;

namespace Burkardt.Laguerre
{
    public static partial class Integrands
    {
        public static double p17_alpha()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P17_ALPHA returns ALPHA for problem 17.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 December 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double P17_ALPHA, the value of ALPHA.
            //
        {
            double alpha;

            alpha = 0.0;

            return alpha;
        }

        public static double p17_exact()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P17_EXACT returns the exact integral for problem 17.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 December 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double P17_EXACT, the value of the integral.
            //
        {
            const double beta = 2.0;
            double exact;
            const double r8_pi = 3.1415926535897932385;

            exact = Math.Sqrt(r8_pi) * Math.Cos(0.5 * Math.Atan(Math.Pow(2.0, beta)))
                    / Math.Sqrt(Math.Sqrt(1.0 + Math.Pow(0.25, beta)));

            return exact;
        }

        public static double[] p17_fun(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P17_FUN evaluates the integrand for problem 17.
            //
            //  Integral:
            //
            //    Integral ( 0 <= x < +oo) exp ( - x / 2^beta ) * cos ( x ) / Math.Sqrt ( x ) dx
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 December 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Robert Piessens, Elise de Doncker-Kapenga, 
            //    Christian Ueberhuber, David Kahaner,
            //    QUADPACK: A Subroutine Package for Automatic Integration,
            //    Springer, 1983, page 84.
            //
            //  Parameters:
            //
            //    Input, int N, the number of evaluation points.
            //
            //    Input, double X[N], the evaluation points.
            //
            //    Output, double P17_FUN[N], the integrand values.
            //
        {
            double beta = 2.0;
            double[] fx;
            int i;

            fx = new double[n];

            for (i = 0; i < n; i++)
            {
                if (x[i] == 0.0)
                {
                    fx[i] = 0.0;
                }
                else
                {
                    fx[i] = Math.Exp(-x[i] / Math.Pow(2.0, beta)) * Math.Cos(x[i])
                            / Math.Sqrt(x[i]);
                }
            }

            return fx;
        }

        public static string p17_title()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P17_TITLE returns the title for problem 17.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 December 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, string P17_TITLE, the title of the problem.
            //
        {
            string title;

            title = "exp ( - x / 2^beta ) * cos ( x ) / Math.Sqrt ( x )";

            return title;
        }
    }
}