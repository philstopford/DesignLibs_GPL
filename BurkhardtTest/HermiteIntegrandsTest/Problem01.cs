﻿using System;

namespace HermiteIntegrandsTest
{
    public static class Problem01
    {
        public static double p01_exact()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P01_EXACT returns the exact integral for problem 1.
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
            //    Output, double P01_EXACT, the value of the integral.
            //
        {
            double exact;
            double omega = 1.0;
            const double r8_pi = 3.141592653589793;

            exact = Math.Sqrt(r8_pi) * Math.Exp(-omega * omega);

            return exact;
        }

        public static void p01_fun(int option, int n, double[] x, ref double[] f )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P01_FUN evaluates the integrand for problem 1.
        //
        //  Discussion:
        //
        //    Squire gives exact value as sqrt(pi) * exp(-w*w).
        //
        //    Integral ( -oo < x < +oo ) exp(-x*x) cos(2*w*x) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Squire,
        //    Comparison of Gauss-Hermite and Midpoint Quadrature with Application
        //    to the Voigt Function,
        //    in Numerical Integration: 
        //    Recent Developments, Software and Applications,
        //    edited by Patrick Keast, Graeme Fairweather,
        //    Reidel, 1987, pages 337-340,
        //    ISBN: 9027725144,
        //    LC: QA299.3.N38.
        //
        //  Parameters:
        //
        //    Input, int OPTION:
        //    0, integrand is f(x).
        //    1, integrand is exp(-x*x) * f(x);
        //    2, integrand is exp(-x*x/2) * f(x);
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[N], the evaluation points.
        //
        //    Output, double F[N], the function values.
        //
        {
            int i;
            double omega = 1.0;

            for (i = 0; i < n; i++)
            {
                f[i] = Math.Cos(2.0 * omega * x[i]);
            }

            if (option == 0)
            {
                for (i = 0; i < n; i++)
                {
                    f[i] = f[i] * Math.Exp(-x[i] * x[i]);
                }
            }
            else if (option == 1)
            {
            }
            else if (option == 2)
            {
                for (i = 0; i < n; i++)
                {
                    f[i] = f[i] * Math.Exp(-x[i] * x[i] / 2.0);
                }
            }

            return;
        }

        public static string p01_title()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P01_TITLE returns the title for problem 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    26 May 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, string P01_TITLE, the title of the problem.
            //
        {
            string title;

            title = "exp(-x*x) * cos(2*omega*x)";

            return title;
        }
    }
}