using System;
using Burkardt.Types;

namespace Burkardt.Laguerre
{
    public static partial class Integrands
    {
        public static double p16_alpha()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P16_ALPHA returns ALPHA for problem 16.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 August 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double P16_ALPHA, the value of ALPHA.
            //
        {
            double alpha;

            alpha = 0.0;

            return alpha;
        }

        public static double p16_exact()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P16_EXACT returns the estimated integral for problem 16.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 August 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double EXACT, the estimated value of the integral.
            //
        {
            double exact;

            exact = 1.0;

            return exact;
        }

        public static double[] p16_fun(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P16_FUN evaluates the integrand for problem 16.
            //
            //  Integral:
            //
            //    Integral ( 0 <= x < +oo ) cos ( pi * x / 2 ) / sqrt ( x ) dx
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 August 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Robert Piessens, Elise deDoncker-Kapenga, 
            //    Christian Ueberhuber, David Kahaner,
            //    QUADPACK: A Subroutine Package for Automatic Integration,
            //    Springer, 1983,
            //    ISBN: 3540125531,
            //    LC: QA299.3.Q36.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double X[N], the evaluation points.
            //
            //    Output, double P16_FUN[N], the function values.
            //
        {
            double[] f;
            int i;
            const double r8_pi = 3.1415926535897932385;

            f = new double[n];

            for (i = 0; i < n; i++)
            {
                if (x[i] == 0.0)
                {
                    f[i] = typeMethods.r8_huge();
                }
                else
                {
                    f[i] = Math.Cos(r8_pi * x[i] / 2.0) / Math.Sqrt(x[i]);
                }
            }

            return f;
        }

        public static string p16_title()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P16_TITLE returns the title for problem 16.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 August 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, string P16_TITLE, the title of the problem.
            //
        {
            string title;

            title = "cos ( pi x / 2 ) / sqrt ( x )";

            return title;
        }

    }
}