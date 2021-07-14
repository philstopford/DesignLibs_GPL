using System;
using Burkardt.Types;

namespace Burkardt.Laguerre
{
    public static partial class Integrands
    {
        public static double p15_alpha()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P15_ALPHA returns ALPHA for problem 15.
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
            //    Output, double P15_ALPHA, the value of ALPHA.
            //
        {
            double alpha;

            alpha = 0.0;

            return alpha;
        }

        public static double p15_exact()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P15_EXACT returns the estimated integral for problem 15.
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
            const double r8_pi = 3.1415926535897932385;

            exact = -r8_pi * Math.Log(10.0) / 20.0;

            return exact;
        }

        public static double[] p15_fun(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P15_FUN evaluates the integrand for problem 15.
            //
            //  Integral:
            //
            //    Integral ( 0 <= x < +oo ) log(x) / (1+100*x*x) dx
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
            //    Output, double P15_FUN[N], the function values.
            //
        {
            double[] f;
            int i;

            f = new double[n];

            for (i = 0; i < n; i++)
            {
                if (x[i] == 0.0)
                {
                    f[i] = -typeMethods.r8_huge();
                }
                else
                {
                    f[i] = Math.Log(x[i]) / (1.0 + 100.0 * x[i] * x[i]);
                }
            }

            return f;
        }

        public static string p15_title()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P15_TITLE returns the title for problem 15.
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
            //    Output, string P15_TITLE, the title of the problem.
            //
        {
            string title;

            title = "log(x) / ( 1 + 100 x^2 )";

            return title;
        }

    }
}