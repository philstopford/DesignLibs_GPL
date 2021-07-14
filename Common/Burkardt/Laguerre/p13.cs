using System;

namespace Burkardt.Laguerre
{
    public static partial class Integrands
    {
        public static double p13_alpha()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P13_ALPHA returns ALPHA for problem 13.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double P13_ALPHA, the value of ALPHA.
            //
        {
            double alpha;

            alpha = 0.0;

            return alpha;
        }

        public static double p13_exact()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P13_EXACT returns the estimated integral for problem 13.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 July 2007
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

            exact = r8_pi / 2.0;

            return exact;
        }

        public static double[] p13_fun(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P13_FUN evaluates the integrand for problem 13.
            //
            //  Discussion:
            //
            //    S&S gives exact value as pi/2 = 1.5707963267948966192...
            //    S&S gives Laguerre(16) as       1.4399523793...
            //    S&S gives EXP_TRANSFORM(16) as  1.3045186595...
            //
            //  Integral:
            //
            //    Integral ( 0 <= x < +oo ) sin ( x ) / x dx
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud, Don Secrest,
            //    Gaussian Quadrature Formulas,
            //    Prentice Hall, 1966,
            //    LC: QA299.4G3S7.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double X[N], the evaluation points.
            //
            //    Output, double P13_FUN[N], the function values.
            //
        {
            double[] f;
            int i;

            f = new double[n];

            for (i = 0; i < n; i++)
            {
                if (x[i] == 0.0)
                {
                    f[i] = 1.0;
                }
                else
                {
                    f[i] = Math.Sin(x[i]) / x[i];
                }
            }

            return f;
        }

        public static string p13_title()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P13_TITLE returns the title for problem 13.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, string P13_TITLE, the title of the problem.
            //
        {
            string title;

            title = "sin(x) / x";

            return title;
        }

    }
}