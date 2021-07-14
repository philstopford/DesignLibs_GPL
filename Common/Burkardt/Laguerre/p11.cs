using System;

namespace Burkardt.Laguerre
{
    public static partial class Integrands
    {
        public static double p11_alpha()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P11_ALPHA returns ALPHA for problem 11.
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
            //    Output, double P11_ALPHA, the value of ALPHA.
            //
        {
            double alpha;

            alpha = 0.0;

            return alpha;
        }

        public static double p11_exact()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P11_EXACT returns the estimated integral for problem 11.
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

            exact = r8_pi;

            return exact;
        }

        public static double[] p11_fun(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P11_FUN evaluates the integrand for problem 11.
            //
            //  Discussion:
            //
            //    S&S gives exact value as pi =  3.1415926535897932385...
            //    S&S gives Laguerre(16) as      2.6652685196...
            //    S&S gives EXP_TRANSFORM(16) as 2.3629036166... 
            //
            //  Integral:
            //
            //    Integral ( 0 <= x < +oo ) 1/((1+x)*sqrt(x)) dx
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
            //    Output, double P11_FUN[N], the function values.
            //
        {
            double[] f;
            int i;

            f = new double[n];

            for (i = 0; i < n; i++)
            {
                if (x[i] == 0.0)
                {
                    f[i] = 0.0;
                }
                else
                {
                    f[i] = 1.0 / ((1.0 + x[i]) * Math.Sqrt(x[i]));
                }
            }

            return f;
        }

        public static string p11_title()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P11_TITLE returns the title for problem 11.
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
            //    Output, string P11_TITLE, the title of the problem.
            //
        {
            string title;

            title = "1 / ( (1+x) * sqrt(x) )";

            return title;
        }

    }
}