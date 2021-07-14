using System;

namespace Burkardt.Laguerre
{
    public static partial class Integrands
    {
        public static double p06_alpha()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P06_ALPHA returns ALPHA for problem 6.
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
            //    Output, double P06_ALPHA, the value of ALPHA.
            //
        {
            double alpha;

            alpha = 2.0;

            return alpha;
        }

        public static double p06_exact()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P06_EXACT returns the exact integral for problem 6.
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
            //    Output, double P06_EXACT, the estimated value of the integral.
            //
        {
            double exact;

            exact = 0.00056103711148387120640;

            return exact;
        }

        public static double[] p06_fun(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P06_FUN evaluates the integrand for problem 6.
            //
            //  Discussion:
            //
            //    D&R gives "exact" value as 0.0005610371...
            //    Mathematica returns        0.00056103711148387120640...
            //    D&R gives Laguerre(16) as  0.00056100775...
            //
            //  Integral:
            //
            //    exp ( -2 ) Integral ( 2 <= x < +oo ) exp ( -x^2 ) dx
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 July 2007
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
            //    Output, double P06_FUN[N], the function values.
            //
        {
            double exponent_min = -80.0;
            double[] f;
            int i;

            f = new double[n];

            for (i = 0; i < n; i++)
            {
                if (-x[i] * x[i] < exponent_min)
                {
                    f[i] = 0.0;
                }
                else
                {
                    f[i] = Math.Exp(-2.0) * Math.Exp(-x[i] * x[i]);
                }
            }

            return f;
        }

        public static string p06_title()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P06_TITLE returns the title for problem 6.
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
            //    Output, string P06_TITLE, the title of the problem.
            //
        {
            string title;

            title = "Complementary error function";

            return title;
        }

    }
}