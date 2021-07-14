using System;

namespace Burkardt.Laguerre
{
    public static partial class Integrands
    {
        public static double p02_alpha()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P02_ALPHA returns ALPHA for problem 2.
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
            //    Output, double P02_ALPHA, the value of ALPHA.
            //
        {
            double alpha;

            alpha = 2.0;

            return alpha;
        }

        public static double p02_exact()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P02_EXACT returns the exact integral for problem 2.
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
            //    Output, double P02_EXACT, the value of the integral.
            //
        {
            double exact;

            exact = 0.32510848278991335198;

            return exact;
        }

        public static double[] p02_fun(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P02_FUN evaluates the integrand for problem 2.
            //
            //  Discussion:
            //
            //    D&R gives "exact" value as 0.32510855.
            //    Mathematica returns        0.32510848278991335198...
            //    D&R gives Laguerre(16) as  0.19142399...
            //
            //  Integral:
            //
            //    exp ( -2 ) Integral ( 2 <= x < +oo ) / ( x * log(x)^(3/2) ) dx
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
            //    Output, double P02_FUN[N], the function values.
            //
        {
            double[] f;
            int i;

            f = new double[n];

            for (i = 0; i < n; i++)
            {
                f[i] = Math.Exp(-2.0) / (x[i] * Math.Pow(Math.Log(x[i]), 1.5));
            }

            return f;
        }

        public static string p02_title()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P02_TITLE returns the title for problem 2.
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
            //    Output, string P02_TITLE, the title of the problem.
            //
        {
            string title;

            title = "1 / ( x * log ( x )^(3/2) )";

            return title;
        }

    }
}