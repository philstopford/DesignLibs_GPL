using System;

namespace HermiteIntegrandsTest
{
    public static class Problem08
    {
        public static double p08_exact()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P08_EXACT returns the exact integral for problem 8.
            //
            //  Discussion:
            //
            //    The 20 digit value of the answer was computed by Mathematica.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 July 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double P08_EXACT, the value of the integral.
            //
        {
            double exact;

            exact = 3.0088235661136433510;

            return exact;
        }

        public static void p08_fun(int option, int n, double[] x, ref double[] f )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P08_FUN evaluates the integrand for problem 8.
        //
        //  Discussion:
        //
        //    The exact value is sqrt ( 2 pi ) * HypergeometricU ( -1/2, 0, 1 ).
        //
        //    Integral ( -oo < x < +oo ) sqrt(1+x*x/2) * exp(-x*x/2) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2010
        //
        //  Author:
        //
        //    John Burkardt
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

            for (i = 0; i < n; i++)
            {
                f[i] = Math.Sqrt(1.0 + 0.5 * x[i] * x[i]);
            }

            if (option == 0)
            {
                for (i = 0; i < n; i++)
                {
                    f[i] = f[i] * Math.Exp(-0.5 * x[i] * x[i]);
                }
            }
            else if (option == 1)
            {
                for (i = 0; i < n; i++)
                {
                    f[i] = f[i] * Math.Exp(+0.5 * x[i] * x[i]);
                }
            }
            else if (option == 2)
            {

            }
        }

        public static string p08_title()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P08_TITLE returns the title for problem 8.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 July 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, string P08_TITLE, the title of the problem.
            //
        {
            string title;

            title = "sqrt(1+x*x/2) * exp(-x*x/2)";

            return title;
        }
    }
}