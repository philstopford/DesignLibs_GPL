using System;

namespace HermiteIntegrandsTest
{
    public static class Problem07
    {
        public static double p07_exact()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P07_EXACT returns the exact integral for problem 7.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 February 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double P07_EXACT, the value of the integral.
            //
        {
            double e_sqrt_sqrt = 1.2840254166877414841;
            double exact;
            const double r8_pi = 3.141592653589793;

            exact = 0.25 * Math.Sqrt(r8_pi) / e_sqrt_sqrt;

            return exact;
        }

        public static void p07_fun(int option, int n, double[] x, ref double[] f )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P07_FUN evaluates the integrand for problem 7.
        //
        //  Discussion:
        //
        //    The exact value is (1/4) sqrt(pi) / sqrt(sqrt(e)).
        //
        //    Integral ( -oo < x < +oo ) x^2 cos(x) e^(-x^2) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 February 2010
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

            if (option == 0)
            {
                for (i = 0; i < n; i++)
                {
                    f[i] = Math.Pow(x[i], 2) * Math.Cos(x[i]) * Math.Exp(-x[i] * x[i]);
                }
            }
            else if (option == 1)
            {
                for (i = 0; i < n; i++)
                {
                    f[i] = Math.Pow(x[i], 2) * Math.Cos(x[i]);
                }
            }
            else if (option == 2)
            {
                for (i = 0; i < n; i++)
                {
                    f[i] = Math.Pow(x[i], 2) * Math.Cos(x[i]) * Math.Exp(-x[i] * x[i] / 2.0);
                }
            }

            return;
        }

        public static string p07_title()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P07_TITLE returns the title for problem 7.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 February 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, string P07_TITLE, the title of the problem.
            //
        {
            string title;

            title = "x^2 cos ( x ) exp(-x*x)";

            return title;
        }
    }
}