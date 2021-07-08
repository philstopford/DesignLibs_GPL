using System;

namespace HermiteIntegrandsTest
{
    public static class Problem05
    {
        public static double p05_exact()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P05_EXACT returns the estimated integral for problem 5.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 October 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double P05_EXACT, the estimated value of the integral.
            //
        {
            double exact;
            const double r8_pi = 3.141592653589793;

            exact = r8_pi / 3.0;

            return exact;
        }

        public static void p05_fun(int option, int n, double[] x, ref double[] f )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P05_FUN evaluates the integrand for problem 5.
        //
        //  Discussion:
        //
        //    The exact value is pi / 3.
        //
        //    Integral ( -oo < x < +oo ) dx / ( (1+x^2) sqrt(4+3x^2) )
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
                f[i] = 1.0 / ((1.0 + x[i] * x[i]) * Math.Sqrt(4.0 + 3.0 * x[i] * x[i]));
            }

            if (option == 0)
            {
            }
            else if (option == 1)
            {
                for (i = 0; i < n; i++)
                {
                    f[i] = f[i] * Math.Exp(+x[i] * x[i]);
                }
            }
            else if (option == 2)
            {
                for (i = 0; i < n; i++)
                {
                    f[i] = f[i] * Math.Exp(+x[i] * x[i] / 2.0);
                }
            }

            return;
        }

        public static string p05_title()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P05_TITLE returns the title for problem 5.
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
            //    Output, string P05_TITLE, the title of the problem.
            //
        {
            string title;

            title = "1/( (1+x^2) sqrt(4+3x^2) )";

            return title;
        }
    }
}