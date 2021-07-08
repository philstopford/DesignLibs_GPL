using System;

namespace HermiteIntegrandsTest
{
    public static class Problem02
    {
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
            //    31 July 2007
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
            const double r8_pi = 3.141592653589793;

            exact = Math.Sqrt(r8_pi);

            return exact;
        }

        public static void p02_fun(int option, int n, double[] x, ref double[] f )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P02_FUN evaluates the integrand for problem 2.
        //
        //  Discussion:
        //
        //    The exact value is sqrt(pi).
        //
        //    Integral ( -oo < x < +oo ) exp(-x*x) dx
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
                f[i] = 1.0;
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
            //    26 May 2009
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

            title = "exp(-x*x)";

            return title;
        }
    }
}