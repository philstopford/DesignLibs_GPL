using System;

namespace HermiteIntegrandsTest
{
    public static class Problem03
    {
        public static double p03_exact()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P03_EXACT returns the exact integral for problem 3.
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
            //    Output, double P03_EXACT, the value of the integral.
            //
        {
            double exact;
            double p = 1.0;
            double q = 3.0;
            const double r8_pi = 3.141592653589793;

            exact = r8_pi / (q * Math.Sin(r8_pi * p / q));

            return exact;
        }

        public static void p03_fun(int option, int n, double[] x, ref double[] f )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P03_FUN evaluates the integrand for problem 3.
        //
        //  Discussion:
        //
        //    The exact value is pi / (q*sin(pi*p/q) ), assuming 0 < p < q.
        //
        //    Integral ( -oo < x < +oo ) exp(-px) / ( 1 + exp ( -qx) ) dx
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
            double p = 1.0;
            double q = 3.0;

            for (i = 0; i < n; i++)
            {
                f[i] = Math.Exp(-p * x[i]) / (1.0 + Math.Exp(-q * x[i]));
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

        public static string p03_title()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P03_TITLE returns the title for problem 3.
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
            //    Output, string P03_TITLE, the title of the problem.
            //
        {
            string title;

            title = "exp(-px) / ( 1 + exp(-qx) )";

            return title;
        }
    }
}