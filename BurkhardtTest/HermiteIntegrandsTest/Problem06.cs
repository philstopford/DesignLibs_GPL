using System;
using Burkardt.Types;

namespace HermiteIntegrandsTest
{
    public static class Problem06
    {
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
            //    26 May 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double P06_EXACT, the value of the integral.
            //
        {
            double exact;
            int m = 0;
            const double r8_pi = 3.141592653589793;

            p06_param('G', 'M', ref m);

            if (m <= -1)
            {
                exact = -typeMethods.r8_huge();
            }
            else if ((m % 2) == 1)
            {
                exact = 0.0;
            }
            else
            {
                exact = (double) (typeMethods.i4_factorial2(m - 1)) * Math.Sqrt(r8_pi)
                        / Math.Pow(2.0, m / 2);
            }

            return exact;
        }

        public static void p06_fun(int option, int n, double[] x, ref double[] f )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P06_FUN evaluates the integrand for problem 6.
        //
        //  Discussion:
        //
        //    The exact value is (m-1)!! * sqrt ( pi ) / sqrt ( 2**m ).
        //
        //    Integral ( -oo < x < +oo ) x^m exp (-x*x) dx
        //
        //    The parameter M is set by calling P06_PARAM.
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
            int m = 0;

            p06_param('G', 'M', ref m);

            for (i = 0; i < n; i++)
            {
                f[i] = Math.Pow(x[i], m);
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

            return;
        }

        public static void p06_param(char action, char name, ref int value)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P06_PARAM gets or sets parameters for problem 6.
            //
            //  Discussion:
            //
            //    The parameter is named "M", and it represents the value of the exponent
            //    in the integrand function:
            //
            //    Integral ( -oo < x < +oo ) x^m exp (-x*x) dx
            //
            //    M must be greater than -1.
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
            //    Input, char ACTION, the action.
            //    'S' to set the value,
            //    'G' to get the value.
            //
            //    Input, char NAME, the parameter name.
            //    'M', the exponent.
            //
            //    Input/output, int *VALUE, the parameter value.
            //    If ACTION = 'S', then VALUE is an input quantity, and M is set to VALUE.
            //    If ACTION = 'G', then VALUE is an output quantity, and VALUE is set to M.
            //
        {
            int m = 0;

            if (action == 'S' || action == 's')
            {
                if (value <= -1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("P06_PARAM - Fatal error!");
                    Console.WriteLine("  Parameter M must be greater than -1.");
                    return;
                }

                m = value;
            }
            else if (action == 'G' || action == 'g')
            {
                value = m;
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("P06_PARAM - Fatal error!");
                Console.WriteLine("  Unrecognized value of ACTION = \"" + action + "\".");
            }
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
            //    26 May 2009
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

            title = "x^m exp(-x*x)";

            return title;
        }
    }
}