using System;

namespace HermiteIntegrandsTest;

public static class Problem04
{
    public static double p04_exact()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P04_EXACT returns the estimated integral for problem 4.
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
        //    Output, double P04_EXACT, the estimated value of the integral.
        //
    {
        double exact;
        const double r8_pi = 3.141592653589793;

        exact = Math.Sqrt(r8_pi / 2.0);

        return exact;
    }

    public static void p04_fun(int option, int n, double[] x, ref double[] f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P04_FUN evaluates the integrand for problem 4.
        //
        //  Discussion:
        //
        //    The exact value is sqrt ( pi / 2 )
        //
        //    Integral ( -oo < x < +oo ) sin ( x^2 ) dx
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
            f[i] = Math.Sin(x[i] * x[i]);
        }

        switch (option)
        {
            case 0:
                break;
            case 1:
            {
                for (i = 0; i < n; i++)
                {
                    f[i] *= Math.Exp(+x[i] * x[i]);
                }

                break;
            }
            case 2:
            {
                for (i = 0; i < n; i++)
                {
                    f[i] *= Math.Exp(+x[i] * x[i] / 2.0);
                }

                break;
            }
        }
    }

    public static string p04_title()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P04_TITLE returns the title for problem 4.
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
        //    Output, string P04_TITLE, the title of the problem.
        //
    {
        string title;

        title = "sin(x^2)";

        return title;
    }
}