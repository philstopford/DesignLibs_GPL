using System;

namespace Burkardt.Laguerre;

public static partial class Integrands
{
    public static double p09_alpha()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P09_ALPHA returns ALPHA for problem 9.
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
        //    Output, double P09_ALPHA, the value of ALPHA.
        //
    {
        double alpha;

        alpha = 0.0;

        return alpha;
    }

    public static double p09_exact()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P09_EXACT returns the estimated integral for problem 9.
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
        //    Output, double P09_EXACT, the estimated value of the integral.
        //
    {
        double exact;

        exact = 24.0;

        return exact;
    }

    public static double[] p09_fun(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P09_FUN evaluates the integrand for problem 9.
        //
        //  Discussion:
        //
        //    The integral is the definition of the Gamma function for
        //    Z = 5, with exact value (Z-1)! = 24.
        //
        //  Integral:
        //
        //    Integral ( 0 <= x < +oo ) x^4 exp ( -x ) dx
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
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[N], the evaluation points.
        //
        //    Output, double P09_FUN[N], the function values.
        //
    {
        double exponent_min = -80.0;
        double[] f;
        int i;

        f = new double[n];

        for (i = 0; i < n; i++)
        {
            if (-x[i] < exponent_min)
            {
                f[i] = 0.0;
            }
            else
            {
                f[i] = Math.Pow(x[i], 4) * Math.Exp(-x[i]);
            }
        }

        return f;
    }

    public static string p09_title()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P09_TITLE returns the title for problem 9.
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
        //    Output, string P09_TITLE, the title of the problem.
        //
    {
        string title;

        title = "Gamma(Z=5) function";

        return title;
    }
}