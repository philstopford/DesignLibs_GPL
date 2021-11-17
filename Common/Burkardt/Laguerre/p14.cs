using System;

namespace Burkardt.Laguerre;

public static partial class Integrands
{
    public static double p14_alpha()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P14_ALPHA returns ALPHA for problem 14.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double P14_ALPHA, the value of ALPHA.
        //
    {
        double alpha;

        alpha = 0.0;

        return alpha;
    }

    public static double p14_exact()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P14_EXACT returns the estimated integral for problem 14.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double EXACT, the estimated value of the integral.
        //
    {
        double exact;

        exact = 1.0634618101722400407;

        return exact;
    }

    public static double[] p14_fun(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P14_FUN evaluates the integrand for problem 14.
        //
        //  Discussion:
        //
        //    S&S gives "exact" value as     1.0634618101...
        //    Mathematica returns            1.0634618101722400407...
        //    S&S gives Laguerre(16) as      1.0634713425...
        //    S&S gives EXP_TRANSFORM(16) as 1.0634618101...
        //
        //  Integral:
        //
        //    Integral ( 0 <= x < +oo ) sin ( exp ( - x ) + exp ( - 4 x ) ) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Arthur Stroud, Don Secrest,
        //    Gaussian Quadrature Formulas,
        //    Prentice Hall, 1966,
        //    LC: QA299.4G3S7.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[N], the evaluation points.
        //
        //    Output, double P14_FUN[N], the function values.
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
            else if (-4.0 * x[i] < exponent_min)
            {
                f[i] = Math.Sin(Math.Exp(-x[i]));
            }
            else
            {
                f[i] = Math.Sin(Math.Exp(-x[i]) + Math.Exp(-4.0 * x[i]));
            }
        }

        return f;
    }

    public static string p14_title()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P14_TITLE returns the title for problem 14.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, string P14_TITLE, the title of the problem.
        //
    {
        string title;

        title = "sin ( exp(-x) + exp(-4x) )";

        return title;
    }

}