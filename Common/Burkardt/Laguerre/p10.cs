using System;

namespace Burkardt.Laguerre;

public static partial class Integrands
{
    public static double p10_alpha()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P10_ALPHA returns ALPHA for problem 10.
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
        //    Output, double P10_ALPHA, the value of ALPHA.
        //
    {
        const double alpha = 0.0;

        return alpha;
    }

    public static double p10_exact()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P10_EXACT returns the estimated integral for problem 10.
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
        const double exact = Math.PI / 2.0;

        return exact;
    }

    public static double[] p10_fun(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P10_FUN evaluates the integrand for problem 10.
        //
        //  Discussion:
        //
        //    S&S gives exact value as pi/2 = 1.5707963267948966192...
        //    S&S gives Laguerre(16) as       1.5537377347...
        //    S&S gives EXP_TRANSFORM(16) as  1.4293043007...
        //
        //  Integral:
        //
        //    Integral ( 0 <= x < +oo ) 1/(1+x*x) dx
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
        //    Output, double P10_FUN[N], the function values.
        //
    {
        int i;

        double[] f = new double[n];

        for (i = 0; i < n; i++)
        {
            f[i] = 1.0 / (1.0 + x[i] * x[i]);
        }

        return f;
    }

    public static string p10_title()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P10_TITLE returns the title for problem 10.
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
        //    Output, string P10_TITLE, the title of the problem.
        //
    {
        const string title = "1 / ( 1 + x*x )";

        return title;
    }
}