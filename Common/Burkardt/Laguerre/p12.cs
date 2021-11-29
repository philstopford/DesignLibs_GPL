using System;

namespace Burkardt.Laguerre;

public static partial class Integrands
{
    public static double p12_alpha()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P12_ALPHA returns ALPHA for problem 12.
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
        //    Output, double P12_ALPHA, the value of ALPHA.
        //
    {
        const double alpha = 0.0;

        return alpha;
    }

    public static double p12_exact()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P12_EXACT returns the estimated integral for problem 12.
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
        const double exact = 0.5;

        return exact;
    }

    public static double[] p12_fun(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P12_FUN evaluates the integrand for problem 12.
        //
        //  Discussion:
        //
        //    S&S gives exact value as Math.PI =  0.5
        //    S&S gives Laguerre(16) as      0.5000000000...
        //    S&S gives EXP_TRANSFORM(16) as 0.5019065783... 
        //
        //  Integral:
        //
        //    Integral ( 0 <= x < +oo ) exp ( -x ) * cos ( x ) dx
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
        //    Output, double P12_FUN[N], the function values.
        //
    {
        const double exponent_min = -80.0;
        int i;

        double[] f = new double[n];

        for (i = 0; i < n; i++)
        {
            if (-x[i] < exponent_min)
            {
                f[i] = 0.0;
            }
            else
            {
                f[i] = Math.Exp(-x[i]) * Math.Cos(x[i]);
            }
        }

        return f;
    }

    public static string p12_title()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P12_TITLE returns the title for problem 12.
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
        //    Output, string P12_TITLE, the title of the problem.
        //
    {
        const string title = "exp ( - x ) * cos ( x )";

        return title;
    }

}