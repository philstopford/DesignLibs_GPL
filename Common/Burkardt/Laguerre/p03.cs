using System;

namespace Burkardt.Laguerre;

public static partial class Integrands
{
    public static double p03_alpha()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P03_ALPHA returns ALPHA for problem 3.
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
        //    Output, double P03_ALPHA, the value of ALPHA.
        //
    {
        double alpha;

        alpha = 2.0;

        return alpha;
    }

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

        exact = 13.628;

        return exact;
    }

    public static double[] p03_fun(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P03_FUN evaluates the integrand for problem 3.
        //
        //  Discussion:
        //
        //    D&R gives "exact" value as 13.628...
        //    Mathematica returns        13.440045415012575106...
        //    D&R gives Laguerre(16) as   0.44996932...
        //
        //  Integral:
        //
        //    exp ( -2 ) Integral ( 2 <= x < +oo ) / ( x^1.01 ) dx
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
        //  Reference:
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Dover, 2007,
        //    ISBN: 0486453391,
        //    LC: QA299.3.D28.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[N], the evaluation points.
        //
        //    Output, double P03_FUN[N], the function values.
        //
    {
        double[] f;
        int i;

        f = new double[n];

        for (i = 0; i < n; i++)
        {
            f[i] = Math.Exp(-2.0) * 1.0 / Math.Pow(x[i], 1.01);
        }

        return f;
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
        //    28 July 2007
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

        title = "1 / ( x^1.01 )";

        return title;
    }

}