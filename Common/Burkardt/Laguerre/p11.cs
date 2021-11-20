﻿using System;

namespace Burkardt.Laguerre;

public static partial class Integrands
{
    public static double p11_alpha()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P11_ALPHA returns ALPHA for problem 11.
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
        //    Output, double P11_ALPHA, the value of ALPHA.
        //
    {
        const double alpha = 0.0;

        return alpha;
    }

    public static double p11_exact()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P11_EXACT returns the estimated integral for problem 11.
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
        const double exact = Math.PI;

        return exact;
    }

    public static double[] p11_fun(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P11_FUN evaluates the integrand for problem 11.
        //
        //  Discussion:
        //
        //    S&S gives exact value as Math.PI =  3.1415926535897932385...
        //    S&S gives Laguerre(16) as      2.6652685196...
        //    S&S gives EXP_TRANSFORM(16) as 2.3629036166... 
        //
        //  Integral:
        //
        //    Integral ( 0 <= x < +oo ) 1/((1+x)*sqrt(x)) dx
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
        //    Output, double P11_FUN[N], the function values.
        //
    {
        int i;

        double[] f = new double[n];

        for (i = 0; i < n; i++)
        {
            f[i] = x[i] switch
            {
                0.0 => 0.0,
                _ => 1.0 / ((1.0 + x[i]) * Math.Sqrt(x[i]))
            };
        }

        return f;
    }

    public static string p11_title()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P11_TITLE returns the title for problem 11.
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
        //    Output, string P11_TITLE, the title of the problem.
        //
    {
        const string title = "1 / ( (1+x) * sqrt(x) )";

        return title;
    }

}