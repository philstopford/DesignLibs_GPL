﻿using System;

namespace Burkardt.Laguerre;

public static partial class Integrands
{
    public static double p08_alpha()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P08_ALPHA returns ALPHA for problem 8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double P08_ALPHA, the value of ALPHA.
        //
    {
        const double alpha = 0.0;

        return alpha;
    }

    public static double p08_exact()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P08_EXACT returns the estimated integral for problem 8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double P08_EXACT, the estimated value of the integral.
        //
    {
        const double exact = Math.PI * Math.PI / 6.0;

        return exact;
    }

    public static double[] p08_fun(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P08_FUN evaluates the integrand for problem 8.
        //
        //  Integral:
        //
        //    Integral ( 0 <= x < +oo ) x / ( exp ( x ) - 1 ) dx
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
        //    Output, double P08_FUN[N], the function values.
        //
    {
        const double exponent_max = 80.0;
        int i;

        double[] f = new double[n];

        for (i = 0; i < n; i++)
        {
            switch (x[i])
            {
                case 0.0:
                    f[i] = 1.0 / Math.Exp(x[i]);
                    break;
                default:
                {
                    if (x[i] < exponent_max)
                    {
                        f[i] = x[i] / (Math.Exp(x[i]) - 1.0);
                    }
                    else
                    {
                        f[i] = 0.0;
                    }

                    break;
                }
            }
        }

        return f;
    }

    public static string p08_title()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P08_TITLE returns the title for problem 8.
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
        //    Output, string P08_TITLE, the title of the problem.
        //
    {
        const string title = "Debye function";

        return title;
    }

}