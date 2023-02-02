﻿using System;
using Burkardt.Types;

namespace Burkardt.Laguerre;

public static partial class Integrands
{
    public static double p20_alpha()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P20_ALPHA returns ALPHA for problem 20.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double P20_ALPHA, the value of ALPHA.
        //
    {
        const double alpha = 0.0;

        return alpha;
    }

    public class p20Data
    {
        public const double beta = 1.0;
    }

    public static double p20_exact(ref p20Data data)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P20_EXACT returns the exact integral for problem 20.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double P20_EXACT, the value of the integral.
        //
    {
        double exact = (
            Math.Log(1.5) / Math.Pow(2.0, p20Data.beta)
            - 1.0 / Math.Pow(2.0, p20Data.beta + 1.0) *
            Math.Log((16.0 + Math.Pow(0.25, p20Data.beta)) / (1.0 + Math.Pow(0.25, p20Data.beta)))
            - Math.Atan(Math.Pow(2.0, p20Data.beta + 2.0)) - Math.Atan(Math.Pow(2.0, p20Data.beta))
        ) / (1.0 + Math.Pow(0.25, p20Data.beta));

        return exact;
    }

    public static double[] p20_fun(ref p20Data data, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P20_FUN evaluates the integrand for problem 20.
        //
        //  Integral:
        //
        //    Integral ( 0 <= x < +oo ) 
        //      1 / ( 2^beta * ( ( x - 1 )^2 + (1/4)^beta ) * ( x - 2 ) ) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Robert Piessens, Elise de Doncker-Kapenga, 
        //    Christian Ueberhuber, David Kahaner,
        //    QUADPACK: A Subroutine Package for Automatic Integration,
        //    Springer, 1983, page 84.
        //
        //  Parameters:
        //
        //    Input, int N, the number of evaluation points.
        //
        //    Input, double X[N], the evaluation points.
        //
        //    Output, double P20_FUN[N], the integrand values.
        //
    {
        int i;

        double[] fx = new double[n];

        for (i = 0; i < n; i++)
        {
            if (Math.Pow(x[i] - 1.0, 2) + Math.Pow(0.25, p20Data.beta) == 0.0 || Math.Abs(x[i] - 2.0) <= typeMethods.r8_epsilon())
            {
                fx[i] = 0.0;
            }
            else
            {
                fx[i] = 1.0 /
                        (Math.Pow(2.0, p20Data.beta)
                         * (Math.Pow(x[i] - 1.0, 2) + Math.Pow(0.25, p20Data.beta))
                         * (x[i] - 2.0));
            }
        }

        return fx;
    }

    public static string p20_title()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P20_TITLE returns the title for problem 20.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, string P20_TITLE, the title of the problem.
        //
    {
        const string title = "1 / ( 2^beta * ( ( x - 1 )^2 + (1/4)^beta ) * ( x - 2 ) )";

        return title;
    }


}