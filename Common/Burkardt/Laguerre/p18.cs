﻿using System;

namespace Burkardt.Laguerre;

public static partial class Integrands
{
    public static double p18_alpha()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P18_ALPHA returns ALPHA for problem 18.
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
        //    Output, double P18_ALPHA, the value of ALPHA.
        //
    {
        const double alpha = 0.0;

        return alpha;
    }

    public class p18Data
    {
        public double beta = 2.0;

    }
        
    public static double p18_exact(ref p18Data data)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P18_EXACT returns the exact integral for problem 18.
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
        //    Output, double P18_EXACT, the value of the integral.
        //
    {
        double exact = Math.Pow(2.0, 3.0 * data.beta + 1.0);

        return exact;
    }

    public static double[] p18_fun(ref p18Data data, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P18_FUN evaluates the integrand for problem 18.
        //
        //  Integral:
        //
        //    Integral ( 0 <= x < +oo ) x^2 * exp ( - x / 2^beta ) dx
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
        //    Output, double P18_FUN[N], the integrand values.
        //
    {
        int i;

        double[] fx = new double[n];

        for (i = 0; i < n; i++)
        {
            fx[i] = x[i] * x[i] * Math.Exp(-x[i] / Math.Pow(2, data.beta));
        }

        return fx;
    }

    public static string p18_title()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P18_TITLE returns the title for problem 18.
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
        //    Output, string P18_TITLE, the title of the problem.
        //
    {
        const string title = "x^2 * exp ( - x / 2^beta )";

        return title;
    }
}