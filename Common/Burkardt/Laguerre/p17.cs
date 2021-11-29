using System;

namespace Burkardt.Laguerre;

public static partial class Integrands
{
    public static double p17_alpha()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P17_ALPHA returns ALPHA for problem 17.
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
        //    Output, double P17_ALPHA, the value of ALPHA.
        //
    {
        const double alpha = 0.0;

        return alpha;
    }

    public static double p17_exact()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P17_EXACT returns the exact integral for problem 17.
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
        //    Output, double P17_EXACT, the value of the integral.
        //
    {
        const double beta = 2.0;


        double exact = Math.Sqrt(Math.PI) * Math.Cos(0.5 * Math.Atan(Math.Pow(2.0, beta)))
                       / Math.Sqrt(Math.Sqrt(1.0 + Math.Pow(0.25, beta)));

        return exact;
    }

    public class p17Data
    {
        public double beta = 2.0;

    }
        
    public static double[] p17_fun(ref p17Data data, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P17_FUN evaluates the integrand for problem 17.
        //
        //  Integral:
        //
        //    Integral ( 0 <= x < +oo) exp ( - x / 2^beta ) * cos ( x ) / Math.Sqrt ( x ) dx
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
        //    Output, double P17_FUN[N], the integrand values.
        //
    {
        int i;

        double[] fx = new double[n];

        for (i = 0; i < n; i++)
        {
            fx[i] = x[i] switch
            {
                0.0 => 0.0,
                _ => Math.Exp(-x[i] / Math.Pow(2.0, data.beta)) * Math.Cos(x[i]) / Math.Sqrt(x[i])
            };
        }

        return fx;
    }

    public static string p17_title()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P17_TITLE returns the title for problem 17.
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
        //    Output, string P17_TITLE, the title of the problem.
        //
    {
        const string title = "exp ( - x / 2^beta ) * cos ( x ) / Math.Sqrt ( x )";

        return title;
    }
}