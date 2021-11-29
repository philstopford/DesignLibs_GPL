using System;

namespace Burkardt.Laguerre;

public static partial class Integrands
{
    public static double p19_alpha()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P19_ALPHA returns ALPHA for problem 19.
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
        //    Output, double P19_ALPHA, the value of ALPHA.
        //
    {
        const double alpha = 0.0;

        return alpha;
    }

    public static double p19_exact()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P19_EXACT returns the exact integral for problem 19.
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
        //    Output, double P19_EXACT, the value of the integral.
        //
    {
        const double beta = 0.5;


        double exact = (1.0 - beta) * Math.PI
                       / (Math.Pow(10.0, beta) * Math.Sin(Math.PI * beta));

        return exact;
    }

    public class p19Data
    {
        public double beta = 0.5;
            
    }
    public static double[] p19_fun(ref p19Data data, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P19_FUN evaluates the integrand for problem 19.
        //
        //  Integral:
        //
        //    Integral ( 0 <= x < +oo ) x^(alpha-1) / ( 1 + 10 x )^2 dx
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
        //    Output, double P61_FUN[N], the integrand values.
        //
    {
        int i;

        double[] fx = new double[n];

        for (i = 0; i < n; i++)
        {
            fx[i] = data.beta switch
            {
                1.0 => 1.0 / Math.Pow(1.0 + 10.0 * x[i], 2),
                < 1.0 when x[i] == 0.0 => 0.0,
                _ => Math.Pow(x[i], data.beta - 1.0) / Math.Pow(1.0 + 10.0 * x[i], 2)
            };
        }

        return fx;
    }

    public static string p19_title()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P19_TITLE returns the title for problem 19.
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
        //    Output, string P19_TITLE, the title of the problem.
        //
    {
        const string title = "x^(beta-1) / ( 1 + 10 x )^2";

        return title;
    }


}