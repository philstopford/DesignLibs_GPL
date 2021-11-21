using System;
using Burkardt.Quadrature;

namespace Burkardt.IntegralNS;

public static class CauchyPrincipalValue
{
    public static double cpv(Func<double,double> f, double a, double b, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CPV estimates the Cauchy Principal Value of an integral.
        //
        //  Location:
        //
        //    http://people.sc.fsu.edu/~jburkardt/c_src/cauchy_principal_value/cauchy_principal_value.cpp
        //
        //  Discussion:
        //
        //    This function can be used to estimate the Cauchy Principal Value of
        //    a singular integral of the form 
        //      Integral f(t)/(t-x) dt 
        //    over an interval which includes the singularity point t=x.
        //
        //    Isolate the singularity at x in a symmetric interval of finite size delta:
        //
        //        CPV ( Integral ( a         <= t <= b         ) p(t) / ( t - x ) dt )
        //      =       Integral ( a         <= t <= x - delta ) p(t) / ( t - x ) dt
        //      + CPV ( Integral ( x - delta <= t <= x + delta ) p(t) / ( t - x ) dt )
        //      +       Integral ( x + delta <= t <= b         ) p(t) / ( t - x ) dt.
        //
        //    We assume the first and third integrals can be handled in the usual way.
        //    The second integral can be rewritten as
        //      Integral ( -1 <= s <= +1 ) ( p(s*delta+x) - p(x) ) / s ds
        //    and approximated by
        //      Sum ( 1 <= i <= N ) w(i) * ( p(xi*delta+x) - p(x) ) / xi(i)
        //    = Sum ( 1 <= i <= N ) w(i) * ( p(xi*delta+x) ) / xi(i)
        //    if we assume that N is even, so that coefficients of p(x) sum to zero.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Julian Noble,
        //    Gauss-Legendre Principal Value Integration,
        //    Computing in Science and Engineering,
        //    Volume 2, Number 1, January-February 2000, pages 92-95.
        //
        //  Parameters:
        //
        //    Input, double F ( double x ), the function that evaluates the
        //    integrand.
        //
        //    Input, double A, B, the endpoints of the symmetric interval, 
        //    which contains a singularity of the form 1/(X-(A+B)/2).
        //
        //    Input, int N, the number of Gauss points to use.
        //    N must be even.
        //
        //    Output, double CPV, the estimate for the Cauchy Principal Value.
        //
    {
        int i;
        //
        //  N must be even.
        //
        if (n % 2 != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("CPV - Fatal error!");
            Console.WriteLine("  N must be even.");
            return 1;
        }

        //
        //  Get the Gauss-Legendre rule.
        //
        double[] x = new double[n];
        double[] w = new double[n];

        LegendreQuadrature.legendre_set(n, ref x, ref w);
        //
        //  Estimate the integral.
        //
        double value = 0.0;
        for (i = 0; i < n; i++)
        {
            double x2 = ((1.0 - x[i]) * a
                         + (1.0 + x[i]) * b)
                        / 2.0;
            value += w[i] * f(x2) / x[i];
        }

        return value;
    }

}