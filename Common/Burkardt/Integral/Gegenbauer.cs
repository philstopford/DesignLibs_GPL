using System;

namespace Burkardt.IntegralNS;

using Chebyshev = Burkardt.ChebyshevNS.Chebyshev;

public static partial class Integral
{
    public static double gegenbauer_integral(int p, double lambda)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEGENBAUER_INTEGRAL evaluates a monomial integral with Gegenbauer weight.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      integral ( -1 <= x < +1 ) x^p * ( 1 - x^2 )^(lambda-1/2) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P, the exponent.
        //    0 <= P.
        //
        //    Input, double LAMBDA, the exponent term.
        //    -1/2 < LAMBDA.
        //
        //    Output, real GEGENBAUER_INTEGRAL, the value of the integral.
        //
    {
        double s = (p % 2) switch
        {
            0 => Helpers.Gamma(p / 2.0 + 0.5) * Helpers.Gamma(lambda + 0.5) / Helpers.Gamma(p / 2.0 + lambda + 1.0),
            _ => 0.0
        };

        return s;
    }

    public static double gegenbauer_cc1(int n, double lambda, Func<double, double> f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEGENBAUER_CC1 estimates the Gegenbauer integral of a function.
        //
        //  Discussion:
        //
        //     value = integral ( -1 <= x <= + 1 ) ( 1 - x^2 )^(lambda-1/2) * f(x) dx
        //
        //     The approximation uses the practical abscissas, that is, the extreme
        //     points of Tn(x).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    D B Hunter, H V Smith,
        //    A quadrature formula of Clenshaw-Curtis type for the Gegenbauer 
        //    weight function,
        //    Journal of Computational and Applied Mathematics,
        //    Volume 177, 2005, pages 389-400.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points to use.
        //    1 <= N.
        //
        //    Input, double LAMBDA, used in the exponent of (1-x^2).
        //    -0.5 < LAMBDA.
        //
        //    Input, double, external, F(x), the function to be integrated 
        //    with the Gegenbauer weight.
        //
        //    Output, double WEIGHT, the estimate for the Gegenbauer 
        //    integral of F.
        //
    {
        double[] a2;
        int rh;
        int s;
        int sigma;
        double u;
        double value = 0;

        value = 0.0;

        s = n / 2;
        sigma = n % 2;

        a2 = Chebyshev.chebyshev_even1(n, f);

        rh = s;
        u = 0.5 * (sigma + 1) * a2[rh];
        for (rh = s - 1; 1 <= rh; rh--)
        {
            u = (rh - lambda)
                / (rh + lambda + 1.0) * u + a2[rh];
        }

        u = -lambda * u / (lambda + 1.0) + 0.5 * a2[0];

        value = Helpers.Gamma(lambda + 0.5) * Math.Sqrt(Math.PI) * u
                / Helpers.Gamma(lambda + 1.0);

        return value;
    }

    public static double gegenbauer_cc2(int n, double lambda, Func<double, double> f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEGENBAUER_CC2 estimates the Gegenbauer integral of a function.
        //
        //  Discussion:
        //
        //     value = integral ( -1 <= x <= + 1 ) ( 1 - x^2 )^(lambda-1/2) * f(x) dx
        //
        //     The approximation uses the classical abscissas, that is, the zeros
        //     of Tn(x).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    D B Hunter, H V Smith,
        //    A quadrature formula of Clenshaw-Curtis type for the Gegenbauer 
        //    weight function,
        //    Journal of Computational and Applied Mathematics,
        //    Volume 177, 2005, pages 389-400.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points to use.
        //    1 <= N.
        //
        //    Input, double LAMBDA, used in the exponent of (1-x^2).
        //    -0.5 < LAMBDA.
        //
        //    Input, double F( double x ), the function to be integrated with
        //    the Gegenbauer weight.
        //
        //    Output, double WEIGHT, the estimate for the Gegenbauer 
        //    integral of F.
        //
    {
        double[] b2;
        int rh;
        int s;
        int sigma;
        double u;
        double value = 0;

        value = 0.0;

        s = n / 2;
        sigma = n % 2;

        b2 = Chebyshev.chebyshev_even2(n, f);

        rh = s;
        u = (sigma + 1) * b2[rh];
        for (rh = s - 1; 1 <= rh; rh--)
        {
            u = (rh - lambda)
                / (rh + lambda + 1.0) * u + b2[rh];
        }

        u = -lambda * u / (lambda + 1.0) + 0.5 * b2[0];

        value = Helpers.Gamma(lambda + 0.5) * Math.Sqrt(Math.PI) * u
                / Helpers.Gamma(lambda + 1.0);

        return value;
    }
}