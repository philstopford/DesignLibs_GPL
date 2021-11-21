using System;
using Burkardt.Types;

namespace Burkardt.IntegralNS;

public static partial class Integral
{
    public static double j_double_product_integral(int i, int j, double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    J_DOUBLE_PRODUCT_INTEGRAL: integral of J(i,x)*J(j,x)*(1-x)^a*(1+x)^b.
        //
        //  Discussion:
        //
        //    VALUE = integral ( -1 <= x <= +1 ) J(i,x)*J(j,x)*(1-x)^a*(1+x)^b dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, J, the polynomial indices.
        //
        //    Input, double A, B, the parameters.
        //    -1 < A, B.
        //
        //    Output, double VALUE, the value of the integral.
        //
    {
        double value;

        if (i != j)
        {
            value = 0.0;
        }
        else
        {
            double i_r8 = i;

            value = Math.Pow(2, a + b + 1.0)
                    / (2.0 * i_r8 + a + b + 1.0)
                    * Helpers.Gamma(i_r8 + a + 1.0)
                    * Helpers.Gamma(i_r8 + b + 1.0)
                    / typeMethods.r8_factorial(i)
                    / Helpers.Gamma(i_r8 + a + b + 1.0);
        }

        return value;
    }
    public static double j_integral(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    J_INTEGRAL evaluates a monomial integral associated with J(n,a,b,x).
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      integral ( -1 <= x < +1 ) x^n dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the exponent.
        //    0 <= N.
        //
        //    Output, double J_INTEGRAL, the value of the integral.
        //
    {
        double value = (n % 2) switch
        {
            1 => 0.0,
            _ => 2.0 / (n + 1)
        };

        return value;
    }

    public static double jacobi_integral ( int expon, double alpha, double beta )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JACOBI_INTEGRAL evaluates the integral of a monomial with Jacobi weight.
        //
        //  Discussion:
        //
        //    VALUE = Integral ( -1 <= X <= +1 ) x^EXPON (1-x)^ALPHA (1+x)^BETA dx
        //
        //  Modified:
        //
        //    08 September 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int EXPON, the exponent.
        //
        //    Input, double ALPHA, the exponent of (1-X) in the weight factor.
        //
        //    Input, double BETA, the exponent of (1+X) in the weight factor.
        //
        //    Output, double JACOBI_INTEGRAL, the value of the integral.
        //
    {
        double c = expon;

        double s = (expon % 2) switch
        {
            0 => +1.0,
            _ => -1.0
        };

        double arg1 = - alpha;
        double arg2 = 1.0 + c;
        double arg3 = 2.0 + beta + c;
        double arg4 = - 1.0;

        double value1 = typeMethods.r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

        arg1 = - beta;
        arg2 =   1.0 + c;
        arg3 =   2.0 + alpha + c;
        arg4 = - 1.0;

        double value2 = typeMethods.r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

        double value = typeMethods.r8_gamma ( 1.0 + c ) * ( 
            s * typeMethods.r8_gamma ( 1.0 + beta  ) * value1 
            / typeMethods.r8_gamma ( 2.0 + beta  + c ) 
            +     typeMethods.r8_gamma ( 1.0 + alpha ) * value2 
            / typeMethods.r8_gamma ( 2.0 + alpha + c ) );

        return value;
    }
}