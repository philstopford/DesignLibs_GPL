using System;
using Burkardt.Types;

namespace Burkardt.IntegralNS;

public static class C1
{
    public static double c1_geg_monomial_integral ( double alpha, int expon )

//****************************************************************************80
//
//  Purpose:
//
//    C1_GEG_MONOMIAL_INTEGRAL: integral of monomial with Gegenbauer weight on C1.
//
//  Discussion:
//
//    C1_GEG is the interval [-1,+1] with the Gegenbauer weight function
//
//      w(alpha;x) = (1-x^2)^alpha
//
//    with -1.0 < alpha.
//
//    value = integral ( -1 <= x <= +1 ) x^expon (1-x^2)^alpha dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the exponent of (1-X^2).
//    - 1.0 < ALPHA.
//
//    Input, int EXPON, the exponent.
//    0 <= EXPON.
//
//    Output, double C1_GEG_MONOMIAL_INTEGRAL, the value of the integral.
//
    {
        double value;

        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("C1_GEG_MONOMIAL_INTEGRAL - Fatal error!");
                Console.WriteLine("  ALPHA <= -1.0");
                return 1;
        }

        switch (expon % 2)
        {
            case 1:
                value = 0.0;
                return value;
        }

        double c = expon;

        double arg1 = - alpha;
        double arg2 = 1.0 + c;
        double arg3 = 2.0 + alpha + c;
        double arg4 = - 1.0;

        double value1 = typeMethods.r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

        value = 2.0 * typeMethods.r8_gamma ( 1.0 + c ) * typeMethods.r8_gamma ( 1.0 + alpha ) 
            * value1 / typeMethods.r8_gamma ( 2.0 + alpha + c );

        return value;
    }

    public static double c1_jac_monomial_integral ( double alpha, double beta, int expon )

//****************************************************************************80
//
//  Purpose:
//
//    C1_JAC_MONOMIAL_INTEGRAL: integral of a monomial with Jacobi weight over C1.
//
//  Discussion:
//
//    value = integral ( -1 <= x <= +1 ) x^expon (1-x)^alpha (1+x)^beta dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the exponent of (1-X) in the weight factor.
//
//    Input, double BETA, the exponent of (1+X) in the weight factor.
//
//    Input, int EXPON, the exponent.
//
//    Output, double C1_JAC_MONOMIAL_INTEGRAL, the value of the integral.
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

    public static double c1_leg_monomial_integral ( int expon )

//****************************************************************************80
//
//  Purpose:
//
//    C1_LEG_MONOMIAL_INTEGRAL: integral of monomial with Legendre weight on C1.
//
//  Discussion:
//
//    C1_LEG is the interval [-1,+1] with the Legendre weight function
//
//      w(x) = 1.
//
//    value = integral ( -1 <= x <= +1 ) x^expon dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON, the exponent.
//    0 <= EXPON.
//
//    Output, double C1_LEG_MONOMIAL_INTEGRAL, the value of the integral.
//
    {
        double value = 0;

        switch (expon)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("C1_LEG_MONOMIAL_INTEGRAL - Fatal error!");
                Console.WriteLine("  EXPON < 0.");
                return 1;
        }

        switch (expon % 2)
        {
            case 1:
                value = 0.0;
                return value;
            default:
                value = 2.0 / (expon + 1);

                return value;
        }
    }

}