using System;

namespace Burkardt.IntegralNS;

public static partial class Integral
{
    public static double chebyshev1_integral ( int expon )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV1_INTEGRAL evaluates a monomial Chebyshev type 1 integral.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      integral ( -1 <= x <= +1 ) x^n / sqrt ( 1 - x^2 ) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int EXPON, the exponent.
        //
        //    Output, double CHEBYSHEV1_INTEGRAL, the value of the exact integral.
        //
    {
        double exact;

        switch (expon % 2)
        {
            //
            //  Get the exact value of the integral.
            //
            case 0:
            {
                double top = 1;
                double bot = 1;
                int i;
                for ( i = 2; i <= expon; i += 2 )
                {
                    top *= i - 1;
                    bot *= i;
                }
	
                exact = Math.PI * top / bot;
                break;
            }
            default:
                exact = 0.0;
                break;
        }

        return exact;
    }
}