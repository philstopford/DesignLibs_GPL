using System;

namespace Burkardt.IntegralNS;

public static partial class Integral
{
    public static double chebyshev2_integral ( int expon )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV2_INTEGRAL evaluates a monomial Chebyshev type 2 integral.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      integral ( -1 <= x <= +1 ) x^n * sqrt ( 1 - x^2 ) dx
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
        //    Output, double CHEBYSHEV2_INTEGRAL, the value of the exact integral.
        //
    {
        double bot;
        double exact;
        int i;
            
        double top;
        switch (expon % 2)
        {
            //
            //  Get the exact value of the integral.
            //
            case 0:
            {
                top = 1;
                bot = 1;
                for ( i = 2; i <= expon; i += 2 )
                {
                    top *= ( i - 1 );
                    bot *= i;
                }

                bot *= expon + 2;

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