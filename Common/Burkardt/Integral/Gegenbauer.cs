namespace Burkardt.IntegralNS
{
    public static partial class Integral
    {
        public static double gegenbauer_integral ( int p, double lambda )

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
            double s;

            if ( ( p % 2 ) == 0 )
            {
                s = Helpers.Gamma ( p / 2.0 + 0.5 ) * Helpers.Gamma ( lambda + 0.5 ) 
                    / Helpers.Gamma ( p / 2.0 + lambda + 1.0 );
            }
            else
            {
                s = 0.0;
            }

            return s;
        }
    }
}