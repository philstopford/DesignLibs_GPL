using Burkardt.Types;

namespace Burkardt.IntegralNS
{
    public static partial class Integral
    {
        public static double laguerre_integral ( int p )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LAGUERRE_INTEGRAL evaluates a monomial Laguerre integral.
            //
            //  Discussion:
            //
            //    The integral being computed is
            //
            //      integral ( 0 <= x < +oo ) x^p * exp ( -x ) dx
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 May 2014
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
            //    Output, double LAGUERRE_INTEGRAL, the value of the integral.
            //
        {
            double s;

            s = typeMethods.r8_factorial ( p );

            return s;
        }
    }
}