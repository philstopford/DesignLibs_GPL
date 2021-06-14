using System;
using Burkardt.Types;

namespace Burkardt.IntegralNS
{
    public static partial class Integral
    {
        public static double hermite_integral ( int p )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_INTEGRAL evaluates a monomial Hermite integral.
            //
            //  Discussion:
            //
            //    Integral ( -oo < x < oo ) x^p exp(-x^2) dx
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
            //    Input, int P, the exponent of the monomial.  
            //    0 <= P.
            //
            //    Output, double HERMITE_INTEGRAL, the value of the integral.
            //
        {
            
            double value;

            if ( ( p % 2 ) == 0 )
            {
                value = typeMethods.r8_factorial2 ( p - 1 ) * Math.Sqrt ( Math.PI ) / Math.Pow ( 2.0, p / 2 );
            }
            else
            {
                value = 0.0;
            }
            return value;
        }
    }
}