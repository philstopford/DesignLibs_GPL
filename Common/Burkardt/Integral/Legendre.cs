namespace Burkardt.IntegralNS
{
    public static partial class Integral
    {
        public static double legendre_integral ( int p )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEGENDRE_INTEGRAL evaluates a monomial Legendre integral.
            //
            //  Discussion:
            //
            //    Integral ( -1 <= x <= +1 ) x^p dx
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 February 2008
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
            //    Output, double LEGENDRE_INTEGRAL, the value of the exact integral.
            //
        {
            double s;

            if ( ( p % 2 ) == 0 )
            {
                s = 2.0 / ( double ) ( p + 1 );
            }
            else
            {
                s = 0.0;
            }

            return s;
        }
    }
}