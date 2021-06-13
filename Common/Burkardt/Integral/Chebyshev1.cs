namespace Burkardt.IntegralNS
{
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
            double bot;
            double exact;
            int i;
            const double r8_pi = 3.141592653589793;
            double top;
            //
            //  Get the exact value of the integral.
            //
            if ( ( expon % 2 ) == 0 )
            {
                top = 1;
                bot = 1;
                for ( i = 2; i <= expon; i = i + 2 )
                {
                    top = top * ( i - 1 );
                    bot = bot *   i;
                }
	
                exact = r8_pi * ( double ) ( top ) / ( double ) ( bot );
            }
            else
            {
                exact = 0.0;	
            }

            return exact;
        }
    }
}