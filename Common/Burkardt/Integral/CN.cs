using System;

namespace Burkardt.IntegralNS
{
    public static class CN_Geg
    {
        public static double cn_leg_monomial_integral ( int n, int[] expon )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CN_LEG_MONOMIAL_INTEGRAL: integral of monomial with Legendre weight on CN.
            //
            //  Discussion:
            //
            //    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
            //
            //      w(x) = 1.
            //
            //    value = integral ( CN ) product ( 1 <= i <= n ) x(I)^expon(i) dx(i)
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
            //    Input, int N, the spatial dimension.
            //
            //    Input, int EXPON(N), the exponents.
            //
            //    Output, double CN_LEG_MONOMIAL_INTEGRAL, the value of the integral.
            //
        {
            int i;
            double value;
            double value2;

            value = 1.0;
            for ( i = 0; i < n; i++ )
            {
                value2 = C1.c1_leg_monomial_integral ( expon[i] );
                value = value * value2;
            }

            return value;
        }
        
        public static double cn_geg_monomial_integral ( int n, double alpha, int[] expon )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CN_GEG_MONOMIAL_INTEGRAL: integral of monomial with Gegenbauer weight on CN.
            //
            //  Discussion:
            //
            //    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
            //
            //      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
            //
            //    with -1.0 < alpha.
            //
            //    value = integral ( CN ) 
            //      product ( 1 <= i <= n ) x(I)^expon(i) (1-x(i)^2)^alpha dx(i)
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
            //    Input, int N, the spatial dimension.
            //
            //    Input, double ALPHA, the exponent of (1-X).
            //    -1.0 < ALPHA.
            //
            //    Input, int EXPON[N], the exponents.
            //
            //    Output, double CN_GEG_MONOMIAL_INTEGRA, the value of the integral.
            //
        {
            int i;
            double value;
            double value2;

            if ( alpha <= -1.0 )
            {
                Console.WriteLine("");
                Console.WriteLine("CN_GEG_MONOMIAL_INTEGRAL - Fatal error!");
                Console.WriteLine("  ALPHA <= -1.0");
                return ( 1 );
            }

            value = 1.0;
            for ( i = 0; i < n; i++ )
            {
                value2 = C1.c1_geg_monomial_integral ( alpha, expon[i] );
                value = value * value2;
            }

            return value;
        }

        public static double cn_jac_monomial_integral ( int n, double alpha, double beta, 
                int[] expon )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CN_JAC_MONOMIAL_INTEGRAL: integral of a monomial with Jacobi weight over CN.
            //
            //  Discussion:
            //
            //    value = integral ( CN ) 
            //      product ( 1 <= i <= n ) x(I)^expon(i) (1-x(i))^alpha (1+x(i))^beta dx(i)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    26 January 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the spatial dimension.
            //
            //    Input, double ALPHA, the exponent of (1-X) in the weight factor.
            //
            //    Input, double BETA, the exponent of (1+X) in the weight factor.
            //
            //    Input, int EXPON[N], the exponents.
            //
            //    Output, double CN_JAC_MONOMIAL_INTEGRAL, the value of the integral.
            //
        {
            int i;
            double value;
            double value2;

            value = 1.0;
            for ( i = 0; i < n; i++ )
            {
                value2 = C1.c1_jac_monomial_integral ( alpha, beta, expon[i] );
                value = value * value2;
            }

            return value;
        }    }
}