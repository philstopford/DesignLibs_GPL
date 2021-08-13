using System;

namespace Burkardt.PolynomialNS
{
    public static class Meixner
    {
        public static void meixner ( int n, double beta, double c, double x, ref double[] v )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MEIXNER evaluates Meixner polynomials at a point.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 March 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Walter Gautschi,
            //    Orthogonal Polynomials: Computation and Approximation,
            //    Oxford, 2004,
            //    ISBN: 0-19-850672-4,
            //    LC: QA404.5 G3555.
            //
            //  Parameters:
            //
            //    Input, int N, the maximum order of the polynomial.  
            //    N must be at least 0.
            //
            //    Input, double BETA, the Beta parameter.  0 < BETA.
            //
            //    Input, double C, the C parameter.  0 < C < 1.
            //
            //    Input, double X, the evaluation point.
            //
            //    Output, double V[N+1], the value of the polynomials at X.
            //
        {
            int i;

            if ( beta <= 0.0 )
            {
                Console.WriteLine("");
                Console.WriteLine("MEIXNER - Fatal error!");
                Console.WriteLine("  Parameter BETA must be positive.");
                return;
            }

            if ( c <= 0.0 || 1.0 <= c )
            {
                Console.WriteLine("");
                Console.WriteLine("MEIXNER - Fatal error!");
                Console.WriteLine("  Parameter C must be strictly between 0 and 1.");
                return;
            }

            if ( n < 0 )
            {
                Console.WriteLine("");
                Console.WriteLine("MEIXNER - Fatal error!");
                Console.WriteLine("  Parameter N must be nonnegative.");
                return;
            }

            v[0] = 1.0;

            if ( n == 0 )
            {
                return;
            }

            v[1] = ( c - 1.0 ) * x / beta / c + 1.0;

            if ( n == 1 )
            {
                return;
            }

            for ( i = 1; i < n; i++ )
            {
                v[i+1] = ( 
                    ( ( c - 1.0 ) * x + ( 1.0 + c ) 
                        * ( double ) ( i ) + beta * c ) * v[i]
                    - ( double ) ( i ) * v[i-1] 
                ) / ( ( double ) ( i ) + beta );
            }

        }

    }
}