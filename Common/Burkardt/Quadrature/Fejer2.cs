using System;

namespace Burkardt.Quadrature
{
    public static class Fejer2
    {
        public static void fejer2_compute_np ( int n, int np, double[] p, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEJER2_COMPUTE_NP computes a Fejer type 2 rule.
        //
        //  Discussion:
        //
        //    Our convention is that the abscissas are numbered from left to right.
        //
        //    The rule is defined on [-1,1].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //    1 <= N.
        //
        //    Input, int NP, the number of parameters.
        //
        //    Input, double P[NP], parameters which are not needed by this function.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
        {
            fejer2_compute ( n, ref x, ref w );
        }
        
        public static void fejer2_compute ( int n, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEJER2_COMPUTE computes a Fejer type 2 rule.
        //
        //  Discussion:
        //
        //    Our convention is that the abscissas are numbered from left to right.
        //
        //    The rule is defined on [-1,1].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //    1 <= N.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
        {
            int i;
            int j;
            double p;
            double pi = 3.141592653589793;
            double theta;

            if ( n < 1 )
            {
                Console.WriteLine("");
                Console.WriteLine("FEJER2_COMPUTE - Fatal error!");
                Console.WriteLine("  Illegal value of N = " + n + "");
                return;
            }
            else if ( n == 1 )
            {
                x[0] = 0.0;
                w[0] = 2.0;
                return;
            }

            for ( i = 0; i < n; i++ )
            {
                x[i] =  Math.Cos ( ( double ) ( n - i ) * pi
                                   / ( double ) ( n + 1 ) );
            }
            if ( ( n % 2 ) == 1 )
            {
                x[(n-1)/2] = 0.0;
            }

            if ( n == 2 )
            {
                w[0] = 1.0;
                w[1] = 1.0;
            }
            else
            {
                for ( i = 0; i < n; i++ )
                {
                    theta = ( double ) ( n - i ) * pi
                            / ( double ) ( n + 1 );

                    w[i] = 1.0;

                    for ( j = 1; j <= ( ( n - 1 ) / 2 ); j++ )
                    {
                        w[i] = w[i] - 2.0 *  Math.Cos ( 2.0 * ( double ) ( j ) * theta )
                            / ( double ) ( 4 * j * j - 1 );
                    }
                    p = 2.0 * ( double ) ( ( ( n + 1 ) / 2 ) ) - 1.0;
                    w[i] = w[i] -  Math.Cos ( ( p + 1.0 ) * theta ) / p;
                }
                for ( i = 0; i < n; i++ )
                {
                    w[i] = 2.0 * w[i] / ( double ) ( n + 1 );
                }
            }
        }
    }
}