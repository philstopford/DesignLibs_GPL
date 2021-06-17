using System;

namespace Burkardt
{
    public static class Monomial
    {
        public static double[] monomial_value ( int m, int n, int[] e, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONOMIAL_VALUE evaluates a monomial.
        //
        //  Discussion:
        //
        //    This routine evaluates a monomial of the form
        //
        //      product ( 1 <= i <= m ) x(i)^e(i)
        //
        //    where the exponents are nonnegative integers.  Note that
        //    if the combination 0^0 is encountered, it should be treated
        //    as 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N, the number of points at which the
        //    monomial is to be evaluated.
        //
        //    Input, int E[M], the exponents.
        //
        //    Input, double X[M*N], the point coordinates.
        //
        //    Output, double MONOMIAL_VALUE[N], the value of the monomial.
        //
        {
            int i;
            int j;
            double[] v;

            v = new double[n];

            for ( j = 0; j < n; j++ )
            {
                v[j] = 1.0;
            }

            for ( i = 0; i < m; i++ )
            {
                if ( 0 != e[i] )
                {
                    for ( j = 0; j < n; j++ )
                    {
                        v[j] = v[j] * Math.Pow ( x[i+j*m], e[i] );
                    }
                }
            }

            return v;
        }
    }
}