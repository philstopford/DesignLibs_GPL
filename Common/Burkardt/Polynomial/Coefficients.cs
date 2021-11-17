using System;

namespace Burkardt.PolynomialNS;

public static class Coefficients
{
    public static int[] binomial_table ( int qs, int m, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BINOMIAL_TABLE computes a table of bionomial coefficients MOD QS.
        //
        //  Discussion:
        //
        //    Thanks to Michael Baudin for pointing out an error in a previous
        //    version of this function, 07 December 2009.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 December 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int QS, the base for the MOD operation.
        //
        //    Input, int M, N, the limits of the binomial table.
        //
        //    Output, int BINOMIAL_TABLE[(M+1)*(N+1)], the table of binomial 
        //    coefficients modulo QS.
        //
    {
        int[] coef;
        int i;
        int j;

        coef = new int[(m+1)*(n+1)];

        for ( j = 0; j <= n; j++ )
        {
            for ( i = 0; i <= m; i++ )
            {
                coef[i+j*(m+1)] = 0;
            }
        }

        coef[0] = 1;

        j = 0;
        for ( i = 1; i <= m; i++ )
        {
            coef[i+j*(m+1)] = 1;
        }

        for ( i = 1; i <= Math.Min ( m, n ); i++ )
        {
            j = i;
            coef[i+j*(m+1)] = 1;
        }

        for( j = 1; j <= n; j++ )
        {
            for ( i = j + 1; i <= m; i++ )
            {
                coef[i+j*(m+1)] = ( coef[i-1+j*(m+1)] + coef[i-1+(j-1)*(m+1)] ) % qs;
            }
        }

        return coef;
    }
}