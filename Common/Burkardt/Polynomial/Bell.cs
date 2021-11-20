namespace Burkardt.PolynomialNS;

public static class Bell
{
    public static int[] bell_poly_coef ( int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BELL_POLY_COEF: Coefficients of a Bell polynomial.
        //
        //  First terms:
        //
        //    N    0    1    2    3    4    5    6    7    8
        //
        //    0    1
        //    1    0    1    
        //    2    0    1    1
        //    3    0    1    3    1
        //    4    0    1    7    6    1
        //    5    0    1   15   25   10    1
        //    6    0    1   31   90   65   15    1
        //    7    0    1   63  301  350  140   21    1
        //    8    0    1  127  966 1701 1050  266   28    1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 March 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the polynomial.
        //
        //    Output, int BELL_POLY_COEF[N+1], the coefficients.
        //
    {
        int i;
        int j;

        int[] c = new int[n+1];

        c[0] = 1;
        for ( j = 1; j <= n; j++ )
        {
            c[j] = 0;
        }
 
        for ( i = 1; i <= n; i++ )
        {
            for ( j = i; 1 <= j; j-- )
            {
                c[j] = j * c[j] + c[j-1];
            }
            c[0] = 0;
        }
 
        return c;
    }
}