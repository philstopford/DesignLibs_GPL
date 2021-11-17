namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] r8vec_legendre_new ( int n, double a_first, double a_last )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_LEGENDRE_NEW creates a vector of Chebyshev spaced values.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 June 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double A_FIRST, A_LAST, the first and last entries.
        //
        //    Output, double R8VEC_LEGENDRE_NEW[N], a vector of Legendre spaced data.
        //
    {
        int i;

        double[] a = PolynomialNS.Legendre.legendre_zeros ( n );

        for ( i = 0; i < n; i++ )
        {
            a[i] = ( ( 1.0 - a[i] ) * a_first  
                     + ( 1.0 + a[i] ) * a_last ) 
                   /   2.0;
        }
        return a;
    }
}