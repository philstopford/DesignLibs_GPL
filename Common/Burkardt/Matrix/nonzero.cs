namespace Burkardt.MatrixNS;

public static partial class Matrix
{
    public static int nonzeros ( int m, int n, double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NONZEROS counts the nonzeros in a matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the matrix.
        //
        //    Output, int NONZEROS, the number of nonzero entries.
        //
    {
        int j;

        int nnz = 0;
        for ( j = 0; j < n; j++ )
        {
            int i;
            for ( i = 0; i < m; i++ )
            {
                if ( a[i+j*m] != 0.0 )
                {
                    nnz += 1;
                }
            }
        }

        return nnz;
    }
}