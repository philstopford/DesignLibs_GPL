namespace Burkardt
{
    public static partial class Matrix
    {
        public static double[] dge_mxv ( int m, int n, double[] a, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGE_MXV multiplies a DGE matrix times a vector.
        //
        //  Discussion:
        //
        //    The DGE storage format is used for a general M by N matrix.  A physical storage
        //    space is made for each logical entry.  The two dimensional logical
        //    array is mapped to a vector, in which storage is by columns.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Input, double A[M*N], the SGE matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double DGE_MXV[M], the product A * x.
        //
        {
            double[] b;
            int i;
            int j;

            b = new double[m];

            for ( i = 0; i < m; i++ )
            {
                b[i] = 0.0;
                for ( j = 0; j < n; j++ )
                {
                    b[i] = b[i] + a[i+j*m] * x[j];
                }
            }

            return b;
        }
    }
}