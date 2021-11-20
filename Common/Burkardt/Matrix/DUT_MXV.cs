﻿namespace Burkardt.MatrixNS;

public static partial class Matrix
{
    public static double[] dut_mxv ( int m, int n, double[] a, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DUT_MXV multiplies an DUT matrix times a vector.
        //
        //  Discussion:
        //
        //    The DUT storage format is used for an M by N upper triangular matrix,
        //    and allocates space even for the zero entries.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 September 2003
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
        //    Input, double A[M*N], the DUT matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double DUT_MXV[M], the product A * x.
        //
    {
        int i;

        double[] b = new double[m];

        for ( i = 0; i < m; i++ )
        {
            b[i] = 0.0;
            int j;
            for ( j = i; j < n; j++ )
            {
                b[i] += a[i+j*m] * x[j];
            }
        }

        return b;
    }
}