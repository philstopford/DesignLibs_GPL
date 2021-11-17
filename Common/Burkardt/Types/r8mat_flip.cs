namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] r8mat_flip_cols_new(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_FLIP_COLS_NEW makes a new copy of an R8MAT with reversed column order.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the matrix to be copied.
        //
        //    Output, double R8MAT_FLIP_COLS_NEW[M*N], the reversed-column-order copy.
        //
    {
        double[] b;
        int i;
        int j;

        b = new double[m * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                b[i + (n - 1 - j) * m] = a[i + j * m];
            }
        }

        return b;
    }

    public static double[] r8mat_flip_rows_new(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_FLIP_ROWS_NEW makes a new copy of an R8MAT with reversed row order.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the matrix to be copied.
        //
        //    Output, double R8MAT_FLIP_ROWS_NEW[M*N], the reversed-rows-order copy.
        //
    {
        double[] b;
        int i;
        int j;

        b = new double[m * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                b[m - 1 - i + j * m] = a[i + j * m];
            }
        }

        return b;
    }
        
}