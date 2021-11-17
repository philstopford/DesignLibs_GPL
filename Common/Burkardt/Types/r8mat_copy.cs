namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8mat_copy(int m, int n, double[] a1, ref double[] a2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_COPY copies one R8MAT to another.
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
        //    16 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A1[M*N], the matrix to be copied.
        //
        //    Output, double A2[M*N], the copy of A1.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                a2[i + j * m] = a1[i + j * m];
            }
        }
    }

    public static double[] r8mat_copy_new(int m, int n, double[] a1)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_COPY_NEW copies one R8MAT to a "new" R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8's, which
        //    may be stored as a vector in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A1[M*N], the matrix to be copied.
        //
        //    Output, double R8MAT_COPY_NEW[M*N], the copy of A1.
        //
    {
        int j;

        double[] a2 = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                a2[i + j * m] = a1[i + j * m];
            }
        }

        return a2;
    }
 
    public static void r8mat_row_copy(int m, int n, int i, double[] v, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_ROW_COPY copies a vector into a row of an R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, int I, the index of the row.
        //    0 <= I <= M-1.
        //
        //    Input, double V[N], the row to be copied.
        //
        //    Input/output, double A[M*N], the matrix into which
        //    the row is to be copied.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {
            a[i + j * m] = v[j];
        }
    }

        
}