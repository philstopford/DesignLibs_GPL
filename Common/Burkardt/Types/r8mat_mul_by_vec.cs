namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8mat_mv(int m, int n, double[] a, double[] x, ref double[] ax)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_MV multiplies a matrix times a vector.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    For this routine, the result is returned as an argument.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Input, double A[M,N], the M by N matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double AX[M], the product A*X.
        //
    {
        int i;

        double[] y = new double[m];

        for (i = 0; i < m; i++)
        {
            y[i] = 0.0;
            int j;
            for (j = 0; j < n; j++)
            {
                y[i] += a[i + j * m] * x[j];
            }
        }

        ax = new double[y.Length];

        r8vec_copy(m, y, ref ax);

    }

    public static double[] r8mat_mv_new(int m, int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_MV_NEW multiplies a matrix times a vector.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    For this routine, the result is returned as the function value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Input, double A[M,N], the M by N matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8MAT_MV_NEW[M], the product A*X.
        //
    {
        int i;

        double[] y = new double[m];

        for (i = 0; i < m; i++)
        {
            y[i] = 0.0;
            int j;
            for (j = 0; j < n; j++)
            {
                y[i] += a[i + j * m] * x[j];
            }
        }

        return y;
    }

    public static void r8mat_mtv ( int m, int n, double[] a, double[] x, ref double[] atx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_MTV multiplies a transposed matrix times a vector.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    For this routine, the result is returned as an argument.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Input, double A[M,N], the M by N matrix.
        //
        //    Input, double X[M], the vector to be multiplied by A.
        //
        //    Output, double ATX[N], the product A'*X.
        //
    {
        int j;

        double[] y = new double[n];

        for ( j = 0; j < n; j++ )
        {
            y[j] = 0.0;
            int i;
            for ( i = 0; i < m; i++ )
            {
                y[j] += a[i+j*m] * x[i];
            }
        }

        r8vec_copy ( n, y, ref atx );

    }

    public static double[] r8mat_mtv_new(int m, int n, double[] a, double[] x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_MTV_NEW multiplies a transposed matrix times a vector.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    For this routine, the result is returned as the function value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Input, double A[M,N], the M by N matrix.
        //
        //    Input, double X[M], the vector to be multiplied by A.
        //
        //    Output, double R8MAT_MTV_NEW[N], the product A'*X.
        //
    {
        double[] y = new double[n];

        for (int j = 0; j < n; j++)
        {
            y[j] = 0.0;
            for (int i = 0; i < m; i++)
            {
                y[j] += a[i + j * m] * x[i];
            }
        }

        return y;
    }
        
}