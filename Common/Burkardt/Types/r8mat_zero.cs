namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] r8mat_zero_new(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_ZERO_NEW returns a new zeroed R8MAT.
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
        //    03 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Output, double R8MAT_ZERO[M*N], the new zeroed matrix.
        //
    {
        double[] a = new double[m * n];

        for (int j = 0; j < n; j++)
        {
            for (int i = 0; i < m; i++)
            {
                a[i + j * m] = 0.0;
            }
        }

        return a;
    }

    public static void r8mat_zeros ( int m, int n, ref double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_ZEROS zeroes an R8MAT.
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
        //    Output, double A[M*N], a matrix of zeroes.
        //
    {
        int j;

        for ( j = 0; j < n; j++ )
        {
            int i;
            for ( i = 0; i < m; i++ )
            {
                a[i+j*m] = 0.0;
            }
        }
    }

    public static double[] r8mat_zeros_new(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_ZEROS_NEW returns a new zeroed R8MAT.
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
        //    03 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Output, double R8MAT_ZEROS_NEW[M*N], the new zeroed matrix.
        //
    {
        int j;

        double[] a = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                a[i + j * m] = 0.0;
            }
        }

        return a;
    }
        
}