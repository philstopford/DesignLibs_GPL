namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8mat_diag_add_scalar(int n, ref double[] a, double s)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_DIAG_ADD_SCALAR adds a scalar to the diagonal of an R8MAT.
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
        //    07 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //
        //    Input/output, double A[N*N], the N by N matrix to be modified.
        //
        //    Input, double S, the value to be added to the diagonal
        //    of the matrix.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            a[i + i * n] += s;
        }
    }

    public static void r8mat_diag_add_vector(int n, ref double[] a, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_DIAG_ADD_VECTOR adds a vector to the diagonal of an R8MAT.
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
        //    07 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //
        //    Input/output, double A[N*N], the N by N matrix.
        //
        //    Input, double V[N], the vector to be added to the diagonal of A.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            a[i + i * n] += v[i];
        }

    }

    public static void r8mat_diag_get_vector(int n, double[] a, ref double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_DIAG_GET_VECTOR gets the value of the diagonal of an R8MAT.
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
        //    15 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //
        //    Input, double A[N*N], the N by N matrix.
        //
        //    Output, double V[N], the diagonal entries
        //    of the matrix.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            v[i] = a[i + i * n];
        }

    }

    public static double[] r8mat_diag_get_vector_new(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_DIAG_GET_VECTOR_NEW gets the value of the diagonal of an R8MAT.
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
        //    15 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //
        //    Input, double A[N*N], the N by N matrix.
        //
        //    Output, double R8MAT_DIAG_GET_VECTOR_NEW[N], the diagonal entries
        //    of the matrix.
        //
    {
        int i;

        double[] v = new double[n];

        for (i = 0; i < n; i++)
        {
            v[i] = a[i + i * n];
        }

        return v;
    }

    public static void r8mat_diag_set_scalar(int n, ref double[] a, double s)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_DIAG_SET_SCALAR sets the diagonal of an R8MAT to a scalar value.
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
        //    07 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //
        //    Input/output, double A[N*N], the N by N matrix to be modified.
        //
        //    Input, double S, the value to be assigned to the diagonal
        //    of the matrix.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            a[i + i * n] = s;
        }
    }

    public static void r8mat_diag_set_vector(int n, ref double[] a, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_DIAG_SET_VECTOR sets the diagonal of an R8MAT to a vector.
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
        //    07 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //
        //    Input/output, double A[N*N], the N by N matrix.
        //
        //    Input, double V[N], the vector to be assigned to the
        //    diagonal of A.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            a[i + i * n] = v[i];
        }
    }

    public static double[] r8mat_diagonal_new(int n, double[] diag)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_DIAGONAL_NEW returns a diagonal matrix.
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
        //    31 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of A.
        //
        //    Input, double DIAG[N], the diagonal entries.
        //
        //    Output, double R8MAT_DIAGONAL_NEW[N*N], the N by N identity matrix.
        //
    {
        int j;

        double[] a = new double[n * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                if (i == j)
                {
                    a[i + j * n] = diag[i];
                }
                else
                {
                    a[i + j * n] = 0.0;
                }
            }
        }

        return a;
    }
        
}