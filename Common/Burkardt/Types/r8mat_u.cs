namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] r8mat_u_inverse(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_U_INVERSE inverts an upper triangular R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    An upper triangular matrix is a matrix whose only nonzero entries
        //    occur on or above the diagonal.
        //
        //    The inverse of an upper triangular matrix is an upper triangular matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 September 2005
        //
        //  Author:
        //
        //    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    Albert Nijenhuis, Herbert Wilf,
        //    Combinatorial Algorithms,
        //    Academic Press, 1978, second edition,
        //    ISBN 0-12-519260-6.
        //
        //  Parameters:
        //
        //    Input, int N, number of rows and columns in the matrix.
        //
        //    Input, double A[N*N], the upper triangular matrix.
        //
        //    Output, double R8MAT_U_INVERSE[N*N], the inverse matrix.
        //
    {
        double[] b;
        int i;
        int j;
        int k;

        b = new double[n * n];

        for (j = n - 1; 0 <= j; j--)
        {
            for (i = n - 1; 0 <= i; i--)
            {
                if (j < i)
                {
                    b[i + j * n] = 0.0;
                }
                else if (i == j)
                {
                    b[i + j * n] = 1.0 / a[i + j * n];
                }
                else
                {
                    b[i + j * n] = 0.0;
                    for (k = i + 1; k <= j; k++)
                    {
                        b[i + j * n] -= a[i + k * n] * b[k + j * n];
                    }

                    b[i + j * n] /= a[i + i * n];
                }
            }
        }

        return b;
    }

    public static double[] r8mat_u_solve(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_U_SOLVE solves an upper triangular linear system.
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
        //    21 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of
        //    the matrix A.
        //
        //    Input, double A[N*N], the N by N upper triangular matrix.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Output, double R8MAT_U_SOLVE[N], the solution of the linear system.
        //
    {
        int i;
        int j;
        double[] x;
        //
        //  Solve U * x = b.
        //
        x = new double[n];

        for (i = n - 1; 0 <= i; i--)
        {
            x[i] = b[i];
            for (j = i + 1; j < n; j++)
            {
                x[i] -= a[i + j * n] * x[j];
            }

            x[i] /= a[i + i * n];
        }

        return x;
    }

    public static double[] r8mat_u1_inverse(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_U1_INVERSE inverts a unit upper triangular R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    A unit upper triangular matrix is a matrix with only 1's on the main
        //    diagonal, and only 0's below the main diagonal.
        //
        //    The inverse of a unit upper triangular matrix is also
        //    a unit upper triangular matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 September 2005
        //
        //  Author:
        //
        //    C++ translation by John Burkardt
        //
        //  Reference:
        //
        //    Albert Nijenhuis, Herbert Wilf,
        //    Combinatorial Algorithms,
        //    Academic Press, 1978, second edition,
        //    ISBN 0-12-519260-6.
        //
        //  Parameters:
        //
        //    Input, int N, number of rows and columns in the matrix.
        //
        //    Input, double A[N*N], the unit upper triangular matrix.
        //
        //    Output, double R8MAT_U1_INVERSE[N*N), the inverse matrix.
        //
    {
        double[] b;
        int i;
        int j;
        int k;

        b = new double[n * n];

        for (j = n - 1; 0 <= j; j--)
        {
            for (i = n - 1; 0 <= i; i--)
            {
                if (j < i)
                {
                    b[i + j * n] = 0.0;
                }
                else if (i == j)
                {
                    b[i + j * n] = 1.0;
                }
                else
                {
                    b[i + j * n] = 0.0;
                    for (k = i + 1; k <= j; k++)
                    {
                        b[i + j * n] -= a[i + k * n] * b[k + j * n];
                    }

                    b[i + j * n] /= a[i + i * n];
                }
            }
        }

        return b;
    }
}