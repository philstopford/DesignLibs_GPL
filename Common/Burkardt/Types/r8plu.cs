namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8plu_det(int n, int[] pivot, double[] lu )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PLU_DET computes the determinant of a real PLU matrix.
        //
        //  Discussion:
        //
        //    The matrix should have been factored by R8MAT_TO_R8PLU.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 August 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, 1979
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int PIVOT[N], the pivot vector computed by R8MAT_TO_R8PLU.
        //
        //    Input, double LU[N*N], the LU factors computed by R8MAT_TO_R8PLU.
        //
        //    Output, double R8PLU_DET, the determinant of the matrix.
        //
    {
        double det;
        int i;

        det = 1.0;

        for (i = 0; i < n; i++)
        {
            det *= lu[i + i * n];
            if (pivot[i] != i + 1)
            {
                det = -det;
            }
        }

        return det;
    }

    public static void r8plu_inverse(int n, int[] pivot, double[] lu, ref double[] a_inverse )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PLU_INVERSE computes the inverse of a real PLU matrix.
        //
        //  Discussion:
        //
        //    The matrix should have been factored by R8MAT_TO_R8PLU.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 August 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix A.
        //
        //    Input, int PIVOT[N], the pivot vector from R8MAT_TO_R8PLU.
        //
        //    Input, double LU[N*N], the LU factors computed by R8MAT_TO_R8PLU.
        //
        //    Output, double A_INVERSE[N*N], the inverse of the original matrix
        //    A that was factored by R8MAT_TO_R8PLU.
        //
    {
        int i;
        int j;
        int k;
        double temp;
        double[] work;
        //
        work = new double[n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                a_inverse[i + j * n] = lu[i + j * n];
            }
        }

        //
        //  Compute Inverse(U).
        //
        for (k = 1; k <= n; k++)
        {
            a_inverse[k - 1 + (k - 1) * n] = 1.0 / a_inverse[k - 1 + (k - 1) * n];
            for (i = 1; i <= k - 1; i++)
            {
                a_inverse[i - 1 + (k - 1) * n] = -a_inverse[i - 1 + (k - 1) * n] * a_inverse[k - 1 + (k - 1) * n];
            }

            for (j = k + 1; j <= n; j++)
            {
                temp = a_inverse[k - 1 + (j - 1) * n];
                a_inverse[k - 1 + (j - 1) * n] = 0.0;
                for (i = 1; i <= k; i++)
                {
                    a_inverse[i - 1 + (j - 1) * n] += temp * a_inverse[i - 1 + (k - 1) * n];
                }
            }
        }

        //
        //  Form Inverse(U) * Inverse(L).
        //
        for (k = n - 1; 1 <= k; k--)
        {
            for (i = k + 1; i <= n; i++)
            {
                work[i - 1] = a_inverse[i - 1 + (k - 1) * n];
                a_inverse[i - 1 + (k - 1) * n] = 0.0;
            }

            for (j = k + 1; j <= n; j++)
            {
                for (i = 1; i <= n; i++)
                {
                    a_inverse[i - 1 + (k - 1) * n] += a_inverse[i - 1 + (j - 1) * n] * work[j - 1];
                }
            }

            if (pivot[k - 1] != k)
            {
                for (i = 1; i <= n; i++)
                {
                    temp = a_inverse[i - 1 + (k - 1) * n];
                    a_inverse[i - 1 + (k - 1) * n] = a_inverse[i - 1 + (pivot[k - 1] - 1) * n];
                    a_inverse[i - 1 + (pivot[k - 1] - 1) * n] = temp;
                }
            }
        }
    }

    public static void r8plu_mul(int n, int[] pivot, double[] lu, double[] x, ref double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PLU_MUL computes A * x using the PLU factors of A.
        //
        //  Discussion:
        //
        //    It is assumed that R8MAT_TO_R8PLU has computed the PLU factors of
        //    the matrix A.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 August 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int PIVOT[N], the pivot vector computed by R8MAT_TO_R8PLU.
        //
        //    Input, double LU[N*N], the matrix factors computed by R8MAT_TO_R8PLU.
        //
        //    Input, double X[N], the vector to be multiplied.
        //
        //    Output, double B[N], the result of the multiplication.
        //
    {
        int i;
        int j;
        int k;
        double temp;
        //
        for (i = 0; i < n; i++)
        {
            b[i] = x[i];
        }

        //
        //  Y = U * X.
        //
        for (j = 1; j <= n; j++)
        {
            for (i = 0; i < j - 1; i++)
            {
                b[i] += lu[i + (j - 1) * n] * b[j - 1];
            }

            b[j - 1] = lu[j - 1 + (j - 1) * n] * b[j - 1];
        }

        //
        //  B = PL * Y = PL * U * X = A * x.
        //
        for (j = n - 1; 1 <= j; j--)
        {
            for (i = j; i < n; i++)
            {
                b[i] -= lu[i + (j - 1) * n] * b[j - 1];
            }

            k = pivot[j - 1];

            if (k != j)
            {
                temp = b[k - 1];
                b[k - 1] = b[j - 1];
                b[j - 1] = temp;
            }
        }
    }

    public static void r8plu_sol(int n, int[] pivot, double[] lu, double[] b, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PLU_SOL solves a linear system A*x=b from the PLU factors.
        //
        //  Discussion:
        //
        //    The PLU factors should have been computed by R8MAT_TO_R8PLU.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 August 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int PIVOT[N], the pivot vector from R8MAT_TO_R8PLU.
        //
        //    Input, double LU[N*N], the LU factors from R8MAT_TO_R8PLU.
        //
        //    Input, double B[N], the right hand side vector.
        //
        //    Output, double X[N], the solution vector.
        //
    {
        int i;
        int j;
        int k;
        double temp;
        //
        //  Solve PL * Y = B.
        //
        for (i = 0; i < n; i++)
        {
            x[i] = b[i];
        }

        for (k = 1; k <= n - 1; k++)
        {
            j = pivot[k - 1];

            if (j != k)
            {
                temp = x[j - 1];
                x[j - 1] = x[k - 1];
                x[k - 1] = temp;
            }

            for (i = k + 1; i <= n; i++)
            {
                x[i - 1] += lu[i - 1 + (k - 1) * n] * x[k - 1];
            }
        }

        //
        //  Solve U * X = Y.
        //
        for (k = n; 1 <= k; k--)
        {
            x[k - 1] /= lu[k - 1 + (k - 1) * n];
            for (i = 1; i <= k - 1; i++)
            {
                x[i - 1] -= lu[i - 1 + (k - 1) * n] * x[k - 1];
            }
        }
    }

    public static void r8plu_to_r8mat(int n, int[] pivot, double[] lu, ref double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PLU_TO_R8MAT recovers the matrix A that was factored by R8MAT_TO_R8PLU.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 August 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int PIVOT[N], the pivot vector computed by R8MAT_TO_R8PLU.
        //
        //    Input, double LU[N*N], the matrix factors computed by R8MAT_TO_R8PLU.
        //
        //    Output, double A[N*N], the matrix whose factors are represented by
        //    LU and PIVOT.
        //
    {
        int i;
        int j;
        int k;
        double temp;

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                if (i == j)
                {
                    a[i + j * n] = 1.0;
                }
                else
                {
                    a[i + j * n] = 0.0;
                }
            }
        }

        for (j = 1; j <= n; j++)
        {
            for (i = 1; i <= n; i++)
            {
                for (k = 1; k <= i - 1; k++)
                {
                    a[k - 1 + (j - 1) * n] += lu[k - 1 + (i - 1) * n] * a[i - 1 + (j - 1) * n];
                }

                a[i - 1 + (j - 1) * n] = lu[i - 1 + (i - 1) * n] * a[i - 1 + (j - 1) * n];
            }

            //
            //  B = PL * Y = PL * U * X = A * x.
            //
            for (i = n - 1; 1 <= i; i--)
            {
                for (k = i + 1; k <= n; k++)
                {
                    a[k - 1 + (j - 1) * n] -= lu[k - 1 + (i - 1) * n] * a[i - 1 + (j - 1) * n];
                }

                k = pivot[i - 1];

                if (k != i)
                {
                    temp = a[k - 1 + (j - 1) * n];
                    a[k - 1 + (j - 1) * n] = a[i - 1 + (j - 1) * n];
                    a[i - 1 + (j - 1) * n] = temp;
                }
            }
        }
    }
}