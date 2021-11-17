using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] r8mat_cholesky_factor(int n, double[] a, ref int flag)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_CHOLESKY_FACTOR computes the Cholesky factor of a symmetric R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    The matrix must be symmetric and positive semidefinite.
        //
        //    For a positive semidefinite symmetric matrix A, the Cholesky factorization
        //    is a lower triangular matrix L such that:
        //
        //      A = L * L'
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix A.
        //
        //    Input, double A[N*N], the N by N matrix.
        //
        //    Output, int &FLAG, an error flag.
        //    0, no error occurred.
        //    1, the matrix is not positive definite.
        //    2, the matrix is not nonnegative definite.
        //
        //    Output, double R8MAT_CHOLESKY_FACTOR[N*N], the N by N lower triangular
        //    Cholesky factor.
        //
    {
        int j;

        flag = 0;
        double tol = Math.Sqrt(r8_epsilon());

        double[] c = r8mat_copy_new(n, n, a);

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < j; i++)
            {
                c[i + j * n] = 0.0;
            }

            for (i = j; i < n; i++)
            {
                double sum2 = c[j + i * n];
                int k;
                for (k = 0; k < j; k++)
                {
                    sum2 -= c[j + k * n] * c[i + k * n];
                }

                if (i == j)
                {
                    switch (sum2)
                    {
                        case > 0.0:
                            c[i + j * n] = Math.Sqrt(sum2);
                            break;
                        default:
                        {
                            if (sum2 < -tol)
                            {
                                flag = 2;
                                Console.WriteLine("");
                                Console.WriteLine("R8MAT_CHOLESKY_FACTOR - Fatal error!");
                                Console.WriteLine("  Matrix is not nonnegative definite.");
                                Console.WriteLine("  Diagonal I = " + i + "");
                                Console.WriteLine("  SUM2 = " + sum2 + "");
                                return null;
                            }
                            else
                            {
                                flag = 1;
                                c[i + j * n] = 0.0;
                            }

                            break;
                        }
                    }
                }
                else
                {

                    if (c[j + j * n] != 0.0)
                    {
                        c[i + j * n] = sum2 / c[j + j * n];
                    }
                    else
                    {
                        c[i + j * n] = 0.0;
                    }
                }
            }
        }

        return c;
    }

    public static double[] r8mat_cholesky_factor_upper(int n, double[] a, ref int flag)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_CHOLESKY_FACTOR_UPPER: upper Cholesky factor of a symmetric R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    The matrix must be symmetric and positive semidefinite.
        //
        //    For a positive semidefinite symmetric matrix A, the Cholesky factorization
        //    is an upper triangular matrix R such that:
        //
        //      A = R' * R
        //
        //    Note that the usual Cholesky factor is a LOWER triangular matrix L
        //    such that
        //
        //      A = L * L'
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix A.
        //
        //    Input, double A[N*N], the N by N matrix.
        //
        //    Output, int &FLAG, an error flag.
        //    0, no error occurred.
        //    1, the matrix is not positive definite.  A NULL factor is returned.
        //
        //    Output, double R8MAT_CHOLESKY_FACTOR[N*N], the N by N upper triangular
        //    Cholesky factor.
        //
    {
        int j;

        flag = 0;

        double[] c = r8mat_copy_new(n, n, a);

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < j; i++)
            {
                c[j + i * n] = 0.0;
            }

            for (i = j; i < n; i++)
            {
                double sum2 = c[i + j * n];
                int k;
                for (k = 0; k < j; k++)
                {
                    sum2 -= c[k + j * n] * c[k + i * n];
                }

                if (i == j)
                {
                    switch (sum2)
                    {
                        case <= 0.0:
                            flag = 1;
                            return null;
                        default:
                            c[j + i * n] = Math.Sqrt(sum2);
                            break;
                    }
                }
                else
                {
                    if (c[j + j * n] != 0.0)
                    {
                        c[j + i * n] = sum2 / c[j + j * n];
                    }
                    else
                    {
                        c[j + i * n] = 0.0;
                    }
                }
            }
        }

        return c;
    }

    public static void r8mat_cholesky_inverse(int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_CHOLESKY_INVERSE computes the inverse of a symmetric matrix.
        //
        //  Discussion:
        //
        //    The matrix must be symmetric and positive semidefinite.
        //
        //    The upper triangular Cholesky factorization R is computed, so that:
        //
        //      A = R' * R
        //
        //    Then the inverse B is computed by
        //
        //      B = inv ( A ) = inv ( R ) * inv ( R' )
        //
        //    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2013
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
        //    Input/output, double A[N*N].  On input, the matrix.
        //    On output, the inverse of the matrix.
        //
    {
        int i;
        int j;
        int k;
        double t;

        for (j = 0; j < n; j++)
        {
            double s = 0.0;

            for (k = 0; k < j; k++)
            {
                t = a[k + j * n];
                for (i = 0; i < k; i++)
                {
                    t -= a[i + k * n] * a[i + j * n];
                }

                t /= a[k + k * n];
                a[k + j * n] = t;
                s += t * t;
            }

            s = a[j + j * n] - s;

            switch (s)
            {
                case <= 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("R8MAT_CHOLESKY_INVERSE - Fatal error!");
                    Console.WriteLine("  The matrix is singular.");
                    return;
            }

            a[j + j * n] = Math.Sqrt(s);

            for (i = j + 1; i < n; i++)
            {
                a[i + j * n] = 0.0;
            }
        }

        //
        //  Compute inverse(R).
        //
        for (k = 0; k < n; k++)
        {
            a[k + k * n] = 1.0 / a[k + k * n];
            for (i = 0; i < k; i++)
            {
                a[i + k * n] = -a[i + k * n] * a[k + k * n];
            }

            for (j = k + 1; j < n; j++)
            {
                t = a[k + j * n];
                a[k + j * n] = 0.0;
                for (i = 0; i <= k; i++)
                {
                    a[i + j * n] += t * a[i + k * n];
                }
            }
        }

        //
        //  Form inverse(R) * (inverse(R))'.
        //
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < j; k++)
            {
                t = a[k + j * n];
                for (i = 0; i <= k; i++)
                {
                    a[i + k * n] += t * a[i + j * n];
                }
            }

            t = a[j + j * n];
            for (i = 0; i <= j; i++)
            {
                a[i + j * n] *= t;
            }
        }

        //
        //  Use reflection.
        //
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < i; j++)
            {
                a[i + j * n] = a[j + i * n];
            }
        }
    }

    public static double[] r8mat_cholesky_solve(int n, double[] l, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_CHOLESKY_SOLVE solves a Cholesky factored linear system A * x = b.
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
        //    Input, int N, the number of rows and columns of the matrix A.
        //
        //    Input, double L[N*N], the N by N Cholesky factor of the
        //    system matrix A.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Output, double R8MAT_CHOLESKY_SOLVE[N], the solution of the linear system.
        //
    {
        //
        //  Solve L * y = b.
        //
        double[] y = r8mat_l_solve(n, l, b);
        //
        //  Solve L' * x = y.
        //
        double[] x = r8mat_lt_solve(n, l, y);

        return x;
    }

    public static double[] r8mat_cholesky_solve_upper(int n, double[] r, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_CHOLESKY_SOLVE_UPPER solves Cholesky factored linear system A * x = b.
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
        //    21 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix A.
        //
        //    Input, double R[N*N], the N by N Cholesky factor of the
        //    system matrix A.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Output, double R8MAT_CHOLESKY_SOLVE_UPPER[N], the solution of the linear system.
        //
    {
        //
        //  Solve U' * y = b.
        //
        double[] y = r8mat_ut_solve(n, r, b);
        //
        //  Solve U * x = y.
        //
        double[] x = r8mat_u_solve(n, r, y);

        return x;
    }
        
}