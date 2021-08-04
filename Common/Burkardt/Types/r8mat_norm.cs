using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double r8mat_norm_eis(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_NORM_EIS returns the EISPACK norm of an R8MAT.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //    The EISPACK norm is defined as:
            //
            //      R8MAT_NORM_EIS =
            //        sum ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows in A.
            //
            //    Input, int N, the number of columns in A.
            //
            //    Input, double A[M*N], the matrix whose EISPACK norm is desired.
            //
            //    Output, double R8MAT_NORM_EIS, the EISPACK norm of A.
            //
        {
            int i;
            int j;
            double value;

            value = 0.0;
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    value = value + Math.Abs(a[i + j * m]);
                }
            }

            return value;
        }


        public static double r8mat_norm_l1(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_NORM_L1 returns the matrix L1 norm of an R8MAT.
            //
            //  Discussion:
            //
            //    An R8MAT is an array of R8 values.
            //
            //    The matrix L1 norm is defined as:
            //
            //      R8MAT_NORM_L1 = max ( 1 <= J <= N )
            //        sum ( 1 <= I <= M ) abs ( A(I,J) ).
            //
            //    The matrix L1 norm is derived from the vector L1 norm, and
            //    satisifies:
            //
            //      r8vec_norm_l1 ( A * x ) <= r8mat_norm_l1 ( A ) * r8vec_norm_l1 ( x ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 December 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows in A.
            //
            //    Input, int N, the number of columns in A.
            //
            //    Input, double A(M,N), the matrix whose L1 norm is desired.
            //
            //    Output, double R8MAT_NORM_L1, the L1 norm of A.
            //
        {
            double col_sum;
            int i;
            int j;
            double value;

            value = 0.0;

            for (j = 0; j < n; j++)
            {
                col_sum = 0.0;
                for (i = 0; i < m; i++)
                {
                    col_sum = col_sum + Math.Abs(a[i + j * m]);
                }

                value = Math.Max(value, col_sum);
            }

            return value;
        }

        public static double r8mat_norm_l2(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_NORM_L2 returns the matrix L2 norm of an R8MAT.
            //
            //  Discussion:
            //
            //    An R8MAT is an array of R8 values.
            //
            //    The matrix L2 norm is defined as:
            //
            //      R8MAT_NORM_L2 = sqrt ( max ( 1 <= I <= M ) LAMBDA(I) )
            //
            //    where LAMBDA contains the eigenvalues of A * A'.
            //
            //    The matrix L2 norm is derived from the vector L2 norm, and
            //    satisifies:
            //
            //      r8vec_norm_l2 ( A * x ) <= r8mat_norm_l2 ( A ) * r8vec_norm_l2 ( x ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 December 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows in A.
            //
            //    Input, int N, the number of columns in A.
            //
            //    Input, double A(M,N), the matrix whose L2 norm is desired.
            //
            //    Output, double R8MAT_NORM_L2, the L2 norm of A.
            //
        {
            double[] at;
            double[] b;
            double[] diag;
            double value;

            at = r8mat_transpose_new(m, n, a);
            //
            //  Compute B = A * A'.
            //
            b = r8mat_mm_new(m, n, m, a, at);
            //
            //  Diagonalize B.
            //
            r8mat_symm_jacobi(m, ref b);
            //
            //  Find the maximum eigenvalue, and take its square root.
            //
            diag = r8mat_diag_get_vector_new(m, b);

            value = Math.Sqrt(r8vec_max(m, diag));

            return value;
        }

        public static double r8mat_norm_li(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_NORM_LI returns the matrix L-oo norm of an R8MAT.
            //
            //  Discussion:
            //
            //    An R8MAT is an array of R8 values.
            //
            //    The matrix L-oo norm is defined as:
            //
            //      R8MAT_NORM_LI =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).
            //
            //    The matrix L-oo norm is derived from the vector L-oo norm,
            //    and satisfies:
            //
            //      r8vec_norm_li ( A * x ) <= r8mat_norm_li ( A ) * r8vec_norm_li ( x ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 December 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows in A.
            //
            //    Input, int N, the number of columns in A.
            //
            //    Input, double A[M*N], the matrix whose L-oo
            //    norm is desired.
            //
            //    Output, double R8MAT_NORM_LI, the L-oo norm of A.
            //
        {
            int i;
            int j;
            double row_sum;
            double value;

            value = 0.0;

            for (i = 0; i < m; i++)
            {
                row_sum = 0.0;
                for (j = 0; j < n; j++)
                {
                    row_sum = row_sum + Math.Abs(a[i + j * m]);
                }

                value = Math.Max(value, row_sum);
            }

            return value;
        }

        public static double[] r8mat_normal_01_new(int m, int n, ref typeMethods.r8vecNormalData data, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_NORMAL_01_NEW returns a unit pseudonormal R8MAT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 April 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Paul Bratley, Bennett Fox, Linus Schrage,
            //    A Guide to Simulation,
            //    Springer Verlag, pages 201-202, 1983.
            //
            //    Bennett Fox,
            //    Algorithm 647:
            //    Implementation and Relative Efficiency of Quasirandom
            //    Sequence Generators,
            //    ACM Transactions on Mathematical Software,
            //    Volume 12, Number 4, pages 362-376, 1986.
            //
            //    Peter Lewis, Allen Goodman, James Miller,
            //    A Pseudo-Random Number Generator for the System/360,
            //    IBM Systems Journal,
            //    Volume 8, pages 136-143, 1969.
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the array.
            //
            //    Input/output, int &SEED, the "seed" value, which should NOT be 0.
            //    On output, SEED has been updated.
            //
            //    Output, double R8MAT_NORMAL_01_NEW[M*N], the array of pseudonormal values.
            //
        {
            double[] r;

            r = r8vec_normal_01_new(m * n, ref data, ref seed);

            return r;
        }

        public static double r8mat_norm_fro(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_NORM_FRO returns the Frobenius norm of an R8MAT.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //    The Frobenius norm is defined as
            //
            //      R8MAT_NORM_FRO = sqrt (
            //        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)**2 )
            //    The matrix Frobenius norm is not derived from a vector norm, but
            //    is compatible with the vector L2 norm, so that:
            //
            //      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows in A.
            //
            //    Input, int N, the number of columns in A.
            //
            //    Input, double A[M*N], the matrix whose Frobenius
            //    norm is desired.
            //
            //    Output, double R8MAT_NORM_FRO, the Frobenius norm of A.
            //
        {
            int i;
            int j;
            double value;

            value = 0.0;
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    value = value + Math.Pow(a[i + j * m], 2);
                }
            }

            value = Math.Sqrt(value);

            return value;
        }

        public static double r8mat_norm_fro_affine(int m, int n, double[] a1, double[] a2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_NORM_FRO_AFFINE returns the Frobenius norm of an R8MAT difference.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //    The Frobenius norm is defined as
            //
            //      R8MAT_NORM_FRO = sqrt (
            //        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)^2 )
            //    The matrix Frobenius norm is not derived from a vector norm, but
            //    is compatible with the vector L2 norm, so that:
            //
            //      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows.
            //
            //    Input, int N, the number of columns.
            //
            //    Input, double A1[M*N], A2[M,N], the matrice for whose difference the 
            //    Frobenius norm is desired.
            //
            //    Output, double R8MAT_NORM_FRO_AFFINE, the Frobenius norm of A1 - A2.
            //
        {
            double value = 0.0;
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < m; i++)
                {
                    value = value + Math.Pow(a1[i + j * m] - a2[i + j * m], 2);
                }
            }

            value = Math.Sqrt(value);

            return value;
        }
        
    }
}