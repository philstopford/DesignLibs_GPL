using System;
using System.Numerics;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void c8mat_add(int m, int n, Complex alpha, Complex[] a,
            Complex beta, Complex[] b, Complex[] c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_ADD combines two C8MAT's using complex scalar factors.
        //
        //  Discussion:
        //
        //    An C8MAT is a doubly dimensioned array of complex double precision values, 
        //    which may be stored as a vector in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 March 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, Complex ALPHA, the first scale factor.
        //
        //    Input, Complex A[M*N], the first matrix.
        //
        //    Input, Complex BETA, the second scale factor.
        //
        //    Input, Complex B[M*N], the second matrix.
        //
        //    Output, Complex C[M*N], the result.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                c[i + j * m] = alpha * a[i + j * m] + beta * b[i + j * m];
            }
        }
    }

    public static void c8mat_add_r8(int m, int n, double alpha, Complex[] a,
            double beta, Complex[] b, ref Complex[] c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_ADD_R8 combines two C8MAT's using real scalar factors.
        //
        //  Discussion:
        //
        //    An C8MAT is a doubly dimensioned array of complex double precision values, 
        //    which may be stored as a vector in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 March 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double ALPHA, the first scale factor.
        //
        //    Input, Complex A[M*N], the first matrix.
        //
        //    Input, double BETA, the second scale factor.
        //
        //    Input, Complex B[M*N], the second matrix.
        //
        //    Output, Complex C[M*N], the result.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                c[i + j * m] = alpha * a[i + j * m] + beta * b[i + j * m];
            }
        }
    }

    public static void c8mat_copy(int m, int n, Complex[] a1, ref Complex[] a2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_COPY copies one C8MAT to another.
        //
        //  Discussion:
        //
        //    An C8MAT is a doubly dimensioned array of complex double precision values, 
        //    which may be stored as a vector in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, Complex A1[M*N], the matrix to be copied.
        //
        //    Output, Complex A2[M*N], the copy of A1.
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

    public static Complex[] c8mat_copy_new(int m, int n, Complex[] a1)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_COPY_NEW copies one C8MAT to a "new" C8MAT.
        //
        //  Discussion:
        //
        //    An C8MAT is a doubly dimensioned array of complex double precision values, 
        //    which may be stored as a vector in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, Complex A1[M*N], the matrix to be copied.
        //
        //    Output, Complex C8MAT_COPY_NEW[M*N], the copy of A1.
        //
    {
        int j;

        Complex[] a2 = new Complex[m * n];

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

    public static void c8mat_fss(int n, ref Complex[] a, int nb, ref Complex[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_FSS factors and solves a system with multiple right hand sides.
        //
        //  Discussion:
        //
        //    This routine uses partial pivoting, but no pivot vector is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 March 2013
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
        //    Input/output, Complex A[N*N].
        //    On input, A is the coefficient matrix of the linear system.
        //    On output, A is in unit upper triangular form, and
        //    represents the U factor of an LU factorization of the
        //    original coefficient matrix.
        //
        //    Input, int NB, the number of right hand sides.
        //
        //    Input/output, Complex X[N*NB], on input, the right hand sides of the
        //    linear systems.  On output, the solutions of the linear systems.
        //
    {
        int i;
        int j;
        int jcol;

        for (jcol = 1; jcol <= n; jcol++)
        {
            //
            //  Find the maximum element in column I.
            //
            double piv = c8_abs(a[jcol - 1 + (jcol - 1) * n]);
            int ipiv = jcol;
            for (i = jcol + 1; i <= n; i++)
            {
                if (piv < c8_abs(a[i - 1 + (jcol - 1) * n]))
                {
                    piv = c8_abs(a[i - 1 + (jcol - 1) * n]);
                    ipiv = i;
                }
            }

            switch (piv)
            {
                case 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("C8MAT_FSS - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + jcol + "");
                    return;
            }

            //
            //  Switch rows JCOL and IPIV, and X.
            //
            Complex t;
            if (jcol != ipiv)
            {
                for (j = 1; j <= n; j++)
                {
                    t = a[jcol - 1 + (j - 1) * n];
                    a[jcol - 1 + (j - 1) * n] = a[ipiv - 1 + (j - 1) * n];
                    a[ipiv - 1 + (j - 1) * n] = t;
                }

                for (j = 0; j < nb; j++)
                {
                    t = x[jcol - 1 + j * n];
                    x[jcol - 1 + j * n] = x[ipiv - 1 + j * n];
                    x[ipiv - 1 + j * n] = t;
                }
            }

            //
            //  Scale the pivot row.
            //
            t = a[jcol - 1 + (jcol - 1) * n];
            a[jcol - 1 + (jcol - 1) * n] = 1.0;
            for (j = jcol + 1; j <= n; j++)
            {
                a[jcol - 1 + (j - 1) * n] /= t;
            }

            for (j = 0; j < nb; j++)
            {
                x[jcol - 1 + j * n] /= t;
            }

            //
            //  Use the pivot row to eliminate lower entries in that column.
            //
            for (i = jcol + 1; i <= n; i++)
            {
                if (a[i - 1 + (jcol - 1) * n] != 0.0)
                {
                    t = -a[i - 1 + (jcol - 1) * n];
                    a[i - 1 + (jcol - 1) * n] = 0.0;
                    for (j = jcol + 1; j <= n; j++)
                    {
                        a[i - 1 + (j - 1) * n] += t * a[jcol - 1 + (j - 1) * n];
                    }

                    for (j = 0; j < nb; j++)
                    {
                        x[i - 1 + j * n] += t * x[jcol - 1 + j * n];
                    }
                }
            }
        }

        //
        //  Back solve.
        //
        for (jcol = n; 2 <= jcol; jcol--)
        {
            for (i = 1; i < jcol; i++)
            {
                for (j = 0; j < nb; j++)
                {
                    x[i - 1 + j * n] -= a[i - 1 + (jcol - 1) * n] * x[jcol - 1 + j * n];
                }
            }
        }

    }

    public static Complex[] c8mat_fss_new(int n, ref Complex[] a, int nb,
            Complex[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_FSS_NEW factors and solves a system with multiple right hand sides.
        //
        //  Discussion:
        //
        //    This routine uses partial pivoting, but no pivot vector is required.
        //
        //    A C8MAT is a doubly dimensioned array of C8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 March 2013
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
        //    Input/output, Complex A[N*N].
        //    On input, A is the coefficient matrix of the linear system.
        //    On output, A is in unit upper triangular form, and
        //    represents the U factor of an LU factorization of the
        //    original coefficient matrix.
        //
        //    Input, int NB, the number of right hand sides.
        //
        //    Input, Complex B[N*NB], the right hand sides of the linear systems.
        //
        //    Output, Complex C8MAT_FSS_NEW[N*NB], the solutions of the 
        //    linear systems.
        //
    {
        int i;
        int j;
        int jcol;

        Complex[] x = new Complex[n * nb];

        for (j = 0; j < nb; j++)
        {
            for (i = 0; i < n; i++)
            {
                x[i + j * n] = b[i + j * n];
            }
        }

        for (jcol = 1; jcol <= n; jcol++)
        {
            //
            //  Find the maximum element in column I.
            //
            double piv = c8_abs(a[jcol - 1 + (jcol - 1) * n]);
            int ipiv = jcol;
            for (i = jcol + 1; i <= n; i++)
            {
                if (piv < c8_abs(a[i - 1 + (jcol - 1) * n]))
                {
                    piv = c8_abs(a[i - 1 + (jcol - 1) * n]);
                    ipiv = i;
                }
            }

            switch (piv)
            {
                case 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("C8MAT_FSS_NEW - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + jcol + "");
                    return null;
            }

            //
            //  Switch rows JCOL and IPIV, and X.
            //
            Complex t;
            if (jcol != ipiv)
            {
                for (j = 1; j <= n; j++)
                {
                    t = a[jcol - 1 + (j - 1) * n];
                    a[jcol - 1 + (j - 1) * n] = a[ipiv - 1 + (j - 1) * n];
                    a[ipiv - 1 + (j - 1) * n] = t;
                }

                for (j = 0; j < nb; j++)
                {
                    t = x[jcol - 1 + j * n];
                    x[jcol - 1 + j * n] = x[ipiv - 1 + j * n];
                    x[ipiv - 1 + j * n] = t;
                }
            }

            //
            //  Scale the pivot row.
            //
            t = a[jcol - 1 + (jcol - 1) * n];
            a[jcol - 1 + (jcol - 1) * n] = 1.0;
            for (j = jcol + 1; j <= n; j++)
            {
                a[jcol - 1 + (j - 1) * n] /= t;
            }

            for (j = 0; j < nb; j++)
            {
                x[jcol - 1 + j * n] /= t;
            }

            //
            //  Use the pivot row to eliminate lower entries in that column.
            //
            for (i = jcol + 1; i <= n; i++)
            {
                if (a[i - 1 + (jcol - 1) * n] != 0.0)
                {
                    t = -a[i - 1 + (jcol - 1) * n];
                    a[i - 1 + (jcol - 1) * n] = 0.0;
                    for (j = jcol + 1; j <= n; j++)
                    {
                        a[i - 1 + (j - 1) * n] += t * a[jcol - 1 + (j - 1) * n];
                    }

                    for (j = 0; j < nb; j++)
                    {
                        x[i - 1 + j * n] += t * x[jcol - 1 + j * n];
                    }
                }
            }
        }

        //
        //  Back solve.
        //
        for (jcol = n; 2 <= jcol; jcol--)
        {
            for (i = 1; i < jcol; i++)
            {
                for (j = 0; j < nb; j++)
                {
                    x[i - 1 + j * n] -= a[i - 1 + (jcol - 1) * n] * x[jcol - 1 + j * n];
                }
            }
        }

        return x;
    }

    public static Complex[] c8mat_identity_new(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_IDENTITY_NEW sets a C8MAT to the identity.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 November 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, Complex C8MAT_IDENTITY_NEW[N*N], the matrix.
        //
    {
        int j;

        Complex[] a = new Complex [n * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                if (i == j)
                {
                    a[i + j * n] = new Complex(1.0, 0.0);
                }
                else
                {
                    a[i + j * n] = new Complex(0.0, 0.0);
                }
            }
        }

        return a;
    }

    public static Complex[] c8mat_indicator_new(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_INDICATOR_NEW returns the C8MAT indicator matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 November 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Output, Complex C8MAT_INDICATOR_NEW[M*N], the matrix.
        //
    {
        int j;

        Complex[] a = new Complex [m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                a[i + j * m] = new Complex(i + 1, j + 1);
            }
        }

        return a;
    }

    public static void c8mat_minvm(int n1, int n2, Complex[] a,
            Complex[] b, ref Complex[] c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_MINVM returns inverse(A) * B for C8MAT's.
        //
        //  Discussion:
        //
        //    A C8MAT is an array of C8 values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 March 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the order of the matrices.
        //
        //    Input, Complex A[N1*N1], B[N1*N2], the matrices.
        //
        //    Output, Complex C[N1*N2], the result, 
        //    C = inverse(A) * B.
        //
    {
        Complex[] alu = c8mat_copy_new(n1, n1, a);

        c8mat_copy(n1, n2, b, ref c);

        c8mat_fss(n1, ref alu, n2, ref c);
    }

    public static Complex[] c8mat_minvm_new(int n1, int n2, Complex[] a,
            Complex[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_MINVM_NEW returns inverse(A) * B for C8MAT's.
        //
        //  Discussion:
        //
        //    A C8MAT is an array of C8 values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 March 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the order of the matrices.
        //
        //    Input, Complex A[N1*N1], B[N1*N2], the matrices.
        //
        //    Output, Complex C8MAT_MINVM_NEW[N1*N2], the result, 
        //    C = inverse(A) * B.
        //
    {
        Complex[] alu = c8mat_copy_new(n1, n1, a);
        Complex[] c = c8mat_fss_new(n1, ref alu, n2, b);

        return c;
    }

    public static void c8mat_mm(int n1, int n2, int n3, Complex[] a,
            Complex[] b, ref Complex[] c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_MM multiplies two matrices.
        //
        //  Discussion:
        //
        //    A C8MAT is a doubly dimensioned array of C8 values, stored as a vector
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
        //    03 March 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, N3, the order of the matrices.
        //
        //    Input, Complex A[N1*N2], Complex B[N2*N3], 
        //    the matrices to multiply.
        //
        //    Output, Complex C[N1*N3], the product matrix C = A * B.
        //
    {
        int i;

        Complex[] c1 = new Complex [n1 * n3];

        for (i = 0; i < n1; i++)
        {
            int j;
            for (j = 0; j < n3; j++)
            {
                c1[i + j * n1] = 0.0;
                int k;
                for (k = 0; k < n2; k++)
                {
                    c1[i + j * n1] += a[i + k * n1] * b[k + j * n2];
                }
            }
        }

        c8mat_copy(n1, n3, c1, ref c);
    }

    public static Complex[] c8mat_mm_new(int n1, int n2, int n3, Complex[] a,
            Complex[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_MM_NEW multiplies two matrices.
        //
        //  Discussion:
        //
        //    A C8MAT is a doubly dimensioned array of C8 values, stored as a vector
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
        //    25 June 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, N3, the order of the matrices.
        //
        //    Input, Complex A[N1*N2], Complex B[N2*N3], 
        //    the matrices to multiply.
        //
        //    Output, Complex C8MAT_MM_NEW[N1*N3], the product matrix C = A * B.
        //
    {
        int i;

        Complex[] c = new Complex [n1 * n3];

        for (i = 0; i < n1; i++)
        {
            int j;
            for (j = 0; j < n3; j++)
            {
                c[i + j * n1] = 0.0;
                int k;
                for (k = 0; k < n2; k++)
                {
                    c[i + j * n1] += a[i + k * n1] * b[k + j * n2];
                }
            }
        }

        return c;
    }

    public static void c8mat_nint(int m, int n, ref Complex[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_NINT rounds the entries of a C8MAT.
        //
        //  Discussion:
        //
        //    A C8MAT is an array of Complex values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 November 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of A.
        //
        //    Input/output, Complex A[M*N], the matrix to be NINT'ed.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                a[i + j * m] = c8_nint(a[i + j * m]);
            }
        }
    }

    public static double c8mat_norm_fro(int m, int n, Complex[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_NORM_FRO returns the Frobenius norm of a C8MAT.
        //
        //  Discussion:
        //
        //    A C8MAT is an array of C8 values.
        //
        //    The Frobenius norm is defined as
        //
        //      C8MAT_NORM_FRO = Math.Sqrt (
        //        sum ( 1 <= I <= M ) Sum ( 1 <= J <= N ) |A(I,J)| )
        //
        //    The matrix Frobenius-norm is not derived from a vector norm, but
        //    is compatible with the vector L2 norm, so that:
        //
        //      c8vec_norm_l2 ( A*x ) <= c8mat_norm_fro ( A ) * c8vec_norm_l2 ( x ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 March 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, Complex A[M*N], the matrix.
        //
        //    Output, double C8MAT_NORM_FRO, the Frobenius norm of A.
        //
    {
        int j;

        double value = 0.0;
        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                value = value + Math.Pow(a[i + j * m].Real, 2)
                              + Math.Pow(a[i + j * m].Imaginary, 2);
            }
        }

        value = Math.Sqrt(value);

        return value;
    }

    public static double c8mat_norm_l1(int m, int n, Complex[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_NORM_L1 returns the matrix L1 norm of a C8MAT.
        //
        //  Discussion:
        //
        //    A C8MAT is an MxN array of C8's, stored by (I,J) -> [I+J*M].
        //
        //    The matrix L1 norm is defined as:
        //
        //      C8MAT_NORM_L1 = max ( 1 <= J <= N )
        //        sum ( 1 <= I <= M ) abs ( A(I,J) ).
        //
        //    The matrix L1 norm is derived from the vector L1 norm, and
        //    satisifies:
        //
        //      c8vec_norm_l1 ( A * x ) <= c8mat_norm_l1 ( A ) * c8vec_norm_l1 ( x ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 March 2013
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
        //    Input, Complex A[M*N], the matrix whose L1 norm is desired.
        //
        //    Output, double C8MAT_NORM_L1, the L1 norm of A.
        //
    {
        int j;

        double value = 0.0;

        for (j = 0; j < n; j++)
        {
            double col_sum = 0.0;
            int i;
            for (i = 0; i < m; i++)
            {
                col_sum += c8_abs(a[i + j * m]);
            }

            value = Math.Max(value, col_sum);
        }

        return value;
    }

    public static double c8mat_norm_li(int m, int n, Complex[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_NORM_LI returns the matrix L-oo norm of a C8MAT.
        //
        //  Discussion:
        //
        //    A C8MAT is an array of C8 values.
        //
        //    The matrix L-oo norm is defined as:
        //
        //      C8MAT_NORM_LI =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).
        //
        //    The matrix L-oo norm is derived from the vector L-oo norm,
        //    and satisifies:
        //
        //      c8vec_norm_li ( A * x ) <= c8mat_norm_li ( A ) * c8vec_norm_li ( x ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 March 2013
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
        //    Input, Complex A[M*N], the matrix whose L-oo
        //    norm is desired.
        //
        //    Output, double C8MAT_NORM_LI, the L-oo norm of A.
        //
    {
        int i;

        double value = 0.0;

        for (i = 0; i < m; i++)
        {
            double row_sum = 0.0;
            int j;
            for (j = 0; j < n; j++)
            {
                row_sum += c8_abs(a[i + j * m]);
            }

            value = Math.Max(value, row_sum);
        }

        return value;
    }

    public static void c8mat_print(int m, int n, Complex[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_PRINT prints a C8MAT.
        //
        //  Discussion:
        //
        //    A C8MAT is an array of Complex values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the matrix.
        //
        //    Input, Complex A[M*N], the matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        c8mat_print_some(m, n, a, 1, 1, m, n, title);
    }

    public static void c8mat_print_some(int m, int n, Complex[] a, int ilo, int jlo,
            int ihi, int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_PRINT_SOME prints some of a C8MAT.
        //
        //  Discussion:
        //
        //    A C8MAT is an array of Complex values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the matrix.
        //
        //    Input, Complex A[M*N], the matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        int i2lo = 0;
        const int incx = 4;

        Console.WriteLine("");
        Console.WriteLine(title + "");

        if (m <= 0 || n <= 0)
        {
            Console.WriteLine("");
            Console.WriteLine("  (None)");
            return;
        }

        //
        //  Print the columns of the matrix, in strips of INCX.
        //
        // for (j2lo = jlo; j2lo <= jhi; j2lo = j2lo + incx)
        {
            int j2hi = jlo + incx - 1;
            if (n < j2hi)
            {
                j2hi = n;
            }

            if (jhi < j2hi)
            {
                j2hi = jhi;
            }

            int inc = j2hi + 1 - jlo;

            Console.WriteLine("");
            string cout = "  Col: ";
            int j;
            int j2;
            for (j = jlo; j <= j2hi; j++)
            {
                j2 = j + 1 - jlo;
                cout += "     " + j.ToString().PadLeft(10) + "     ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Row");
            Console.WriteLine("  ---");
            i2lo = i2lo switch
            {
                < 1 => 1,
                //
                //  Determine the range of the rows in this strip.
                //
                _ => ilo
            };

            int i2hi = ihi;
            if (m < i2hi)
            {
                i2hi = m;
            }

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                cout = i.ToString().PadLeft(5) + ":";
                //
                //  Print out (up to) INCX entries in row I, that lie in the current strip.
                //
                for (j2 = 1; j2 <= inc; j2++)
                {
                    j = jlo - 1 + j2;
                    Complex c = a[i - 1 + (j - 1) * m];
                    cout += "  " + c.Real.ToString().PadLeft(8)
                                 + "  " + c.Imaginary.ToString().PadLeft(8);
                }

                Console.WriteLine(cout);
            }
        }

    }

    public static void c8mat_scale(int m, int n, Complex alpha, ref Complex[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_SCALE scales a C8MAT by a complex scalar factor.
        //
        //  Discussion:
        //
        //    An C8MAT is a doubly dimensioned array of complex double precision values, 
        //    which may be stored as a vector in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 March 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, Complex ALPHA, the scale factor.
        //
        //    Input/output, Complex A[M*N], the matrix to be scaled.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                a[i + j * m] *= alpha;
            }
        }
    }

    public static void c8mat_scale_r8(int m, int n, double alpha, ref Complex[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_SCALE_R8 scales a C8MAT by a real scalar factor.
        //
        //  Discussion:
        //
        //    An C8MAT is a doubly dimensioned array of complex double precision values, 
        //    which may be stored as a vector in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 March 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double ALPHA, the scale factor.
        //
        //    Input/output, Complex A[M*N], the matrix to be scaled.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                a[i + j * m] *= alpha;
            }
        }
    }

    public static void c8mat_uniform_01(int m, int n, ref int seed, ref Complex[] c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_UNIFORM_01 returns a unit pseudorandom C8MAT.
        //
        //  Discussion:
        //
        //    A C8MAT is an array of Complex values.
        //
        //    The angles should be uniformly distributed between 0 and 2 * PI,
        //    the square roots of the radius uniformly distributed between 0 and 1.
        //
        //    This results in a uniform distribution of values in the unit circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the matrix.
        //
        //    Input/output, int &SEED, the "seed" value, which should NOT be 0.
        //    On output, SEED has been updated.
        //
        //    Output, Complex C[M*N], the pseudorandom complex matrix.
        //
    {
        int j;

        switch (seed)
        {
            case 0:
                Console.WriteLine("");
                Console.WriteLine("C8MAT_UNIFORM_01 - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return;
        }

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                int k = seed / 127773;

                seed = 16807 * (seed - k * 127773) - k * 2836;

                switch (seed)
                {
                    case < 0:
                        seed += i4_huge();
                        break;
                }

                double r = Math.Sqrt(seed * 4.656612875E-10);

                k = seed / 127773;

                seed = 16807 * (seed - k * 127773) - k * 2836;

                switch (seed)
                {
                    case < 0:
                        seed += i4_huge();
                        break;
                }

                double theta = 2.0 * Math.PI * (seed * 4.656612875E-10);

                c[i + j * m] = r * new Complex(Math.Cos(theta), Math.Sin(theta));
            }
        }
    }

    public static Complex[] c8mat_uniform_01_new(int m, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_UNIFORM_01_NEW returns a unit pseudorandom C8MAT.
        //
        //  Discussion:
        //
        //    A C8MAT is an array of Complex values.
        //
        //    The angles should be uniformly distributed between 0 and 2 * PI,
        //    the square roots of the radius uniformly distributed between 0 and 1.
        //
        //    This results in a uniform distribution of values in the unit circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the matrix.
        //
        //    Input/output, int &SEED, the "seed" value, which should NOT be 0.
        //    On output, SEED has been updated.
        //
        //    Output, Complex C8MAT_UNIFORM_01_NEW[M*N], the pseudorandom 
        //    complex matrix.
        //
    {
        int j;

        Complex[] c = new Complex [m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                int k = seed / 127773;

                seed = 16807 * (seed - k * 127773) - k * 2836;

                switch (seed)
                {
                    case < 0:
                        seed += 2147483647;
                        break;
                }

                double r = Math.Sqrt(seed * 4.656612875E-10);

                k = seed / 127773;

                seed = 16807 * (seed - k * 127773) - k * 2836;

                switch (seed)
                {
                    case < 0:
                        seed += 2147483647;
                        break;
                }

                double theta = 2.0 * Math.PI * (seed * 4.656612875E-10);

                c[i + j * m] = r * new Complex(Math.Cos(theta), Math.Sin(theta));
            }
        }

        return c;
    }

    public static Complex[] c8mat_zero_new(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_ZERO_NEW returns a new zeroed C8MAT.
        //
        //  Discussion:
        //
        //    An C8MAT is a doubly dimensioned array of complex double precision values, 
        //    which may be stored as a vector in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 March 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Output, Complex C8MAT_ZERO_NEW[M*N], the zeroed matrix.
        //
    {
        int j;

        Complex[] a = new Complex[m * n];

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