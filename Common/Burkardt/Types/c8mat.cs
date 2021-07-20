using System;
using System.Numerics;

namespace Burkardt.Types
{

    public static partial class typeMethods
    {
        public static void c8mat_add_r8(int m, int n, double alpha, Complex[] a,
                double beta, Complex[] b, Complex[] c)

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
            int i;
            int j;

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    c[i + j * m] = alpha * a[i + j * m] + beta * b[i + j * m];
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
            Complex[] a2;
            int i;
            int j;

            a2 = new Complex[m * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    a2[i + j * m] = a1[i + j * m];
                }
            }

            return a2;
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
            int j;
            double row_sum;
            double value;

            value = 0.0;

            for (i = 0; i < m; i++)
            {
                row_sum = 0.0;
                for (j = 0; j < n; j++)
                {
                    row_sum = row_sum + Complex.Abs(a[i + j * m]);
                }

                value = Math.Max(value, row_sum);
            }

            return value;
        }

        public static void c8mat_mm(int n1, int n2, int n3, Complex[] a,
                Complex[] b, Complex[] c)

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
            Complex[] c1;
            int i;
            int j;
            int k;

            c1 = new Complex [n1 * n3];

            for (i = 0; i < n1; i++)
            {
                for (j = 0; j < n3; j++)
                {
                    c1[i + j * n1] = 0.0;
                    for (k = 0; k < n2; k++)
                    {
                        c1[i + j * n1] = c1[i + j * n1] + a[i + k * n1] * b[k + j * n2];
                    }
                }
            }

            c8mat_copy(n1, n3, c1, ref c);
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
            Complex[] alu;

            alu = c8mat_copy_new(n1, n1, a);

            c8mat_copy(n1, n2, b, ref c);

            c8mat_fss(n1, ref alu, n2, ref c);
        }
        
        public static void c8mat_scale_r8 ( int m, int n, double alpha, ref Complex[] a )

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
            //    Input/output, complex <double> A[M*N], the matrix to be scaled.
            //
        {
            int i;
            int j;

            for ( j = 0; j < n; j++ )
            {
                for ( i = 0; i < m; i++ )
                {
                    a[i+j*m] = a[i+j*m] * alpha;
                }
            }
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
            int ipiv;
            int j;
            int jcol;
            double piv;
            Complex t;

            for (jcol = 1; jcol <= n; jcol++)
            {
                //
                //  Find the maximum element in column I.
                //
                piv = Complex.Abs(a[jcol - 1 + (jcol - 1) * n]);
                ipiv = jcol;
                for (i = jcol + 1; i <= n; i++)
                {
                    if (piv < Complex.Abs(a[i - 1 + (jcol - 1) * n]))
                    {
                        piv = Complex.Abs(a[i - 1 + (jcol - 1) * n]);
                        ipiv = i;
                    }
                }

                if (piv == 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("C8MAT_FSS - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + jcol + "");
                    return;
                }

                //
                //  Switch rows JCOL and IPIV, and X.
                //
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
                    a[jcol - 1 + (j - 1) * n] = a[jcol - 1 + (j - 1) * n] / t;
                }

                for (j = 0; j < nb; j++)
                {
                    x[jcol - 1 + j * n] = x[jcol - 1 + j * n] / t;
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
                            a[i - 1 + (j - 1) * n] = a[i - 1 + (j - 1) * n] + t * a[jcol - 1 + (j - 1) * n];
                        }

                        for (j = 0; j < nb; j++)
                        {
                            x[i - 1 + j * n] = x[i - 1 + j * n] + t * x[jcol - 1 + j * n];
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
                        x[i - 1 + j * n] = x[i - 1 + j * n] - a[i - 1 + (jcol - 1) * n] * x[jcol - 1 + j * n];
                    }
                }
            }

            return;
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
            int i;
            int j;

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    a2[i + j * m] = a1[i + j * m];
                }
            }
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
            Complex[] a;
            int i;
            int j;

            a = new Complex [n * n];

            for (j = 0; j < n; j++)
            {
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
            Complex c;
            int i;
            int i2hi;
            int i2lo;
            int inc;
            int incx = 4;
            int j;
            int j2;
            int j2hi;
            int j2lo;

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
            for (j2lo = jlo; j2lo <= Math.Min(jhi, n); j2lo = j2lo + incx)
            {
                j2hi = j2lo + incx - 1;
                j2hi = Math.Min(j2hi, n);
                j2hi = Math.Min(j2hi, jhi);

                inc = j2hi + 1 - j2lo;

                Console.WriteLine("");
                string cout = "  Col: ";
                for (j = j2lo; j <= j2hi; j++)
                {
                    j2 = j + 1 - j2lo;
                    cout += "     " + j.ToString().PadLeft(10) + "     ";
                }

                Console.WriteLine(cout);
                Console.WriteLine("  Row");
                Console.WriteLine("  ---");
                //
                //  Determine the range of the rows in this strip.
                //
                i2lo = Math.Max(ilo, 1);
                i2hi = Math.Min(ihi, m);

                for (i = i2lo; i <= i2hi; i++)
                {
                    cout = i.ToString().PadLeft(5) + ":";
                    //
                    //  Print out (up to) INCX entries in row I, that lie in the current strip.
                    //
                    for (j2 = 1; j2 <= inc; j2++)
                    {
                        j = j2lo - 1 + j2;
                        c = a[i - 1 + (j - 1) * m];
                        cout += "  " + c.Real.ToString().PadLeft(8)
                                     + "  " + c.Imaginary.ToString().PadLeft(8);
                    }

                    Console.WriteLine(cout);
                }
            }
        }

    }
}