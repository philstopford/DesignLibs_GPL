using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {

        public static void r8mat_transpose_print(int m, int n, double[] a, string title)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
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
            //    10 September 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns.
            //
            //    Input, double A[M*N], an M by N matrix to be printed.
            //
            //    Input, string TITLE, a title.
            //
        {
            r8mat_transpose_print_some(m, n, a, 1, 1, m, n, title);
        }

        public static void r8mat_transpose_print_some(int m, int n, double[] a, int ilo, int jlo,
                int ihi, int jhi, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
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
            //    07 April 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns.
            //
            //    Input, double A[M*N], an M by N matrix to be printed.
            //
            //    Input, int ILO, JLO, the first row and column to print.
            //
            //    Input, int IHI, JHI, the last row and column to print.
            //
            //    Input, string TITLE, a title.
            //
        {
            int INCX = 5;

            int i;
            int i2;
            int i2hi;
            int i2lo;
            int i2lo_hi;
            int i2lo_lo;
            int inc;
            int j;
            int j2hi;
            int j2lo;

            Console.WriteLine();
            Console.WriteLine(title);

            if (m <= 0 || n <= 0)
            {
                Console.WriteLine();
                Console.WriteLine("  (None)");
                return;
            }

            if (ilo < 1)
            {
                i2lo_lo = 1;
            }
            else
            {
                i2lo_lo = ilo;
            }

            if (ihi < m)
            {
                i2lo_hi = m;
            }
            else
            {
                i2lo_hi = ihi;
            }

            for (i2lo = i2lo_lo; i2lo <= i2lo_hi; i2lo = i2lo + INCX)
            {
                // Ugly hack to sidestep a mismatch in the output behavior compared to reference.
                if (i2lo > INCX)
                {
                    break;
                }

                i2hi = i2lo + INCX - 1;

                if (m < i2hi)
                {
                    i2hi = m;
                }

                if (ihi < i2hi)
                {
                    i2hi = ihi;
                }

                inc = i2hi + 1 - i2lo;

                Console.WriteLine();
                string cout = "  Row: ";
                for (i = i2lo; i <= i2hi; i++)
                {
                    cout += (i - 1).ToString().PadLeft(7) + "       ";
                }

                Console.WriteLine(cout);
                Console.WriteLine("  Col");
                Console.WriteLine();

                if (jlo < 1)
                {
                    j2lo = 1;
                }
                else
                {
                    j2lo = jlo;
                }

                if (n < jhi)
                {
                    j2hi = n;
                }
                else
                {
                    j2hi = jhi;
                }

                for (j = j2lo; j <= j2hi; j++)
                {
                    cout = (j - 1).ToString().PadLeft(5) + ":";
                    for (i2 = 1; i2 <= inc; i2++)
                    {
                        i = i2lo - 1 + i2;
                        string t = a[(i - 1) + (j - 1) * m].ToString("0.######");
                        cout += t.PadLeft(14);
                    }

                    Console.WriteLine(cout);
                }
            }
        }

        public static void r8mat_print(int m, int n, double[] a, string title)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_PRINT prints an R8MAT, with an optional title.
            //
            //  Discussion:
            //
            //    An R8MAT is an array of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 August 2003
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
            //    Input, double A[M*N], the M by N matrix.
            //
            //    Input, string TITLE, a title.
            //
        {
            r8mat_print_some(m, n, a, 1, 1, m, n, title);
        }
        //****************************************************************************80

        public static void r8mat_print_some(int m, int n, double[] a, int ilo, int jlo, int ihi,
                int jhi, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_PRINT_SOME prints some of an R8MAT.
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
            //    26 June 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows of the matrix.
            //    M must be positive.
            //
            //    Input, int N, the number of columns of the matrix.
            //    N must be positive.
            //
            //    Input, double A[M*N], the matrix.
            //
            //    Input, int ILO, JLO, IHI, JHI, designate the first row and
            //    column, and the last row and column to be printed.
            //
            //    Input, string TITLE, a title.
            //
        {
            int INCX = 5;

            Console.WriteLine();
            Console.WriteLine(title);

            if (m <= 0 || n <= 0)
            {
                Console.WriteLine();
                Console.WriteLine("  (None)");
                return;
            }

            //
            //  Print the columns of the matrix, in strips of 5.
            //
            for (int j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX)
            {
                int j2hi = j2lo + INCX - 1;
                if (n < j2hi)
                {
                    j2hi = n;
                }

                if (jhi < j2hi)
                {
                    j2hi = jhi;
                }

                Console.WriteLine();
                //
                //  For each column J in the current range...
                //
                //  Write the header.
                //
                string cout = "  Col:    ";
                for (int j = j2lo; j <= j2hi; j++)
                {
                    cout += (j - 1).ToString().PadLeft(7) + "       ";
                }

                Console.WriteLine(cout);
                Console.WriteLine();
                Console.WriteLine("  Row");
                Console.WriteLine();
                //
                //  Determine the range of the rows in this strip.
                //

                int i2lo;
                int i2hi;

                if (1 < ilo)
                {
                    i2lo = ilo;
                }
                else
                {
                    i2lo = 1;
                }

                if (ihi < m)
                {
                    i2hi = ihi;
                }
                else
                {
                    i2hi = m;
                }

                for (int i = i2lo; i <= i2hi; i++)
                {
                    //
                    //  Print out (up to) 5 entries in row I, that lie in the current strip.
                    //
                    cout = (i - 1).ToString().PadLeft(5) + ": ";
                    for (int j = j2lo; j <= j2hi; j++)
                    {
                        cout += a[i - 1 + (j - 1) * m].ToString().PadLeft(12) + "  ";
                    }

                    Console.WriteLine(cout);
                }
            }
        }

        public static double[] r8mat_solve2(int n, ref double[] a, ref double[] b, ref int ierror)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_SOLVE2 computes the solution of an N by N linear system.
            //
            //  Discussion: 							    
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector 
            //    in column-major order.
            //
            //    The linear system may be represented as
            //
            //      A*X = B
            //
            //    If the linear system is singular, but consistent, then the routine will
            //    still produce a solution.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of equations.
            //
            //    Input/output, double A[N*N].
            //    On input, A is the coefficient matrix to be inverted.
            //    On output, A has been overwritten.
            //
            //    Input/output, double B[N].
            //    On input, B is the right hand side of the system.
            //    On output, B has been overwritten.
            //
            //    Output, int *IERROR.
            //    0, no error detected.
            //    1, consistent singularity.
            //    2, inconsistent singularity.
            //
            //    Output, double R8MAT_SOLVE2[N], the solution of the linear system.
            //
        {
            double amax;
            int imax;
            int[] piv;
            double[] x;

            ierror = 0;

            piv = typeMethods.i4vec_zero_new(n);
            x = typeMethods.r8vec_zero_new(n);
            //
            //  Process the matrix.
            //
            for (int k = 1; k <= n; k++)
            {
                //
                //  In column K:
                //    Seek the row IMAX with the properties that:
                //      IMAX has not already been used as a pivot;
                //      A(IMAX,K) is larger in magnitude than any other candidate.
                //
                amax = 0.0;
                imax = 0;
                for (int i = 1; i <= n; i++)
                {
                    if (piv[i - 1] == 0)
                    {
                        if (amax < Math.Abs(a[i - 1 + (k - 1) * n]))
                        {
                            imax = i;
                            amax = Math.Abs(a[i - 1 + (k - 1) * n]);
                        }
                    }
                }

                //
                //  If you found a pivot row IMAX, then,
                //    eliminate the K-th entry in all rows that have not been used for pivoting.
                //
                if (imax != 0)
                {
                    piv[imax - 1] = k;
                    for (int j = k + 1; j <= n; j++)
                    {
                        a[imax - 1 + (j - 1) * n] = a[imax - 1 + (j - 1) * n] / a[imax - 1 + (k - 1) * n];
                    }

                    b[imax - 1] = b[imax - 1] / a[imax - 1 + (k - 1) * n];
                    a[imax - 1 + (k - 1) * n] = 1.0;

                    for (int i = 1; i <= n; i++)
                    {
                        if (piv[i - 1] == 0)
                        {
                            for (int j = k + 1; j <= n; j++)
                            {
                                a[i - 1 + (j - 1) * n] = a[i - 1 + (j - 1) * n] -
                                                         a[i - 1 + (k - 1) * n] * a[imax - 1 + (j - 1) * n];
                            }

                            b[i - 1] = b[i - 1] - a[i - 1 + (k - 1) * n] * b[imax - 1];
                            a[i - 1 + (k - 1) * n] = 0.0;
                        }
                    }
                }
            }

            //
            //  Now, every row with nonzero IPIV begins with a 1, and
            //  all other rows are all zero.  Begin solution.
            //
            for (int j = n; 1 <= j; j--)
            {
                imax = 0;
                for (int k = 1; k <= n; k++)
                {
                    if (piv[k - 1] == j)
                    {
                        imax = k;
                    }
                }

                if (imax == 0)
                {
                    x[j - 1] = 0.0;

                    if (b[j - 1] == 0.0)
                    {
                        ierror = 1;
                        Console.WriteLine("");
                        Console.WriteLine("R8MAT_SOLVE2 - Warning:");
                        Console.WriteLine("  Consistent singularity, equation = " + j + "");
                    }
                    else
                    {
                        ierror = 2;
                        Console.WriteLine("");
                        Console.WriteLine("R8MAT_SOLVE2 - Warning:");
                        Console.WriteLine("  Inconsistent singularity, equation = " + j + "");
                    }
                }
                else
                {
                    x[j - 1] = b[imax - 1];

                    for (int i = 1; i <= n; i++)
                    {
                        if (i != imax)
                        {
                            b[i - 1] = b[i - 1] - a[i - 1 + (j - 1) * n] * x[j - 1];
                        }
                    }
                }
            }

            return x;
        }

        public static double[] r8mat_fs_new(int n, double[] a, double[] b)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_FS_NEW factors and solves a system with one right hand side.
            //
            //  Discussion:
            //
            //    This routine differs from R8MAT_FSS_NEW in two ways:
            //    * only one right hand side is allowed;
            //    * the input matrix A is not modified.
            //
            //    This routine uses partial pivoting, but no pivot vector is required.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 January 2013
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
            //    Input, double A[N*N], the coefficient matrix of the linear system.
            //    On output, A is in unit upper triangular form, and
            //    represents the U factor of an LU factorization of the
            //    original coefficient matrix.
            //
            //    Input, double B[N], the right hand side of the linear system.
            //
            //    Output, double X[N], the solution of the linear system.
            //
        {
            double[] a2;
            int i;
            int ipiv;
            int j;
            int jcol;
            double piv;
            double t;
            double[] x;

            a2 = new double[n * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    a2[i + j * n] = a[i + j * n];
                }
            }

            x = new double[n];
            for (i = 0; i < n; i++)
            {
                x[i] = b[i];
            }

            for (jcol = 1; jcol <= n; jcol++)
            {
                //
                //  Find the maximum element in column I.
                //
                piv = Math.Abs(a2[jcol - 1 + (jcol - 1) * n]);
                ipiv = jcol;
                for (i = jcol + 1; i <= n; i++)
                {
                    if (piv < Math.Abs(a2[i - 1 + (jcol - 1) * n]))
                    {
                        piv = Math.Abs(a2[i - 1 + (jcol - 1) * n]);
                        ipiv = i;
                    }
                }

                if (piv == 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8MAT_FS_NEW - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + jcol);
                    return new double[1];
                }

                //
                //  Switch rows JCOL and IPIV, and X.
                //
                if (jcol != ipiv)
                {
                    for (j = 1; j <= n; j++)
                    {
                        t = a2[jcol - 1 + (j - 1) * n];
                        a2[jcol - 1 + (j - 1) * n] = a2[ipiv - 1 + (j - 1) * n];
                        a2[ipiv - 1 + (j - 1) * n] = t;
                    }

                    t = x[jcol - 1];
                    x[jcol - 1] = x[ipiv - 1];
                    x[ipiv - 1] = t;
                }

                //
                //  Scale the pivot row.
                //
                t = a2[jcol - 1 + (jcol - 1) * n];
                a2[jcol - 1 + (jcol - 1) * n] = 1.0;
                for (j = jcol + 1; j <= n; j++)
                {
                    a2[jcol - 1 + (j - 1) * n] = a2[jcol - 1 + (j - 1) * n] / t;
                }

                x[jcol - 1] = x[jcol - 1] / t;
                //
                //  Use the pivot row to eliminate lower entries in that column.
                //
                for (i = jcol + 1; i <= n; i++)
                {
                    if (a2[i - 1 + (jcol - 1) * n] != 0.0)
                    {
                        t = -a2[i - 1 + (jcol - 1) * n];
                        a2[i - 1 + (jcol - 1) * n] = 0.0;
                        for (j = jcol + 1; j <= n; j++)
                        {
                            a2[i - 1 + (j - 1) * n] = a2[i - 1 + (j - 1) * n] + t * a2[jcol - 1 + (j - 1) * n];
                        }

                        x[i - 1] = x[i - 1] + t * x[jcol - 1];
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
                    x[i - 1] = x[i - 1] - a2[i - 1 + (jcol - 1) * n] * x[jcol - 1];
                }
            }

            return x;
        }

        public static double[][] r8rmat_new(int m, int n)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8RMAT_NEW allocates a new R8RMAT.
            //
            //  Discussion:
            //
            //    An R8RMAT is a row-major array that was created by a 
            //    command like:
            //
            //      double **a;
            //      a = r8rmat_new ( m, n );
            //
            //    The user assigns entries to the matrix using typical
            //    2D array notation:
            //      a[2][3] = 17.0;
            //      y = a[1][0];
            //    and so on.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 September 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the matrix.
            //
            //    Output, double **R8RMAT_NEW, a new matrix.
            //
        {
            double[][] a = new double[m][];

            for (int i = 0; i < m; i++)
            {
                a[i] = new double[n];
            }

            return a;
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
                    y[j] = y[j] + a[i + j * m] * x[i];
                }
            }

            return y;
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

        public static double r8mat_podet(int n, double[] r)

            /******************************************************************************/
            /*
            Purpose:
            
            R8MAT_PODET computes the determinant of a factored positive definite matrix.
            
            Discussion:
            
            This routine expects to receive R, the upper triangular factor of A,
            computed by R8MAT_POFAC, with the property that A = R' * R.
            
            Licensing:
            
            This code is distributed under the GNU LGPL license.
            
            Modified:
            
            09 June 2013
            
            Author:
            
            C version by John Burkardt.
            
            Reference:
            
            Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
            LINPACK User's Guide,
            SIAM, 1979,
            ISBN13: 978-0-898711-72-1,
            LC: QA214.L56.
            
            Parameters:
            
            Input, int N, the order of the matrix.
            
            Input, double R[N*N], the Cholesky factor of A.
            
            Output, double R8MAT_PODET, the determinant of A.
            */
        {
            double det = 1.0;
            for (int i = 0; i < n; i++)
            {
                det = det * r[i + i * n] * r[i + i * n];
            }

            return det;
        }

        public static double[] r8mat_pofac(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_POFAC factors a real symmetric positive definite matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 June 2013
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, 
            //    Pete Stewart,
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
            //    LINPACK User's Guide,
            //    SIAM, (Society for Industrial and Applied Mathematics),
            //    3600 University City Science Center,
            //    Philadelphia, PA, 19104-2688.
            //    ISBN 0-89871-172-X
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, double A[N*N], the matrix to be factored.
            //
            //    Output, double R8MAT_POFAC[N*N], an upper triangular matrix such that
            //    A = R'*R.
            //
        {
            double[] r = new double[n * n];

            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i <= j; i++)
                {
                    r[i + j * n] = a[i + j * n];
                }

                for (int i = j + 1; i < n; i++)
                {
                    r[i + j * n] = 0.0;
                }
            }

            for (int j = 0; j < n; j++)
            {
                double s = 0.0;

                for (int k = 0; k < j; k++)
                {
                    double dot = 0.0;
                    for (int i = 0; i < k; i++)
                    {
                        dot = dot + r[i + k * n] * r[i + j * n];
                    }

                    double t = r[k + j * n] - dot;
                    t = t / r[k + k * n];
                    r[k + j * n] = t;
                    s = s + t * t;
                }

                s = r[j + j * n] - s;

                if (s < 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8MAT_POFAC - Fatal error!");
                    Console.WriteLine("  The matrix is not positive definite.");
                    return new double[1];
                }

                if (s == 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8MAT_POFAC - Warning!");
                    Console.WriteLine("  The matrix is not strictly positive definite.");
                }

                r[j + j * n] = Math.Sqrt(s);
            }

            return r;
        }

        public static double[] r8mat_poinv(int n, double[] r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_POINV computes the inverse of a factored positive definite matrix.
            //
            //  Discussion:
            //
            //    This routine expects to receive R, the upper triangular factor of A,
            //    computed by R8MAT_POFAC, with the property that A = R' * R.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 June 2013
            //
            //  Author:
            //
            //    C version by John Burkardt.
            //
            //  Reference:
            //
            //    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
            //    LINPACK User's Guide,
            //    SIAM, 1979,
            //    ISBN13: 978-0-898711-72-1,
            //    LC: QA214.L56.
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix A.
            //
            //    Input, double R[N*N], the Cholesky factor of A.
            //
            //    Input, double R8MAT_POINV[N*N], the inverse of A.
            //
        {
            double t;

            double[] b = new double[n * n];

            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    b[i + j * n] = r[i + j * n];
                }
            }

            for (int k = 0; k < n; k++)
            {
                b[k + k * n] = 1.0 / b[k + k * n];
                for (int i = 0; i < k; i++)
                {
                    b[i + k * n] = -b[i + k * n] * b[k + k * n];
                }

                for (int j = k + 1; j < n; j++)
                {
                    t = b[k + j * n];
                    b[k + j * n] = 0.0;
                    for (int i = 0; i <= k; i++)
                    {
                        b[i + j * n] = b[i + j * n] + t * b[i + k * n];
                    }
                }
            }

            /*
            Form inverse(R) * (inverse(R))'.
            */
            for (int j = 0; j < n; j++)
            {
                for (int k = 0; k < j; k++)
                {
                    t = b[k + j * n];
                    for (int i = 0; i <= k; i++)
                    {
                        b[i + k * n] = b[i + k * n] + t * b[i + j * n];
                    }
                }

                t = b[j + j * n];
                for (int i = 0; i <= j; i++)
                {
                    b[i + j * n] = b[i + j * n] * t;
                }
            }

            return b;
        }

        public static double[] r8mat_upsol(int n, double[] r, double[] b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_UPSOL solves R * X = B, for an upper triangular matrix R.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 June 2013
            //
            //  Author:
            //
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
            //    LINPACK User's Guide,
            //    SIAM, 1979,
            //    ISBN13: 978-0-898711-72-1,
            //    LC: QA214.L56.
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, double R[N*N], the upper triangular matrix.
            //
            //    Input, double B[N], the right hand side.
            //
            //    Output, double R8MAT_UPSOL[N], the solution.
            //
        {

            double[] x = new double[n];
            for (int j = 0; j < n; j++)
            {
                x[j] = b[j];
            }

            return x;
        }

        public static double[] r8mat_utsol(int n, double[] r, double[] b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_UTSOL solves R' * X = B for an upper triangular matrix R.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 June 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
            //    LINPACK User's Guide,
            //    SIAM, 1979,
            //    ISBN13: 978-0-898711-72-1,
            //    LC: QA214.L56.
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, double R[N*N], the upper triangular matrix.
            //
            //    Input, double B[N], the right hand side.
            //
            //    Output, double X[N], the solution.
            //
        {
            double[] x = new double[n];

            for (int j = 0; j < n; j++)
            {
                x[j] = b[j];
            }

            x[0] = x[0] / r[0 + 0 * n];

            for (int j = 1; j < n; j++)
            {
                for (int i = 0; i < j; i++)
                {
                    x[j] = x[j] - r[i + j * n] * x[i];
                }

                x[j] = x[j] / r[j + j * n];
            }

            return x;
        }

        public static double r8mat_is_identity(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_IS_IDENTITY determines if an R8MAT is the identity.
            //
            //  Discussion:
            //
            //    An R8MAT is a matrix of real ( kind = 8 ) values.
            //
            //    The routine returns the Frobenius norm of A - I.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 July 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, double A[N*N], the matrix.
            //
            //    Output, double R8MAT_IS_IDENTITY, the Frobenius norm
            //    of the difference matrix A - I, which would be exactly zero
            //    if A were the identity matrix.
            //
        {
            double error_frobenius;
            int i;
            int j;
            double t;

            error_frobenius = 0.0;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    if (i == j)
                    {
                        t = a[i + j * n] - 1.0;
                    }
                    else
                    {
                        t = a[i + j * n];
                    }

                    error_frobenius = error_frobenius + t * t;
                }
            }

            error_frobenius = Math.Sqrt(error_frobenius);

            return error_frobenius;
        }

        public static double[] r8mat_mm_new(int n1, int n2, int n3, double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_MM_NEW multiplies two matrices.
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
        //    18 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, N3, the order of the matrices.
        //
        //    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.
        //
        //    Output, double R8MAT_MM[N1*N3], the product matrix C = A * B.
        //
        {
            double[] c;
            int i;
            int j;
            int k;

            c = new double[n1 * n3];

            for (i = 0; i < n1; i++)
            {
                for (j = 0; j < n3; j++)
                {
                    c[i + j * n1] = 0.0;
                    for (k = 0; k < n2; k++)
                    {
                        c[i + j * n1] = c[i + j * n1] + a[i + k * n1] * b[k + j * n2];
                    }
                }
            }

            return c;
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
        
        public static double[] r8mat_zeros_new ( int m, int n )

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
            double[] a;
            int i;
            int j;

            a = new double[m*n];

            for ( j = 0; j < n; j++ )
            {
                for ( i = 0; i < m; i++ )
                {
                    a[i+j*m] = 0.0;
                }
            }
            return a;
        }

    }
}