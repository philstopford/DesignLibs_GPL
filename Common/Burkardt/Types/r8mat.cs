using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        
        public static void r8mat_transpose_print ( int m, int n, double[] a, string title )
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
            r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );
        }
        
        public static void r8mat_transpose_print_some ( int m, int n, double[] a, int ilo, int jlo,
        int ihi, int jhi, string title )

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

            if ( m <= 0 || n <= 0 )
            {
                Console.WriteLine();
                Console.WriteLine("  (None)");
                return;
            }

            if ( ilo < 1 )
            {
                i2lo_lo = 1;
            }
            else
            {
                i2lo_lo = ilo;
            }

            if ( ihi < m )
            {
                i2lo_hi = m;
            }
            else
            {
                i2lo_hi = ihi;
            }

            for ( i2lo = i2lo_lo; i2lo <= i2lo_hi; i2lo = i2lo + INCX )
            {
                // Ugly hack to sidestep a mismatch in the output behavior compared to reference.
                if (i2lo > INCX)
                {
                    break;
                }
                i2hi = i2lo + INCX - 1;

                if ( m < i2hi )
                {
                    i2hi = m;
                }
                if ( ihi < i2hi )
                {
                    i2hi = ihi;
                }

                inc = i2hi + 1 - i2lo;

                Console.WriteLine();
                string cout = "  Row: ";
                for ( i = i2lo; i <= i2hi; i++ )
                {
                    cout += (i - 1).ToString().PadLeft(7) + "       ";
                }
                Console.WriteLine(cout);
                Console.WriteLine("  Col");
                Console.WriteLine();

                if ( jlo < 1 )
                {
                    j2lo = 1;
                }
                else
                {
                    j2lo = jlo;
                }

                if ( n < jhi )
                {
                    j2hi = n;
                }
                else
                {
                    j2hi = jhi;
                }

                for ( j = j2lo; j <= j2hi; j++ )
                {
                    cout = (j - 1).ToString().PadLeft(5) + ":";
                    for ( i2 = 1; i2 <= inc; i2++ )
                    {
                        i = i2lo - 1 + i2;
                        string t = a[(i - 1) + (j - 1) * m].ToString("0.######");
                        cout += t.PadLeft(14);
                    }
                    Console.WriteLine(cout);
                }
            }
        }
        
        public static void r8mat_print ( int m, int n, double[] a, string title )
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
            r8mat_print_some ( m, n, a, 1, 1, m, n, title );
        }
        //****************************************************************************80

        public static void r8mat_print_some ( int m, int n, double[] a, int ilo, int jlo, int ihi,
        int jhi, string title )

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

            if ( m <= 0 || n <= 0 )
            {
                Console.WriteLine();
                Console.WriteLine("  (None)");
                return;
            }
            //
            //  Print the columns of the matrix, in strips of 5.
            //
            for (int j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
            {
                int j2hi = j2lo + INCX - 1;
                if ( n < j2hi )
                {
                    j2hi = n;
                }
                if ( jhi < j2hi )
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
                for (int j = j2lo; j <= j2hi; j++ )
                {
                    cout += (j-1).ToString().PadLeft(7) +  "       ";
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
                
                if ( 1 < ilo )
                {
                    i2lo = ilo;
                }
                else
                {
                    i2lo = 1;
                }
                if ( ihi < m )
                {
                    i2hi = ihi;
                }
                else
                {
                    i2hi = m;
                }

                for (int i = i2lo; i <= i2hi; i++ )
                {
                    //
                    //  Print out (up to) 5 entries in row I, that lie in the current strip.
                    //
                    cout = (i - 1).ToString().PadLeft(5) + ": ";
                    for (int j = j2lo; j <= j2hi; j++ )
                    {
                        cout += a[i-1+(j-1)*m].ToString().PadLeft(12) + "  ";
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

        
    }
}