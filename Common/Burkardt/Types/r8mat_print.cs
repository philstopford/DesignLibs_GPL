using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8mat_l_print(int m, int n, double[] a, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_L_PRINT prints a lower triangular R8MAT.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //  Example:
            //
            //    M = 5, N = 5
            //    A = (/ 11, 21, 31, 41, 51, 22, 32, 42, 52, 33, 43, 53, 44, 54, 55 /)
            //
            //    11
            //    21 22
            //    31 32 33
            //    41 42 43 44
            //    51 52 53 54 55
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 January 2016
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
            //    Input, double A[*], the M by N matrix.  Only the lower
            //    triangular elements are stored, in column major order.
            //
            //    Input, string TITLE, a title.
            //
        {
            int i;
            int[] indx = new int[10];
            int j;
            int jhi;
            int jlo;
            int jmax;
            int nn;
            int size;

            Console.WriteLine("");
            Console.WriteLine(title + "");

            jmax = Math.Min(n, m);

            if (m <= n)
            {
                size = (m * (m + 1)) / 2;
            }
            else
            {
                size = (n * (n + 1)) / 2 + (m - n) * n;
            }

            if (r8vec_is_integer(size, a))
            {
                nn = 10;
                for (jlo = 1; jlo <= jmax; jlo = jlo + nn)
                {
                    jhi = Math.Min(jlo + nn - 1, Math.Min(m, jmax));
                    Console.WriteLine("");
                    string cout = "  Col   ";
                    for (j = jlo; j <= jhi; j++)
                    {
                        cout += j.ToString().PadLeft(6);
                    }

                    Console.WriteLine(cout);
                    Console.WriteLine("  Row  ");
                    for (i = jlo; i <= m; i++)
                    {
                        jhi = Math.Min(jlo + nn - 1, Math.Min(i, jmax));
                        for (j = jlo; j <= jhi; j++)
                        {
                            indx[j - jlo] = (j - 1) * m + i - (j * (j - 1)) / 2;
                        }

                        cout = "  " + i.ToString().PadLeft(6);
                        for (j = 0; j <= jhi - jlo; j++)
                        {
                            cout += a[indx[j] - 1].ToString().PadLeft(6);
                        }

                        Console.WriteLine(cout);
                    }
                }
            }
            else if (r8vec_amax(size, a) < 1000000.0)
            {
                nn = 5;
                for (jlo = 1; jlo <= jmax; jlo = jlo + nn)
                {
                    jhi = Math.Min(jlo + nn - 1, Math.Min(m - 1, jmax));
                    Console.WriteLine("");
                    string cout = "  Col   ";
                    for (j = jlo; j <= jhi; j++)
                    {
                        cout += j.ToString().PadLeft(14);
                    }

                    Console.WriteLine(cout);
                    Console.WriteLine("  Row  ");
                    for (i = jlo; i <= m; i++)
                    {
                        jhi = Math.Min(jlo + nn - 1, Math.Min(i, jmax));
                        for (j = jlo; j <= jhi; j++)
                        {
                            indx[j - jlo] = (j - 1) * m + i - (j * (j - 1)) / 2;
                        }

                        cout += "  " + i.ToString().PadLeft(6);
                        for (j = 0; j <= jhi - jlo; j++)
                        {
                            cout += a[indx[j] - 1].ToString().PadLeft(14);
                        }

                        Console.WriteLine(cout);
                    }
                }
            }
            else
            {
                nn = 5;

                for (jlo = 1; jlo <= jmax; jlo = jlo + nn)
                {
                    jhi = Math.Min(jlo + nn - 1, Math.Min(m - 1, jmax));
                    Console.WriteLine("");
                    string cout = "  Col ";
                    for (j = jlo; j <= jhi; j++)
                    {
                        cout += j.ToString().PadLeft(7) + "       ";
                    }

                    Console.WriteLine(cout);
                    Console.WriteLine("  Row ");
                    for (i = jlo; i <= m; i++)
                    {
                        jhi = Math.Min(jlo + nn - 1, Math.Min(i, jmax));
                        for (j = jlo; j <= jhi; j++)
                        {
                            indx[j - jlo] = (j - 1) * m + i - (j * (j - 1)) / 2;
                        }

                        cout = i.ToString().PadLeft(6);
                        for (j = 0; j <= jhi - jlo; j++)
                        {
                            cout += a[indx[j] - 1].ToString().PadLeft(14);
                        }

                        Console.WriteLine(cout);
                    }
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
        
        
    }
}