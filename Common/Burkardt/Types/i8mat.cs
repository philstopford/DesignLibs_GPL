using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void i8mat_transpose_print ( int m, int n, long[] a, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8MAT_TRANSPOSE_PRINT prints an I8MAT, transposed.
        //
        //  Discussion:
        //
        //    An I8MAT is an MxN array of I8's, stored by (I,J) -> [I+J*M].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 January 2005
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
        //    Input, long A[M*N], the M by N matrix.
        //
        //    Input, string TITLE, a title.
        //
        {
            i8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );
        }

        static void i8mat_transpose_print_some ( int m, int n, long[] a, int ilo, int jlo,
        int ihi, int jhi, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8MAT_TRANSPOSE_PRINT_SOME prints some of an I8MAT, transposed.
        //
        //  Discussion:
        //
        //    An I8MAT is an MxN array of I8's, stored by (I,J) -> [I+J*M].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 October 2014
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
        //    Input, long A[M*N], the matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title.
        //
        {
            int INCX = 10;

            Console.WriteLine();
            Console.WriteLine(title);

            if ( m <= 0 || n <= 0 )
            {
                Console.WriteLine();
                Console.WriteLine("  (None)");
                return;
            }
            //
            //  Print the columns of the matrix, in strips of INCX.
            //
            for (int i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
            {
                int i2hi = i2lo + INCX - 1;
                if ( m < i2hi )
                {
                    i2hi = m;
                }
                if ( ihi < i2hi )
                {
                    i2hi = ihi;
                }
                Console.WriteLine();
                //
                //  For each row I in the current range...
                //
                //  Write the header.
                //
                string line = "  Row: ";
                for (int  i = i2lo; i <= i2hi; i++ )
                {
                    string t = (i - 1) + "  ";
                    line += t.PadLeft(6);
                }
                Console.WriteLine(line);
                Console.WriteLine();
                Console.WriteLine("  Col");
                Console.WriteLine();
                //
                //  Determine the range of the rows in this strip.
                //
                int j2lo = jlo;
                if ( j2lo < 1 )
                {
                    j2lo = 1;
                }
                int j2hi = jhi;
                if ( n < j2hi )
                {
                    j2hi = n;
                }
                for (int  j = j2lo; j <= j2hi; j++ )
                {
                    //
                    //  Print out (up to INCX) entries in column J, that lie in the current strip.
                    //
                    string t = (j - 1) + ":";
                    line = "";
                    line += t.PadLeft(5);

                    for (int  i = i2lo; i <= i2hi; i++ )
                    {
                        t = a[i-1+(j-1)*m] + "  ";
                        line += t.PadLeft(6);
                    }
                    Console.WriteLine(line);

                    Console.WriteLine();
                }
            }
        }        
        
        public static void i8mat_print ( int m, int n, long[] a, string title )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8MAT_PRINT prints an I8MAT, with an optional title.
        //
        //  Discussion:
        //
        //    An I8MAT is an array of I8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 April 2003
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
        //    Input, long A[M*N], the M by N matrix.
        //
        //    Input, string TITLE, a title.
        //
        {
            i8mat_print_some ( m, n, a, 1, 1, m, n, title );
        }

        static void i8mat_print_some ( int m, int n, long[] a, int ilo, int jlo, int ihi,
        int jhi, string title )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8MAT_PRINT_SOME prints some of an I8MAT.
        //
        //  Discussion:
        //
        //    An I8MAT is an array of I8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 April 2004
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
        //    Input, long A[M*N], the matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title.
        {
            int INCX = 10;

            int i;
            int i2hi;
            int i2lo;
            int j;
            int j2hi;
            int j2lo;

            Console.WriteLine();
            Console.WriteLine(title);
            //
            //  Print the columns of the matrix, in strips of INCX.
            //
            for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
            {
                j2hi = j2lo + INCX - 1;
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
                string cout = "  Col: ";
                for ( j = j2lo; j <= j2hi; j++ )
                {
                    cout += j.ToString().PadLeft(6) + "  ";
                }
                Console.WriteLine(cout);
                Console.WriteLine("  Row");
                Console.WriteLine("  ---");
                //
                //  Determine the range of the rows in this strip.
                //
                if ( 1 < ilo )
                {
                    i2lo = 1;
                }
                else
                {
                    i2lo = ilo;
                }
                if ( ihi < m )
                {
                    i2hi = ihi;
                }
                else
                {
                    i2hi = m;
                }

                cout = "";
                for ( i = i2lo; i <= i2hi; i++ )
                {
                    //
                    //  Print out (up to INCX) entries in row I, that lie in the current strip.
                    //
                    cout += i.ToString().PadLeft(5) + "  ";
                    for ( j = j2lo; j <= j2hi; j++ )
                    {
                        cout += a[i-1+(j-1)*m].ToString().PadLeft(6) + "  ";
                    }
                    Console.WriteLine(cout);
                }
            }
        }
    }
}