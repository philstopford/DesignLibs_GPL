using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void i4mat_transpose_print ( int m, int n, int[] a, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
        //
        //  Discussion:
        //
        //    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
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
        //    Input, int A[M*N], the M by N matrix.
        //
        //    Input, string TITLE, a title.
        //
        {
            i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );
        }
        //****************************************************************************80

        public static void i4mat_transpose_print_some ( int m, int n, int[] a, int ilo, int jlo,
        int ihi, int jhi, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.
        //
        //  Discussion:
        //
        //    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
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
        //    Input, int A[M*N], the matrix.
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
        
        
        public static void i4mat_print ( int m, int n, int[] a, string title )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4MAT_PRINT prints an I4MAT, with an optional title.
        //
        //  Discussion:
        //
        //    An I4MAT is an array of I4's.
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
        //    Input, int A[M*N], the M by N matrix.
        //
        //    Input, string TITLE, a title.
        //
        {
            i4mat_print_some ( m, n, a, 1, 1, m, n, title );
        }

        public static void i4mat_print_some ( int m, int n, int[] a, int ilo, int jlo, int ihi,
        int jhi, string title )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4MAT_PRINT_SOME prints some of an I4MAT.
        //
        //  Discussion:
        //
        //    An I4MAT is an array of I4's.
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
        //    Input, int A[M*N], the matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title.
        {
            int INCX = 10;


            Console.WriteLine();
            Console.WriteLine(title);
            //
            //  Print the columns of the matrix, in strips of INCX.
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
                string cout = "  Col: ";
                for (int j = j2lo; j <= j2hi; j++ )
                {
                    cout += j.ToString().PadLeft(6) + "  ";
                }
                Console.WriteLine(cout);
                Console.WriteLine("  Row");
                Console.WriteLine("  ---");
                //
                //  Determine the range of the rows in this strip.
                //
                int i2lo;
                if ( 1 < ilo )
                {
                    i2lo = 1;
                }
                else
                {
                    i2lo = ilo;
                }

                int i2hi;
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
                    //  Print out (up to INCX) entries in row I, that lie in the current strip.
                    //
                    cout = i.ToString().PadLeft(5) + "  ";
                    for (int j = j2lo; j <= j2hi; j++ )
                    {
                        cout += a[i-1+(j-1)*m].ToString().PadLeft(6) + "  ";
                    }
                    Console.WriteLine(cout);
                }
            }
        }
        
        
        public static int[] i4mat_indicator_new ( int m, int n )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4MAT_INDICATOR_NEW sets up an "indicator" I4MAT.
        //
        //  Discussion:
        //
        //    An I4MAT is an array of I4's.
        //
        //    The value of each entry suggests its location, as in:
        //
        //      11  12  13  14
        //      21  22  23  24
        //      31  32  33  34
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 January 2005
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
        //    Output, int I4MAT_INDICATOR_NEW[M*N], the indicator matrix.
        //
        {
            int[] table = new int[m*n];

            int fac = ( int ) Math.Pow ( 10.0, Math.Floor( Math.Log10 ( n ) + 1 ) );

            for (int i = 1; i <= m; i++ )
            {
                for (int j = 1; j <= n; j++ )
                {
                    table[i-1+(j-1)*m] = fac * i + j;
                }
            }

            return table;
        }
    }
}