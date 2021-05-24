using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Burkardt.Types;

namespace Burkardt.Table
{
    public class TableHeader
    {
        public int m { get; set; }
        public int n { get; set; }
        public int code { get; set; }

        public TableHeader()
        {
            code = 1;
        }
    }

    public static class TableSummary
    {
        static void r8block_print ( int l, int m, int n, double[] a, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BLOCK_PRINT prints a double precision block (a 3D matrix).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int L, M, N, the dimensions of the block.
        //
        //    Input, double A[L*M*N], the matrix to be printed.
        //
        //    Input, char *TITLE, a title to be printed first.
        //    TITLE may be blank.
        //
        {
            if ( 0 < typeMethods.s_len_trim ( title ) )
            {
                Console.WriteLine();
                Console.WriteLine(title);
            }

            for (int k = 1; k <= n; k++ )
            {
                Console.WriteLine();
                Console.WriteLine("  K = " + k);
                Console.WriteLine();
                for (int jlo = 1; jlo <= m; jlo = jlo + 5 )
                {
                    int jhi = Math.Min( jlo + 4, m );
                    Console.WriteLine();
                    string cout = "      ";
                    for (int j = jlo; j <= jhi; j++ )
                    {
                        cout += j.ToString().PadLeft(7) + "       ";
                    }
                    Console.WriteLine(cout);
                    Console.WriteLine();
                    cout = "";
                    for (int i = 1; i <= l; i++ )
                    {
                        string t = i.ToString().PadLeft(4);
                        cout += "  " + t;
                        for (int j = jlo; j <= jhi; j++ )
                        {
                            t = a[i-1+(j-1)*l+(k-1)*l*m].ToString().PadLeft(12);
                            cout += "  " +t;
                        }
                        Console.WriteLine(cout);
                    }
                }
            }
        }


        static void r8mat_print ( int m, int n, double[] a, string title )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_PRINT prints an R8MAT, with an optional title.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector
        //    in column-major order.
        //
        //    Entry A(I,J) is stored as A[I+J*M]
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
        //    Input, char *TITLE, a title to be printed.
        //
        {
            r8mat_print_some ( m, n, a, 1, 1, m, n, title );
        }


        static void r8mat_print_some ( int m, int n, double[] a, int ilo, int jlo, int ihi,
  int jhi, string title )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_PRINT_SOME prints some of an R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector
        //    in column-major order.
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
        //    Input, double A[M*N], the matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, char *TITLE, a title for the matrix.
        //
        {
            int INCX = 5;
            
            if ( 0 < typeMethods.s_len_trim ( title ) )
            {
                Console.WriteLine();
                Console.WriteLine(title);
            }
            //
            //  Print the columns of the matrix, in strips of 5.
            //
            for (int j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
            {
                int j2hi = j2lo + INCX - 1;
                j2hi = Math.Min( j2hi, n );
                j2hi = Math.Min( j2hi, jhi );

                Console.WriteLine();
                //
                //  For each column J in the current range...
                //
                //  Write the header.
                //
                string cout = "  Col:    ";
                for (int  j = j2lo; j <= j2hi; j++ )
                {
                    string t = j.ToString().PadLeft(7) + "       ";
                    cout += t;
                }
                Console.WriteLine(cout);
                Console.WriteLine("  Row");
                Console.WriteLine();
                //
                //  Determine the range of the rows in this strip.
                //
                int i2lo = Math.Max ( ilo, 1 );
                int i2hi = Math.Min ( ihi, m );

                cout = "";
                for (int  i = i2lo; i <= i2hi; i++ )
                {
                    //
                    //  Print out (up to) 5 entries in row I, that lie in the current strip.
                    //
                    cout += i.ToString().PadLeft(5) + "  ";
                    for (int  j = j2lo; j <= j2hi; j++ )
                    {
                        cout += a[i - 1 + (j - 1) * m].ToString().PadLeft(12) + "  ";
                    }
                    Console.WriteLine(cout);
                }
            }

            Console.WriteLine();
        }

        static void r8vec_print ( int n, double[] a, string title )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_PRINT prints an R8VEC.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of components of the vector.
            //
            //    Input, double A[N], the vector to be printed.
            //
            //    Input, string TITLE, a title.
            //
        {
            Console.WriteLine();
            Console.WriteLine(title);
            Console.WriteLine();
            for (int i = 0; i < n; i++ )
            {
                string cout = "  ";
                
                string t = i.ToString().PadLeft(8);

                string t2 = a[i].ToString().PadLeft(14);
                
                cout += t + ": " + t2;  
                
                Console.WriteLine(cout);;
            }
        }

    }
}