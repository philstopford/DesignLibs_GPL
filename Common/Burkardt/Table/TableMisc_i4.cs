using System;
using System.IO;
using System.Linq;
using Burkardt.Types;

namespace Burkardt.Table
{
    public static partial class TableWriter
    {
        public static void i4mat_write ( string output_filename, int m, int n, int[] table )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4MAT_WRITE writes an I4MAT file with no header.
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
        //    01 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string OUTPUT_FILENAME, the output filename.
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, int TABLE[M*N], the data.
        //
        {
            string[] outData = new string[n];
            for (int j = 0; j < n; j++ )
            {
                string line = "";
                for (int i = 0; i < m; i++ )
                {
                    line += "  " + table[i+j*m].ToString().PadLeft(10);
                }

                outData[j] = line;
            }

            try
            {
                File.WriteAllLines(output_filename, outData);
            }
            catch (Exception e)
            {
                Console.WriteLine();
                Console.WriteLine("I4MAT_WRITE - Fatal error!");
                Console.WriteLine("  Could not open the output file: \"" + output_filename + "\"");
                throw;
            }
        }
    }
    
    public static partial class TableReader
    {
        public static TableHeader i4mat_header_read ( string input_filename )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4MAT_HEADER_READ reads the header from an I4MAT file.
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
            //    23 February 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string INPUT_FILENAME, the name of the input file.
            //
            //    Output, int &M, the number of spatial dimensions.
            //
            //    Output, int &N, the number of points
            //
        {
            TableHeader ret = readHeader(input_filename);

            if ( ret.m <= 0 )
            {
                Console.WriteLine();
                Console.WriteLine("I4MAT_HEADER_READ - Fatal error!");
                Console.WriteLine("  FILE_COLUMN_COUNT failed.");
                ret.code = 1;
            }
            
            if ( ret.n <= 0 )
            {
                Console.WriteLine();
                Console.WriteLine("I4MAT_HEADER_READ - Fatal error!");
                Console.WriteLine("  FILE_ROW_COUNT failed.");
                ret.code = 1;
            }

            return ret;
        }
        
        
        public static int[] i4mat_data_read ( string input_filename, int m, int n )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4MAT_DATA_READ reads data from an I4MAT file.
        //
        //  Discussion:
        //
        //    An I4MAT is an array of I4's.
        //
        //    The file is assumed to contain one record per line.
        //
        //    Records beginning with '#' are comments, and are ignored.
        //    Blank lines are also ignored.
        //
        //    Each line that is not ignored is assumed to contain exactly (or at least)
        //    M real numbers, representing the coordinates of a point.
        //
        //    There are assumed to be exactly (or at least) N such records.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string INPUT_FILENAME, the name of the input file.
        //
        //    Input, int M, the number of spatial dimensions.
        //
        //    Input, int N, the number of points.  The program
        //    will stop reading data once N values have been read.
        //
        //    Output, int I4MAT_DATA_READ[M*N], the data.
        //
        {
            string[] lines;

            try
            {
                lines = File.ReadLines(input_filename).ToArray();
            }
            catch (Exception e)
            {
                Console.WriteLine();
                Console.WriteLine("R8MAT_DATA_READ - Fatal error!");
                Console.WriteLine("  Could not open the input file: \"" + input_filename + "\"");
                throw;
            }

            int[] table = new int[m*n];
            
            int j = 0;
            int l = 0;

            while ( j < n )
            {
                string line = lines[l];
                l++;
                if (line[0] == '#' || typeMethods.s_len_trim(line) == 0)
                {
                    continue;
                }

                i4vec res = typeMethods.s_to_i4vec(line, m);

                bool error = res.error;
                int[] x = res.ivec;

                if (!error)
                {
                    int i;
                    for (i = 0; i < m; i++)
                    {
                        table[i + (j * m)] = x[i];
                    }
                }
                j += 1;
            }
            
            return table;
        }
        
    }

    public static partial class TableMisc
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

        static void i4mat_transpose_print_some ( int m, int n, int[] a, int ilo, int jlo,
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