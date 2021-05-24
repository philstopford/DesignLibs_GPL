using System;
using System.IO;
using System.Linq;

namespace Burkardt.Table
{
    public class i8
    {
        public bool error { get; set; }
        public long val { get; set; }
        public int lchar { get; set; }
    }

    public class i8vec
    {
        public bool error { get; set; }
        public long[] ivec { get; set; }
    }

    public static partial class TableWriter
    {
        public static void i8mat_write ( string output_filename, int m, int n, long[] table )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8MAT_WRITE writes an I8MAT file with no header.
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
        //    Input, long TABLE[M*N], the data.
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
                Console.WriteLine("I8MAT_WRITE - Fatal error!");
                Console.WriteLine("  Could not open the output file: \"" + output_filename + "\"");
                throw;
            }
        }
    }
    
    public static partial class TableReader
    {
        public static TableHeader i8mat_header_read ( string input_filename, int m, int n )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I8MAT_HEADER_READ reads the header from an I8MAT file.
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
                Console.WriteLine("I8MAT_HEADER_READ - Fatal error!");
                Console.WriteLine("  FILE_COLUMN_COUNT failed.");
                ret.code = 1;
            }
            
            if ( ret.n <= 0 )
            {
                Console.WriteLine();
                Console.WriteLine("I8MAT_HEADER_READ - Fatal error!");
                Console.WriteLine("  FILE_ROW_COUNT failed.");
                ret.code = 1;
            }

            return ret;
        }
        
        
        public static long[] i8mat_data_read ( string input_filename, int m, int n )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8MAT_DATA_READ reads data from an I8MAT file.
        //
        //  Discussion:
        //
        //    An I8MAT is an array of I8's.
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
        //    Output, long I8MAT_DATA_READ[M*N], the data.
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
                Console.WriteLine("I8MAT_DATA_READ - Fatal error!");
                Console.WriteLine("  Could not open the input file: \"" + input_filename + "\"");
                throw;
            }

            long[] table = new long[m*n];
            
            int j = 0;
            int l = 0;

            while ( j < n )
            {
                string line = lines[l];
                l++;
                if (line[0] == '#' || TableMisc.s_len_trim(line) == 0)
                {
                    continue;
                }

                i8vec res = TableMisc.s_to_i8vec(line, m);

                bool error = res.error;
                long[] x = res.ivec;

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

        
        public static i8 s_to_i8 ( string s )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    S_TO_I8 reads an I8 from a string.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string S, a string to be examined.
        //
        //    Output, int *LAST, the last character of S used to make IVAL.
        //
        //    Output, bool *ERROR is TRUE if an error occurred.
        //
        //    Output, int *S_TO_I4, the integer value read from the string.
        //    If the string is blank, then IVAL will be returned 0.
        //
        {
            i8 ival = new i8() {val = 0};

            int istate = 0;
            long isgn = 1;
            int i = 0;

            for ( ; ; )
            {
                char c = s[i];
                i = i + 1;
                //
                //  Haven't read anything.
                //
                if ( istate == 0 )
                {
                    if ( c == ' ' )
                    {
                    }
                    else if ( c == '-' )
                    {
                        istate = 1;
                        isgn = -1;
                    }
                    else if ( c == '+' )
                    {
                        istate = 1;
                        isgn = + 1;
                    }
                    else if ( '0' <= c && c <= '9' )
                    {
                        istate = 2;
                        ival.val = c - '0';
                    }
                    else
                    {
                        ival.error = true;
                        return ival;
                    }
                }
                //
                //  Have read the sign, expecting digits.
                //
                else if ( istate == 1 )
                {
                    if ( c == ' ' )
                    {
                    }
                    else if ( '0' <= c && c <= '9' )
                    {
                        istate = 2;
                        ival.val = c - '0';
                    }
                    else
                    {
                        ival.error = true;
                        return ival;
                    }
                }
                //
                //  Have read at least one digit, expecting more.
                //
                else if ( istate == 2 )
                {
                    if ( '0' <= c && c <= '9' )
                    {
                        ival.val = 10 * (ival.val) + c - '0';
                    }
                    else
                    {
                        ival.val = isgn * ival.val;
                        ival.lchar = i - 1;
                        return ival;
                    }
                }
            }
            //
            //  If we read all the characters in the string, see if we're OK.
            //
            if ( istate == 2 )
            {
                ival.val = isgn * ival.val;
                ival.lchar = s_len_trim ( s );
            }
            else
            {
                ival.error = true;
                ival.lchar = 0;
            }

            return ival;
        }
     
        
        public static i8vec s_to_i8vec ( string s, int n )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    S_TO_I4VEC reads an I4VEC from a string.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 July 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string S, the string to be read.
            //
            //    Input, int N, the number of values expected.
            //
            //    Output, long IVEC[N], the values read from the string.
            //
            //    Output, bool S_TO_I4VEC, is TRUE if an error occurred.
            //
        {
            i8vec ret = new i8vec() {ivec = new long[n]};

            int begin = 0;
            int length = s.Length;

            for ( int i = 0; i < n; i++ )
            {
                i8 res = s_to_i8 ( s.Substring(begin,length) );

                ret.ivec[i] = res.val;
                int lchar = res.lchar;
                
                if ( res.error )
                {
                    return ret;
                }
                begin = begin + lchar;
                length = length - lchar;
            }

            return ret;
        }

        
    }
}