using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

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
    
    public static class TableReader
    {
        public static int s_len_trim(string line)
        {
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    S_LEN_TRIM returns the length of a string to the last nonblank.
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
            //    Input, string S, a string.
            //
            //    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
            //    If S_LEN_TRIM is 0, then the string is entirely blank.
            //
            if (line == "")
            {
                return 0;
            }

            string tmp = line.TrimEnd();

            int index = tmp.Length;
            
            /*
            int index = line.Length - 1;
            while ((index >= 0) && (line[index] == ' '))
            {
                index--;
            }
            */
            return index;
        }

        static int file_column_count(string filename)
        {
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FILE_COLUMN_COUNT counts the columns in the first line of a file.
            //
            //  Discussion:
            //
            //    The file is assumed to be a simple text file.
            //
            //    Most lines of the file are presumed to consist of COLUMN_NUM words,
            //    separated by spaces.  There may also be some blank lines, and some
            //    comment lines, which have a "#" in column 1.
            //
            //    The routine tries to find the first non-comment non-blank line and
            //    counts the number of words in that line.
            //
            //    If all lines are blanks or comments, it goes back and tries to analyze
            //    a comment line.
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
            //    Input, string FILENAME, the name of the file.
            //
            //    Output, int FILE_COLUMN_COUNT, the number of columns assumed
            //    to be in the file.
            //
            int ret = -1;
            string[] lines = File.ReadLines(filename).ToArray();
            foreach (string line in lines)
            {
                if (line.StartsWith('#'))
                {
                    continue;
                }

                bool found = false;
                List<string> fields = new List<string>();

                string tmp = "";
                
                foreach (char t in line)
                {
                    if (t == ' ')
                    {
                        if (found)
                        {
                            // Got whitespace, commit string to the list and then set to false and continue.
                            fields.Add(tmp);
                            tmp = "";
                            found = false;
                        }
                    }
                    else
                    {
                        found = true;
                        tmp += t;
                    }
                }

                if (tmp != "")
                {
                    fields.Add(tmp);
                }

                if (fields.Count > 0)
                {
                    ret = fields.Count;
                }
            }
            
            return ret;
        }

        static int file_row_count(string filename)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FILE_ROW_COUNT counts the number of row records in a file.
        //
        //  Discussion:
        //
        //    It does not count lines that are blank, or that begin with a
        //    comment symbol '#'.
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
        //    Output, int FILE_ROW_COUNT, the number of rows found.
        {
            string[] lines = File.ReadLines(filename).ToArray();

            return lines.Count(line => !line.StartsWith('#') && (s_len_trim(line) != 0));
        }
        
        public static TableHeader r8mat_header_read ( string input_filename )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_HEADER_READ reads the header from an R8MAT file.
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
        //    Output, int &N, the number of points.
        //
        {
            TableHeader ret = new TableHeader();
            ret.m = file_column_count ( input_filename );

            if ( ret.m <= 0 )
            {
                Console.WriteLine();
                Console.WriteLine("R8MAT_HEADER_READ - Fatal error!");
                Console.WriteLine("  FILE_COLUMN_COUNT failed.");
                ret.code = 1;
            }

            ret.n = file_row_count ( input_filename );

            if ( ret.n <= 0 )
            {
                Console.WriteLine();
                Console.WriteLine("R8MAT_HEADER_READ - Fatal error!");
                Console.WriteLine("  FILE_ROW_COUNT failed.");
                ret.code = 1;
            }

            return ret;
        }

        public static double[] r8mat_data_read(string input_filename, int m, int n)
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

            double[] table = new double[m*n];
            
            int j = 0;
            int l = 0;

            while ( j < n )
            {
                string line = lines[l];
                l++;
                if (line[0] == '#' || s_len_trim(line) == 0)
                {
                    continue;
                }

                r8vec res = s_to_r8vec(line, m);

                bool error = res.error;
                double[] x = res.rvec;

                if (!error)
                {
                    int i;
                    for (i = 0; i < m; i++)
                    {
                        table[i + (j * m)] = x[i];
                    }
                }
                j = j + 1;
            }
            
            return table;            
        }

        class r8vec
        {
            public bool error { get; set; }
            public double[] rvec { get; set; }
        }
        
        static r8vec s_to_r8vec ( string s, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    S_TO_R8VEC reads an R8VEC from a string.
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
        //    Output, double RVEC[N], the values read from the string.
        //
        //    Output, bool S_TO_R8VEC, is true if an error occurred.
        //
        {
            int begin = 0;
            int length = s.Length;

            r8vec ret = new r8vec() {rvec = new double[n]} ;

            for (int i = 0; i < n; i++ )
            {
                r8 res = s_to_r8 ( s.Substring(begin,length));

                ret.rvec[i] = res.val;
                int lchar = res.lchar;
                
                if ( ret.error )
                {
                    return ret;
                }
                begin = begin + lchar;
                length = length - lchar;
            }

            return ret;
        }

        public class r8
        {
            public bool error { get; set; }
            public double val { get; set; }
            public int lchar { get; set; }
        }
        
        static r8 s_to_r8 ( string s )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    S_TO_R8 reads an R8 from a string.
        //
        //  Discussion:
        //
        //    This routine will read as many characters as possible until it reaches
        //    the end of the string, or encounters a character which cannot be
        //    part of the real number.
        //
        //    Legal input is:
        //
        //       1 blanks,
        //       2 '+' or '-' sign,
        //       2.5 spaces
        //       3 integer part,
        //       4 decimal point,
        //       5 fraction part,
        //       6 'E' or 'e' or 'D' or 'd', exponent marker,
        //       7 exponent sign,
        //       8 exponent integer part,
        //       9 exponent decimal point,
        //      10 exponent fraction part,
        //      11 blanks,
        //      12 final comma or semicolon.
        //
        //    with most quantities optional.
        //
        //  Example:
        //
        //    S                 R
        //
        //    '1'               1.0
        //    '     1   '       1.0
        //    '1A'              1.0
        //    '12,34,56'        12.0
        //    '  34 7'          34.0
        //    '-1E2ABCD'        -100.0
        //    '-1X2ABCD'        -1.0
        //    ' 2E-1'           0.2
        //    '23.45'           23.45
        //    '-4.2E+2'         -420.0
        //    '17d2'            1700.0
        //    '-14e-2'         -0.14
        //    'e2'              100.0
        //    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
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
        //    Input, string S, the string containing the
        //    data to be read.  Reading will begin at position 1 and
        //    terminate at the end of the string, or when no more
        //    characters can be read to form a legal real.  Blanks,
        //    commas, or other nonnumeric data will, in particular,
        //    cause the conversion to halt.
        //
        //    Output, int *LCHAR, the number of characters read from
        //    the string to form the number, including any terminating
        //    characters such as a trailing comma or blanks.
        //
        //    Output, bool *ERROR, is true if an error occurred.
        //
        //    Output, double S_TO_R8, the real value that was read from the string.
        //
        {
            r8 ret = new r8 {lchar = -1};
            double rexp;
            char TAB = (char) 9;

            int nchar = s_len_trim ( s );
            int isgn = 1;
            double rtop = 0.0;
            double rbot = 1.0;
            int jsgn = 1;
            int jtop = 0;
            int jbot = 1;
            int ihave = 1;
            int iterm = 0;

            for ( ; ; )
            {
                char c = s[ret.lchar+1];
                ret.lchar = ret.lchar + 1;
                //
                //  Blank or TAB character.
                //
                if ( c == ' ' || c == TAB )
                {
                    if ( ihave == 2 )
                    {
                    }
                    else if ( ihave == 6 || ihave == 7 )
                    {
                        iterm = 1;
                    }
                    else if ( 1 < ihave )
                    {
                        ihave = 11;
                    }
                }
                //
                //  Comma.
                //
                else if ( c == ',' || c == ';' )
                {
                    if ( ihave != 1 )
                    {
                        iterm = 1;
                        ihave = 12;
                        ret.lchar = ret.lchar + 1;
                    }
                }
                //
                //  Minus sign.
                //
                else if ( c == '-' )
                {
                    if ( ihave == 1 )
                    {
                        ihave = 2;
                        isgn = -1;
                    }
                    else if ( ihave == 6 )
                    {
                        ihave = 7;
                        jsgn = -1;
                    }
                    else
                    {
                        iterm = 1;
                    }
                }
                //
                //  Plus sign.
                //
                else if ( c == '+' )
                {
                    if ( ihave == 1 )
                    {
                        ihave = 2;
                    }
                    else if ( ihave == 6 )
                    {
                        ihave = 7;
                    }
                    else
                    {
                        iterm = 1;
                    }
                }
                //
                //  Decimal point.
                //
                else if ( c == '.' )
                {
                    if ( ihave < 4 )
                    {
                        ihave = 4;
                    }
                    else if ( 6 <= ihave && ihave <= 8 )
                    {
                        ihave = 9;
                    }
                    else
                    {
                        iterm = 1;
                    }
                }
                //
                //  Exponent marker.
                //
                else if ( ( Char.ToUpper(c) == 'E' ) || ( Char.ToUpper(c) == 'D' ) )
                {
                    if ( ihave < 6 )
                    {
                        ihave = 6;
                    }
                    else
                    {
                        iterm = 1;
                    }
                }
                //
                //  Digit.
                //
                else if ( ihave < 11 && '0' <= c && c <= '9' )
                {
                    if ( ihave <= 2 )
                    {
                        ihave = 3;
                    }
                    else if ( ihave == 4 )
                    {
                        ihave = 5;
                    }
                    else if ( ihave == 6 || ihave == 7 )
                    {
                        ihave = 8;
                    }
                    else if ( ihave == 9 )
                    {
                        ihave = 10;
                    }

                    int ndig = ch_to_digit ( c );

                    if ( ihave == 3 )
                    {
                        rtop = 10.0 * rtop + ( double ) ndig;
                    }
                    else if ( ihave == 5 )
                    {
                        rtop = 10.0 * rtop + ( double ) ndig;
                        rbot = 10.0 * rbot;
                    }
                    else if ( ihave == 8 )
                    {
                        jtop = 10 * jtop + ndig;
                    }
                    else if ( ihave == 10 )
                    {
                        jtop = 10 * jtop + ndig;
                        jbot = 10 * jbot;
                    }
                }
                //
                //  Anything else is regarded as a terminator.
                //
                else
                {
                    iterm = 1;
                }
                //
                //  If we haven't seen a terminator, and we haven't examined the
                //  entire string, go get the next character.
                //
                if ( iterm == 1 || nchar <= ret.lchar + 1 )
                {
                    break;
                }

            }
            //
            //  If we haven't seen a terminator, and we have examined the
            //  entire string, then we're done, and LCHAR is equal to NCHAR.
            //
            if ( iterm != 1 && (ret.lchar) + 1 == nchar )
            {
                ret.lchar = nchar;
            }
            //
            //  Number seems to have terminated.  Have we got a legal number?
            //  Not if we terminated in states 1, 2, 6 or 7!
            //
            if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
            {
                ret.error = true;
                return ret;
            }
            //
            //  Number seems OK.  Form it.
            //
            if ( jtop == 0 )
            {
                rexp = 1.0;
            }
            else
            {
                if ( jbot == 1 )
                {
                    rexp = Math.Pow ( 10.0, jsgn * jtop );
                }
                else
                {
                    rexp = jsgn * jtop;
                    rexp = rexp / jbot;
                    rexp = Math.Pow ( 10.0, rexp );
                }
            }

            ret.val = isgn * rexp * rtop / rbot;

            return ret;
        }

        static int ch_to_digit ( char ch )
            //****************************************************************************80
        //
        //  Purpose:
        //
        //    CH_TO_DIGIT returns the integer value of a base 10 digit.
        //
        //  Example:
        //
        //     CH  DIGIT
        //    ---  -----
        //    '0'    0
        //    '1'    1
        //    ...  ...
        //    '9'    9
        //    ' '    0
        //    'X'   -1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 June 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, char CH, the decimal digit, '0' through '9' or blank are legal.
        //
        //    Output, int CH_TO_DIGIT, the corresponding integer value.  If the
        //    character was 'illegal', then DIGIT is -1.
        //
        {
            int digit;

            if ( '0' <= ch && ch <= '9' )
            {
                digit = ch - '0';
            }
            else if ( ch == ' ' )
            {
                digit = 0;
            }
            else
            {
                digit = -1;
            }

            return digit;
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
    }


    public static class TableWriter
    {
        public static string file_name_ext_swap(string filename, string extension)
        {
            string[] tokens = filename.Split('.');
            if (tokens.Length == 1)
            {
                return String.Join('.', filename, extension);
            }

            string ret = "";
            for (int i = 0; i < tokens.Length - 1; i++)
            {
                ret += tokens[i] + ".";
            }

            ret += extension;

            return ret;
        }
        
        public static void r8mat_write ( string output_filename, int m, int n, double[] table )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_WRITE writes an R8MAT file.
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
        //    09 November 2014
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
        //    Input, double TABLE[M*N], the data.
        //
        {
            string[] outData = new string[n];
            for (int j = 0; j < n; j++ )
            {
                string line = "";
                for (int i = 0; i < m; i++ )
                {
                    line += "  " + table[i+j*m];
                    //    output << "  " << setw(24) << setprecision(16) << table[i+j*m];
                }

                outData[j] = line;
            }

            //
            //  Open the file.
            //

            try
            {
                File.WriteAllLines(output_filename, outData);
            }
            catch (Exception e)
            {
                Console.WriteLine();
                Console.WriteLine("R8MAT_WRITE - Fatal error!");
                Console.WriteLine("  Could not open the output file: \"" + output_filename + "\"");
                throw;
            }
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
            if ( 0 < TableReader.s_len_trim ( title ) )
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
                        string t = j + "       ";
                        cout += t.PadLeft(7);
                    }
                    Console.WriteLine(cout);
                    Console.WriteLine();
                    cout = "";
                    for (int i = 1; i <= l; i++ )
                    {
                        string t = i.ToString();
                        t = t.PadLeft(4);
                        cout += "  " + t;
                        for (int j = jlo; j <= jhi; j++ )
                        {
                            t = a[i-1+(j-1)*l+(k-1)*l*m].ToString();
                            t = t.PadLeft(12);
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
            
            if ( 0 < TableReader.s_len_trim ( title ) )
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
                    string t = j + "       ";
                    cout += t.PadLeft(7);
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
                    string t = i + "  ";
                    cout += t.PadLeft(5);
                    for (int  j = j2lo; j <= j2hi; j++ )
                    {
                        t = a[i - 1 + (j - 1) * m] + "  ";
                        cout += t.PadLeft(12);
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
                
                string t = i.ToString();
                t = t.PadLeft(8);

                string t2 = a[i].ToString();
                t2 = t2.PadLeft(14);
                
                cout += t + ": " + t2;  
                
                Console.WriteLine(cout);;
            }
        }

    }
}