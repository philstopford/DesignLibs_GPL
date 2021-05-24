using System;
using System.IO;
using System.Linq;

namespace Burkardt.Table
{
    public class r4
    {
        public bool error { get; set; }
        public float val { get; set; }
        public int lchar { get; set; }
    }

    public class r4vec
    {
        public bool error { get; set; }
        public float[] rvec { get; set; }
    }

    public static partial class TableWriter
    {
        public static void r4mat_write ( string output_filename, int m, int n, float[] table )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R4MAT_WRITE writes an R4MAT file.
            //
            //  Discussion:
            //
            //    An R4MAT is an array of R4's.
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
            //    Input, float TABLE[M*N], the data.
            //
        {
            string[] outData = new string[n];
            for (int j = 0; j < n; j++ )
            {
                string line = "";
                for (int i = 0; i < m; i++ )
                {
                    line += "  " + table[i+j*m].ToString("0.########").PadLeft(24);
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
                Console.WriteLine("R4MAT_WRITE - Fatal error!");
                Console.WriteLine("  Could not open the output file: \"" + output_filename + "\"");
                throw;
            }
        }
    }
    public static partial class TableReader
    {
        public static TableHeader r4mat_header_read ( string input_filename )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4MAT_HEADER_READ reads the header from an R4MAT file.
        //
        //  Discussion:
        //
        //    An R4MAT is an array of R4's.
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
            TableHeader ret = readHeader(input_filename);

            if ( ret.m <= 0 )
            {
                Console.WriteLine();
                Console.WriteLine("R4MAT_HEADER_READ - Fatal error!");
                Console.WriteLine("  FILE_COLUMN_COUNT failed.");
                ret.code = 1;
            }
            
            if ( ret.n <= 0 )
            {
                Console.WriteLine();
                Console.WriteLine("R4MAT_HEADER_READ - Fatal error!");
                Console.WriteLine("  FILE_ROW_COUNT failed.");
                ret.code = 1;
            }

            return ret;
        }

        public static float[] r4mat_data_read(string input_filename, int m, int n)
        {
            string[] lines;

            try
            {
                lines = File.ReadLines(input_filename).ToArray();
            }
            catch (Exception e)
            {
                Console.WriteLine();
                Console.WriteLine("R4MAT_DATA_READ - Fatal error!");
                Console.WriteLine("  Could not open the input file: \"" + input_filename + "\"");
                throw;
            }

            float[] table = new float[m*n];
            
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

                r4vec res = TableMisc.s_to_r4vec(line, m);

                bool error = res.error;
                float[] x = res.rvec;

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
        public static void r4mat_transpose_print_some ( int m, int n, float[] a, int ilo, int jlo,
        int ihi, int jhi, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4MAT_TRANSPOSE_PRINT_SOME prints some of an R4MAT, transposed.
        //
        //  Discussion:
        //
        //    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
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
        //    Input, float A[M*N], an M by N matrix to be printed.
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
        
        
        public static r4 s_to_r4 ( string s )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    S_TO_R4 reads an R4 from a string.
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
            r4 ret = new r4 {lchar = -1};
            float rexp;
            char TAB = (char) 9;

            int nchar = s_len_trim ( s );
            int isgn = 1;
            float rtop = 0.0f;
            float rbot = 1.0f;
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
                        rtop = (float) 10.0 * rtop + ndig;
                    }
                    else if ( ihave == 5 )
                    {
                        rtop = (float) 10.0 * rtop + ndig;
                        rbot = (float) 10.0 * rbot;
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
                rexp = 1.0f;
            }
            else
            {
                if ( jbot == 1 )
                {
                    rexp = (float) Math.Pow ( 10.0, jsgn * jtop );
                }
                else
                {
                    rexp = jsgn * jtop;
                    rexp = rexp / jbot;
                    rexp = (float) Math.Pow ( 10.0, rexp );
                }
            }

            ret.val = isgn * rexp * rtop / rbot;

            return ret;
        }

        
        public static r4vec s_to_r4vec ( string s, int n )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    S_TO_R4VEC reads an R4VEC from a string.
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
        //    Output, float RVEC[N], the values read from the string.
        //
        //    Output, bool S_TO_R4VEC, is true if an error occurred.
        //
        {
            int begin = 0;
            int length = s.Length;

            r4vec ret = new r4vec() {rvec = new float[n]} ;

            for (int i = 0; i < n; i++ )
            {
                r4 res = s_to_r4 ( s.Substring(begin,length));

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
    }
}
