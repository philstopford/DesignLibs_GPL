using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Burkardt.Table
{
    
    public static partial class TableMisc
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

        public static int file_column_count(string filename)
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

        public static int file_row_count(string filename)
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

        public static int ch_to_digit ( char ch )
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

        

        

    }
}