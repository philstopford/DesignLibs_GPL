using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Burkardt.Types;

namespace Burkardt.Table
{
    
    public static partial class TableMisc
    {
        public static TableHeader readHeader(string input_filename)
        {
            TableHeader ret = new TableHeader {m = TableMisc.file_column_count(input_filename), n = TableMisc.file_row_count ( input_filename )};

            return ret;
        }
        
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

            return lines.Count(line => !line.StartsWith('#') && (typeMethods.s_len_trim(line) != 0));
        }


        

        

    }
}