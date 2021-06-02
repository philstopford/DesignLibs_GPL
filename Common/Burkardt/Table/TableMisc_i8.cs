using System;
using System.IO;
using System.Linq;
using Burkardt.Types;

namespace Burkardt.Table
{
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
                if (line[0] == '#' || typeMethods.s_len_trim(line) == 0)
                {
                    continue;
                }

                i8vec res = typeMethods.s_to_i8vec(line, m);

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

}