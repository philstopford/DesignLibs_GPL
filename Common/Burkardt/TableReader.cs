using System;
using System.IO;
using System.Linq;

namespace Burkardt.Table
{
    public static class TableReader
    {
        
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
            ret.m = TableMisc.file_column_count ( input_filename );

            if ( ret.m <= 0 )
            {
                Console.WriteLine();
                Console.WriteLine("R8MAT_HEADER_READ - Fatal error!");
                Console.WriteLine("  FILE_COLUMN_COUNT failed.");
                ret.code = 1;
            }

            ret.n = TableMisc.file_row_count ( input_filename );

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
                if (line[0] == '#' || TableMisc.s_len_trim(line) == 0)
                {
                    continue;
                }

                r8vec res = TableMisc.s_to_r8vec(line, m);

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
                j += 1;
            }
            
            return table;            
        }


        

    }
}