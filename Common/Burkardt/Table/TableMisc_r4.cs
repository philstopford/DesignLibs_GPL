using System;
using System.IO;
using System.Linq;
using Burkardt.Types;

namespace Burkardt.Table
{ public static partial class TableWriter
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
                if (line[0] == '#' || typeMethods.s_len_trim(line) == 0)
                {
                    continue;
                }

                r4vec res = typeMethods.s_to_r4vec(line, m);

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

}
