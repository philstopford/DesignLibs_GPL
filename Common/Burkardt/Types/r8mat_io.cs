using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Burkardt.Table;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
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

            double[] table = new double[m * n];

            int j = 0;
            int l = 0;

            while (j < n)
            {
                string line = lines[l];
                l++;
                if (line[0] == '#' || typeMethods.s_len_trim(line) == 0)
                {
                    continue;
                }

                r8vec res = typeMethods.s_to_r8vec(line, m);

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

        public static void r8mat_write(string output_filename, int m, int n, double[] table)

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
            for (int j = 0; j < n; j++)
            {
                string line = "";
                for (int i = 0; i < m; i++)
                {
                    line += "  " + table[((i + j * m) + table.Length) % table.Length].ToString("0.################").PadLeft(24);
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
                Console.WriteLine("R8MAT_WRITE - Fatal error!");
                Console.WriteLine("  Could not open the output file: \"" + output_filename + "\"");
                throw;
            }
        }
        
        public static void r8mat_write0 ( string output_filename, int m, int n, double[] table )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_WRITE0 writes an R8MAT file with no header.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 June 2009
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
            //    Input, double TABLE[M*N], the table data.
            //
        {
            int i;
            int j;
            List<string> output = new List<string>();

            //
            //  Write the data.
            //  For greater precision, try
            //
            //    output << "  " << setw(24) << setprecision(16) << table[i+j*m];
            //
            for ( j = 0; j < n; j++ )
            {
                string cout = "";
                for ( i = 0; i < m; i++ )
                {
                    cout += "  " + table[i+j*m].ToString().PadLeft(10);
                }
                output.Add(cout);
            }

            try
            {
                File.WriteAllLines(output_filename, output);
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("R8MAT_WRITE0 - Fatal error!");
                Console.WriteLine("  Could not open the output file.");
            }
        }

        public static TableHeader r8mat_header_read(string input_filename)

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
            TableHeader ret = TableMisc.readHeader(input_filename);

            if (ret.m <= 0)
            {
                Console.WriteLine();
                Console.WriteLine("R8MAT_HEADER_READ - Fatal error!");
                Console.WriteLine("  FILE_COLUMN_COUNT failed.");
                ret.code = 1;
            }

            if (ret.n <= 0)
            {
                Console.WriteLine();
                Console.WriteLine("R8MAT_HEADER_READ - Fatal error!");
                Console.WriteLine("  FILE_ROW_COUNT failed.");
                ret.code = 1;
            }

            return ret;
        }
        
    }
}