using System;
using System.IO;

namespace Burkardt.Table
{
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
                    line += "  " + table[i+j*m].ToString("0.################").PadLeft(24);
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
}