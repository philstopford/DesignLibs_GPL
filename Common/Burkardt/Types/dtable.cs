using System;
using System.Collections.Generic;
using System.IO;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void dtable_write0 ( string output_filename, int m, int n, double[] table )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DTABLE_WRITE0 writes a DTABLE file with no header.
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

            try
            {
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
                        cout += "  " + table[i+j*m].ToString("0.################").PadLeft(24);
                    }
                    output.Add(cout);
                }
                //
                //  Close the file.
                //
                File.WriteAllLines(output_filename, output);
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("DTABLE_WRITE0 - Fatal error!");
                Console.WriteLine("  Could not open the output file.");
            }
        }
    }
}