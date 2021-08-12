using System;
using System.Collections.Generic;
using System.IO;

namespace Burkardt.IO
{
    public static class DTable
    {
        public static void dtable_data_write(ref List<string> output, int m, int n, double[] table )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DTABLE_DATA_WRITE writes data to a DTABLE file.
        //
        //  Discussion:
        //
        //    The file should already be open.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 December 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, ofstream &OUTPUT, a pointer to the output stream.
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

            for (j = 0; j < n; j++)
            {
                string cout = "";
                for (i = 0; i < m; i++)
                {
                    cout += table[i + j * m].ToString().PadLeft(10) + "  ";
                }

                output.Add(cout);
            }
        }

        public static void dtable_write(string output_filename, int m, int n, double[] table,
        bool header )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DTABLE_WRITE writes information to a DTABLE file.
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
        //    Input, string OUTPUT_FILENAME, the output filename.
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double TABLE[M*N], the table data.
        //
        //    Input, bool HEADER, is TRUE if the header is to be included.
        //
        {

            try
            {
                List<string> output = new List<string>();
                if (header)
                {
                    //  dtable_header_write ( output_filename, output, m, n );
                }

                dtable_data_write(ref output, m, n, table);
            
                File.WriteAllLines(output_filename, output);

            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("DTABLE_WRITE - Fatal error!");
                Console.WriteLine("  Could not open the output file.");
            }
        }
    }
}