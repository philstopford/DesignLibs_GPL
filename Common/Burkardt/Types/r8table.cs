using System;
using System.IO;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] r8table_data_read ( string input_filename, int m, int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8TABLE_DATA_READ reads the data from an R8TABLE file.
            //
            //  Discussion:
            //
            //    The file is assumed to contain one record per line.
            //
            //    Records beginning with the '#' character are comments, and are ignored.
            //    Blank lines are also ignored.
            //
            //    Each line that is not ignored is assumed to contain exactly (or at least)
            //    M real numbers, representing the coordinates of a point.
            //
            //    There are assumed to be exactly (or at least) N such records.
            //
            //  Modified:
            //
            //    27 January 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, char *INPUT_FILENAME, the name of the input file.
            //
            //    Input, int M, the number of spatial dimensions.
            //
            //    Input, int N, the number of points.  The program
            //    will stop reading data once N values have been read.
            //
            //    Output, double DTABLE_DATA_READ[M*N], the table data.
            //
        {
            string[] input;
            int i;
            int j;
            string line;
            double[] table;

            try
            {
                input = File.ReadAllLines(input_filename);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("R8TABLE_DATA_READ - Fatal error!");
                Console.WriteLine("  Could not open the input file: \"" + input_filename + "\"");
                return null;
            }

            table = new double[m*n];
            
            j = 0;

            int index = 0;

            while ( j < n )
            {
                try
                {
                    line = input[index];
                    index++;
                }
                catch
                {
                    break;
                }

                if ( line[0] == '#' || s_len_trim ( line ) == 0 )
                {
                    continue;
                }

                r8vec res = s_to_r8vec ( line, m );
                

                if ( res.error )
                {
                    continue;
                }

                for ( i = 0; i < m; i++ )
                {
                    table[i+j*m] = res.rvec[i];
                }
                j = j + 1;

            }

            return table;
        }
    }
}