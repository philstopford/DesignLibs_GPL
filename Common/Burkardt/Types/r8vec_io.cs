using System;
using System.Globalization;
using System.IO;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8vec_data_read ( string input_filename, int n, ref double[] table )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_DATA_READ reads the data from an R8VEC file.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The file is assumed to contain one record per line.
        //
        //    Records beginning with '#' are comments, and are ignored.
        //    Blank lines are also ignored.
        //
        //    There are assumed to be exactly (or at least) N such records.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string INPUT_FILENAME, the name of the input file.
        //
        //    Input, int N, the number of points.  The program
        //    will stop reading data once N values have been read.
        //
        //    Output, double TABLE[N], the data.
        //
    {
        try
        {
            string[] lines = File.ReadAllLines(input_filename);
            int j = 0;
            foreach (string line in lines)
            {
                if ( line[0] == '#' || s_len_trim ( line ) == 0 )
                {
                    continue;
                }

                table[j] = Convert.ToDouble( line );
                j += 1;
            }
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("R8VEC_DATA_READ - Fatal error!");
            Console.WriteLine("  Could not open the input file: \"" + input_filename + "\"");
        }
    }

    public static void r8vec_write(string output_filename, int n, double[] x)
    {
        string[] x_ = new string[n];
        for (int i = 0; i < n; i++)
        {
            x_[i] = x[i].ToString(CultureInfo.InvariantCulture);
        }
        File.WriteAllLines(output_filename, x_);
    }
        
}