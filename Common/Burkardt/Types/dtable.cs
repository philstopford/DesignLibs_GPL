using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Table;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] dtable_data_read(string input_filename, int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DTABLE_DATA_READ reads the data from a real TABLE file.
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
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2003
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

        try
        {
            input = File.ReadAllLines(input_filename);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("DTABLE_DATA_READ - Fatal error!");
            Console.WriteLine("  Could not open the input file: \"" + input_filename + "\"");
            return null;
        }

        double[] table = new double[m * n];

        double[] x = new double[m];

        int j = 0;
        int index = 0;

        while (j < n)
        {
            string line;
            try
            {
                line = input[index];
                index++;
            }
            catch (Exception)
            {
                break;
            }

            if (line[0] == '#' || s_len_trim(line) == 0)
            {
                continue;
            }

            r8vec res = s_to_r8vec(line, m);

            switch (res.error)
            {
                case true:
                    continue;
            }

            x = res.rvec;

            int i;
            for (i = 0; i < m; i++)
            {
                table[i + j * m] = x[i];
            }

            j += 1;

        }

        return table;
    }

    public static TableHeader dtable_header_read(string input_filename)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DTABLE_HEADER_READ reads the header from a real TABLE file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 June 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, char *INPUT_FILENAME, the name of the input file.
        //
        //    Output, int *M, the number of spatial dimensions.
        //
        //    Output, int *N, the number of points
        //
    {
        TableHeader h = new()
        {
            m = TableMisc.file_column_count(input_filename)
        };

        switch (h.m)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("DTABLE_HEADER_READ - Fatal error!");
                Console.WriteLine("  FILE_COLUMN_COUNT failed.");
                h.n = -1;
                return h;
        }

        h.n = TableMisc.file_row_count(input_filename);

        switch (h.n)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("DTABLE_HEADER_READ - Fatal error!");
                Console.WriteLine("  FILE_ROW_COUNT failed.");
                break;
        }

        return h;
    }

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
        List<string> output = new();

        try
        {
            //
            //  Write the data.
            //  For greater precision, try
            //
            //    output << "  " << setw(24) << setprecision(16) << table[i+j*m];
            //
            int j;
            for ( j = 0; j < n; j++ )
            {
                string cout = "";
                int i;
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
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("DTABLE_WRITE0 - Fatal error!");
            Console.WriteLine("  Could not open the output file.");
        }
    }
}