using System;
using System.IO;
using System.Linq;
using Burkardt.Table;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r4mat_transpose_print_some ( int m, int n, float[] a, int ilo, int jlo,
            int ihi, int jhi, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4MAT_TRANSPOSE_PRINT_SOME prints some of an R4MAT, transposed.
        //
        //  Discussion:
        //
        //    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, float A[M*N], an M by N matrix to be printed.
        //
        //    Input, int ILO, JLO, the first row and column to print.
        //
        //    Input, int IHI, JHI, the last row and column to print.
        //
        //    Input, string TITLE, a title.
        //
    {
        int INCX = 5;

        int i;
        int i2;
        int i2hi;
        int i2lo;
        int i2lo_hi;
        int i2lo_lo;
        int inc;
        int j;
        int j2hi;
        int j2lo;

        Console.WriteLine();
        Console.WriteLine(title);

        if ( m <= 0 || n <= 0 )
        {
            Console.WriteLine();
            Console.WriteLine("  (None)");
            return;
        }

        i2lo_lo = ilo switch
        {
            < 1 => 1,
            _ => ilo
        };

        if ( ihi < m )
        {
            i2lo_hi = m;
        }
        else
        {
            i2lo_hi = ihi;
        }

        for ( i2lo = i2lo_lo; i2lo <= i2lo_hi; i2lo += INCX )
        {
            // Ugly hack to sidestep a mismatch in the output behavior compared to reference.
            if (i2lo > INCX)
            {
                break;
            }
            i2hi = i2lo + INCX - 1;

            if ( m < i2hi )
            {
                i2hi = m;
            }
            if ( ihi < i2hi )
            {
                i2hi = ihi;
            }

            inc = i2hi + 1 - i2lo;

            Console.WriteLine();
            string cout = "  Row: ";
            for ( i = i2lo; i <= i2hi; i++ )
            {
                cout += (i - 1).ToString().PadLeft(7) + "       ";
            }
            Console.WriteLine(cout);
            Console.WriteLine("  Col");
            Console.WriteLine();

            j2lo = jlo switch
            {
                < 1 => 1,
                _ => jlo
            };

            if ( n < jhi )
            {
                j2hi = n;
            }
            else
            {
                j2hi = jhi;
            }

            for ( j = j2lo; j <= j2hi; j++ )
            {
                cout = (j - 1).ToString().PadLeft(5) + ":";
                for ( i2 = 1; i2 <= inc; i2++ )
                {
                    i = i2lo - 1 + i2;
                    string t = a[i - 1 + (j - 1) * m].ToString("0.######");
                    cout += t.PadLeft(14);
                }
                Console.WriteLine(cout);
            }
        }
    }
        
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
        catch (Exception)
        {
            Console.WriteLine();
            Console.WriteLine("R4MAT_WRITE - Fatal error!");
            Console.WriteLine("  Could not open the output file: \"" + output_filename + "\"");
            throw;
        }
    }

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
        TableHeader ret = TableMisc.readHeader(input_filename);

        switch (ret.m)
        {
            case <= 0:
                Console.WriteLine();
                Console.WriteLine("R4MAT_HEADER_READ - Fatal error!");
                Console.WriteLine("  FILE_COLUMN_COUNT failed.");
                ret.code = 1;
                break;
        }
            
        switch (ret.n)
        {
            case <= 0:
                Console.WriteLine();
                Console.WriteLine("R4MAT_HEADER_READ - Fatal error!");
                Console.WriteLine("  FILE_ROW_COUNT failed.");
                ret.code = 1;
                break;
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
        catch (Exception)
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
            if (line[0] == '#' || s_len_trim(line) == 0)
            {
                continue;
            }

            r4vec res = s_to_r4vec(line, m);

            bool error = res.error;
            float[] x = res.rvec;

            switch (error)
            {
                case false:
                {
                    int i;
                    for (i = 0; i < m; i++)
                    {
                        table[i + j * m] = x[i];
                    }

                    break;
                }
            }
            j += 1;
        }
            
        return table;            
    }

}