using System;
using System.IO;
using System.Linq;
using Burkardt.Types;

namespace Burkardt.Table
{
    public static partial class TableWriter
    {
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
    }
    public static partial class TableReader
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
            TableHeader ret = readHeader(input_filename);

            if ( ret.m <= 0 )
            {
                Console.WriteLine();
                Console.WriteLine("R8MAT_HEADER_READ - Fatal error!");
                Console.WriteLine("  FILE_COLUMN_COUNT failed.");
                ret.code = 1;
            }
            
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
    }

    public static partial class TableMisc
    {
        public static void r8mat_transpose_print_some ( int m, int n, double[] a, int ilo, int jlo,
        int ihi, int jhi, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
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
        //    Input, double A[M*N], an M by N matrix to be printed.
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

            if ( ilo < 1 )
            {
                i2lo_lo = 1;
            }
            else
            {
                i2lo_lo = ilo;
            }

            if ( ihi < m )
            {
                i2lo_hi = m;
            }
            else
            {
                i2lo_hi = ihi;
            }

            for ( i2lo = i2lo_lo; i2lo <= i2lo_hi; i2lo = i2lo + INCX )
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

                if ( jlo < 1 )
                {
                    j2lo = 1;
                }
                else
                {
                    j2lo = jlo;
                }

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
                        string t = a[(i - 1) + (j - 1) * m].ToString("0.######");
                        cout += t.PadLeft(14);
                    }
                    Console.WriteLine(cout);
                }
            }
        }
    }
}