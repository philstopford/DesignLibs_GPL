using System;
using System.IO;
using System.Text;
using Burkardt.Types;

namespace Burkardt.IO
{
    public static class Files
    {
        public static void data_read(string file_in_name, int dim_num, int n, ref double[] r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DATA_READ reads generator coordinate data from a file.
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
            //    03 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, char *FILE_IN_NAME, the name of the input file.
            //
            //    Input, int DIM_NUM, the number of spatial dimensions.
            //
            //    Input, int N, the number of points.  The program
            //    will stop reading data once N values have been read.
            //
            //    Output, double R[DIM_NUM*N], the point coordinates.
            //
        {
            bool error;
            string[] file_in;
            int i;
            int j;
            double[] x;

            try
            {
                file_in = File.ReadAllLines(file_in_name);

                j = 0;

                foreach (string line in file_in)
                {

                    if (line[0] == '#' || typeMethods.s_len_trim(line) == 0)
                    {
                        continue;
                    }

                    r8vec t = typeMethods.s_to_r8vec(line, dim_num);
                    x = t.rvec;
                    error = t.error;

                    if (error)
                    {
                        continue;
                    }

                    for (i = 0; i < dim_num; i++)
                    {
                        r[i + j * dim_num] = x[i];
                    }

                    j = j + 1;

                }

                Console.WriteLine("");
                Console.WriteLine("DATA_READ:");
                Console.WriteLine("  Read coordinate data from file.");

            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("DATA_READ - Fatal error!");
                Console.WriteLine("  Could not open the input file: \"" + file_in_name + "\"");
            }

        }

        public static void filename_inc(ref string filename)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FILENAME_INC increments a partially numeric file name.
            //
            //  Discussion:
            //
            //    It is assumed that the digits in the name, whether scattered or
            //    connected, represent a number that is to be increased by 1 on
            //    each call.  If this number is all 9's on input, the output number
            //    is all 0's.  Non-numeric letters of the name are unaffected.
            //
            //    If the name is empty, then the routine stops.
            //
            //    If the name contains no digits, the empty string is returned.
            //
            //  Example:
            //
            //      Input            Output
            //      -----            ------
            //      "a7to11.txt"     "a7to12.txt"  (typical case.  Last digit incremented)
            //      "a7to99.txt"     "a8to00.txt"  (last digit incremented, with carry.)
            //      "a9to99.txt"     "a0to00.txt"  (wrap around)
            //      "cat.txt"        " "           (no digits to increment)
            //      " "              STOP!         (error)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 November 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input/output, string *FILENAME, the filename to be incremented.
            //
        {
            char c;
            int change;
            int i;
            int lens;

            lens = filename.Length;

            if (lens <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("FILENAME_INC - Fatal error!");
                Console.WriteLine("  The input string is empty.");
                return;
            }

            change = 0;

            for (i = lens - 1; 0 <= i; i--)
            {
                c = (filename)[i];

                if ('0' <= c && c <= '9')
                {
                    change = change + 1;

                    if (c == '9')
                    {
                        c = '0';
                        StringBuilder sb = new StringBuilder(filename);
                        sb[i] = c;
                        filename = sb.ToString();
                    }
                    else
                    {
                        c = Convert.ToChar(Convert.ToInt32(c) + 1);
                        StringBuilder sb = new StringBuilder(filename);
                        sb[i] = c;
                        filename = sb.ToString();
                        return;
                    }
                }
            }

            //
            //  No digits were found.  Return blank.
            //
            if (change == 0)
            {
                StringBuilder sb = new StringBuilder(filename);
                for (i = lens - 1; 0 <= i; i--)
                {
                    sb[i] = ' ';
                }
                filename = sb.ToString();
            }
        }
    }
}