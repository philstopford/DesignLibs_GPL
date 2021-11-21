using System;
using System.IO;
using System.Text;
using Burkardt.Types;

namespace Burkardt.IO;

public static class Files
{
    public static string file_name_ext_swap ( string filename, string ext )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FILE_NAME_EXT_SWAP replaces the current "extension" of a file name.
        //
        //  Discussion:
        //
        //    The "extension" of a file name is the string of characters
        //    that appears after the LAST period in the name.  A file
        //    with no period, or with a period as the last character
        //    in the name, has a "null" extension.
        //
        //  Example:
        //
        //          Input           Output
        //    ================     ==================
        //    FILENAME     EXT     FILE_NAME_EXT_SWAP
        //
        //    bob.for      obj     bob.obj
        //    bob.bob.bob  txt     bob.bob.txt
        //    bob          yak     bob.yak
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string FILENAME, a file name.
        //
        //    Input, string EXT, the extension to be added to the file name.
        //
        //    Output, string FILE_NAME_EXT_SWAP, the file name with the new extension.
        //
    {
        //
        //  Look for the LAST occurrence of a period.
        //
        int i = filename.LastIndexOf( ".", StringComparison.Ordinal);

        string filename2 = i switch
        {
            -1 => filename + "." + ext,
            _ => string.Concat(filename.AsSpan(0, i + 1), ext)
        };

        return filename2;
    }
        
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
        try
        {
            string[] file_in = File.ReadAllLines(file_in_name);

            int j = 0;

            foreach (string line in file_in)
            {

                if (line[0] == '#' || typeMethods.s_len_trim(line) == 0)
                {
                    continue;
                }

                r8vec t = typeMethods.s_to_r8vec(line, dim_num);
                double[] x = t.rvec;
                bool error = t.error;

                switch (error)
                {
                    case true:
                        continue;
                }

                int i;
                for (i = 0; i < dim_num; i++)
                {
                    r[i + j * dim_num] = x[i];
                }

                j += 1;

            }

            Console.WriteLine("");
            Console.WriteLine("DATA_READ:");
            Console.WriteLine("  Read coordinate data from file.");

        }
        catch (Exception)
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
        int i;

        int lens = filename.Length;

        switch (lens)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("FILENAME_INC - Fatal error!");
                Console.WriteLine("  The input string is empty.");
                return;
        }

        int change = 0;

        for (i = lens - 1; 0 <= i; i--)
        {
            char c = filename[i];

            switch (c)
            {
                case >= '0' and <= '9':
                {
                    change += 1;

                    switch (c)
                    {
                        case '9':
                        {
                            c = '0';
                            StringBuilder sb = new(filename)
                            {
                                [i] = c
                            };
                            filename = sb.ToString();
                            break;
                        }
                        default:
                        {
                            c = Convert.ToChar(Convert.ToInt32(c) + 1);
                            StringBuilder sb = new(filename)
                            {
                                [i] = c
                            };
                            filename = sb.ToString();
                            return;
                        }
                    }

                    break;
                }
            }
        }

        switch (change)
        {
            //
            //  No digits were found.  Return blank.
            //
            case 0:
            {
                StringBuilder sb = new(filename);
                for (i = lens - 1; 0 <= i; i--)
                {
                    sb[i] = ' ';
                }
                filename = sb.ToString();
                break;
            }
        }
    }

    public static void file_name_inc_nowrap(ref string filename_)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FILE_NAME_INC_NOWRAP increments a partially numeric file name.
        //
        //  Discussion:
        //
        //    It is assumed that the digits in the name, whether scattered or
        //    connected, represent a number that is to be increased by 1 on
        //    each call.  If this number is all 9's on input, the output number
        //    is all 0's.  Non-numeric letters of the name are unaffected.
        //
        //    If the (nonempty) name contains no digits, or all the digits are
        //    9, then the empty string is returned.
        //
        //    If the empty string is input, the routine stops.
        //
        //  Example:
        //
        //      Input            Output
        //      -----            ------
        //      "a7to11.txt"     "a7to12.txt"  (typical case.  Last digit incremented)
        //      "a7to99.txt"     "a8to00.txt"  (last digit incremented, with carry.)
        //      "a8to99.txt"     "a9to00.txt"
        //      "a9to99.txt"     " "
        //      "cat.txt"        " "
        //      " "              STOP!
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
        int i;

        char[] filename = filename_.ToCharArray();

        int lens = filename.Length;


        switch (lens)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("FILE_NAME_INC_NOWRAP - Fatal error!");
                Console.WriteLine("  The input string is empty.");
                return;
        }

        int change = 0;
        int carry = 0;

        for (i = lens - 1; 0 <= i; i--)
        {
            char c = filename[i];

            switch (c)
            {
                case >= '0' and <= '9':
                {
                    change += 1;
                    carry = 0;

                    switch (c)
                    {
                        case '9':
                            carry = 1;
                            c = '0';
                            filename[i] = c;
                            break;
                        default:
                            c = (char) (Convert.ToInt32(c) + 1);
                            filename[i] = c;
                            return;
                    }

                    break;
                }
            }
        }

        switch (carry)
        {
            //
            //  Unsatisfied carry.  The input digits were all 9.  Return blank.
            //
            case 1:
            {
                for (i = lens - 1; 0 <= i; i--)
                {
                    filename[i] = ' ';
                }

                break;
            }
        }

        switch (change)
        {
            //
            //  No digits were found.  Return blank.
            //
            case 0:
            {
                for (i = lens - 1; 0 <= i; i--)
                {
                    filename[i] = ' ';
                }

                break;
            }
        }

    }
}