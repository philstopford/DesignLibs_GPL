﻿using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt.IO;

public static class PBMA
{
    public static void pbma_check_data(int xsize, int ysize, int[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PBMA_CHECK_DATA checks the data for an ASCII PBM file.
        //
        //  Discussion:
        //
        //    XSIZE and YSIZE must be positive, the pointers must not be null,
        //    and the data must be 0 or 1.
        //
        //  Example:
        //
        //    P1
        //    # feep.pbm
        //    24 7
        //    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        //    0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0
        //    0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 1 0
        //    0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 1 0
        //    0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0
        //    0 1 0 0 0 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 0 0 0 0
        //    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
        //
        //    Input, int *B, the array of XSIZE by YSIZE data values.
        //
    {
        int index = 0;
        int j;

        switch (xsize)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("PBMA_CHECK_DATA: Error!");
                Console.WriteLine("  XSIZE <= 0.");
                Console.WriteLine("  XSIZE = " + xsize + "");
                return;
        }

        switch (ysize)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("PBMA_CHECK_DATA: Error!");
                Console.WriteLine("  YSIZE <= 0.");
                Console.WriteLine("  YSIZE = " + ysize + "");
                return;
        }

        switch (b)
        {
            case null:
                Console.WriteLine("");
                Console.WriteLine("PBMA_CHECK_DATA: Error!");
                Console.WriteLine("  Null pointer to B.");
                return;
        }
            
        for (j = 0; j < ysize; j++)
        {
            int i;
            for (i = 0; i < xsize; i++)
            {
                switch (index)
                {
                    case < 0:
                        Console.WriteLine("");
                        Console.WriteLine("PBMA_CHECK_DATA - Fatal error!");
                        Console.WriteLine("  Negative data.");
                        Console.WriteLine("  B(" + i + "," + j + ")=" + b[index] + "");
                        return;
                    case > 1:
                        Console.WriteLine("");
                        Console.WriteLine("PBMA_CHECK_DATA - Fatal error!");
                        Console.WriteLine("  Data exceeds 1");
                        Console.WriteLine("  B(" + i + "," + j + ")=" + b[index] + "");
                        return;
                    default:
                        index++;
                        break;
                }
            }
        }
    }

    public static void pbma_example ( int xsize, int ysize, ref int[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PBMA_EXAMPLE sets up some ASCII PBM data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
        //
        //    Output, int *B, the array of XSIZE by YSIZE gray values.
        //
    {
        int i;
        int indexb = 0;
        double r;

        if ( xsize < ysize )
        {
            r = xsize / 3.0;
        }
        else
        {
            r = ysize / 3.0;
        }

        double xc = xsize / 2.0;
        double yc = ysize / 2.0;

        for ( i = 0; i < ysize; i++ )
        {
            double y = i;
            int j;
            for ( j = 0; j < xsize; j++ )
            {
                double x = j;
                double test = r - Math.Sqrt ( ( x - xc ) * ( x - xc )
                                              + 0.75 * ( y - yc ) * ( y - yc ) );
                b[indexb] = Math.Abs(test) switch
                {
                    <= 3.0 => 1,
                    _ => 0
                };
                indexb++;
            }
        }
    }

    public static void pbma_read(string file_in_name, ref int xsize, ref int ysize, ref int[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PBMA_READ reads the header and data from an ASCII PBM file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string FILE_IN_NAME, the name of the file.
        //
        //    Output, int &XSIZE, &YSIZE, the number of rows and columns of data.
        //
        //    Output, int **B, the array of XSIZE by YSIZE data values.
        //
    {
        string[] file_in;
        int file_index = 0;

        try
        {
            file_in = File.ReadAllLines(file_in_name);
        }
        catch (Exception)
        {
            Console.WriteLine();
            Console.WriteLine("PBMA_READ - Fatal error!");
            Console.WriteLine("  Cannot open the input file \"" + file_in_name + "\".");
            return;
        }

        //
        //  Read the header.
        //
        pbma_read_header(ref file_in, ref file_index, ref xsize, ref ysize);
        //
        //  Allocate storage for the data.
        //
        int numbytes = xsize * ysize * sizeof(int);

        b = new int[numbytes];
        //
        //  Read the data.
        //
        pbma_read_data(ref file_in, ref file_index, xsize, ysize, ref b);
    }

    public static void pbma_read_data(ref string[] file_in, ref int file_index, int xsize, int ysize, ref int[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PBMA_READ_DATA reads the data in an ASCII PBM file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, ifstream &FILE_IN, a pointer to the file.
        //
        //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
        //
        //    Output, int *B, the array of XSIZE by YSIZE data values.
        //
    {
        int index = 0;

        while (file_index < file_in.Length)
        {
            string[] tokens = Helpers.splitStringByWhitespace(file_in[file_index]);
            foreach (string t in tokens)
            {
                b[index] = Convert.ToInt32(t);
                index++;
            }

            file_index++;
        }
    }

    public static void pbma_read_header(ref string[] file_in, ref int file_index, ref int xsize, ref int ysize)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PBMA_READ_HEADER reads the header of an ASCII PBM file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, ifstream &FILE_IN, a pointer to the file.
        //
        //    Output, int &XSIZE, &YSIZE, the number of rows and columns of data.
        //
    {
        string rest = "";
        int step = 0;
        string word = "";
            
        while (true)
        {
            string 
                line;
            try
            {
                line = file_in[file_index];
                file_index++;
            }
            catch (Exception)
            {
                Console.WriteLine();
                Console.WriteLine("PBMA_READ_HEADER - Fatal error!");
                Console.WriteLine("  End of file.");
                return;
            }
                
            switch (line[0])
            {
                case '#':
                    continue;
            }

            switch (step)
            {
                case 0:
                {
                    typeMethods.s_word_extract_first(line, ref word, ref rest);

                    if (typeMethods.s_len_trim(word) == 0)
                    {
                        continue;
                    }

                    line = rest;

                    if (word[0] != 'P' && word[0] != 'p' ||
                        word[1] != '1')
                    {
                        Console.WriteLine();
                        Console.WriteLine("PBMA_READ_HEADER - Fatal error.");
                        Console.WriteLine("  Bad magic number = \"" + word + "\".");
                        return;
                    }

                    step = 1;
                    break;
                }
            }

            switch (step)
            {
                case 1:
                {
                    typeMethods.s_word_extract_first(line, ref word, ref rest);

                    if (typeMethods.s_len_trim(word) == 0)
                    {
                        continue;
                    }

                    xsize = Convert.ToInt32(word);
                    line = rest;
                    step = 2;
                    break;
                }
            }

            typeMethods.s_word_extract_first(line, ref word, ref rest);

            if (typeMethods.s_len_trim(word) == 0)
            {
                continue;
            }

            ysize = Convert.ToInt32(word);
            line = rest;
            break;
        }
    }

    public static void pbma_read_test(string file_in_name)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PBMA_READ_TEST tests the ASCII PBM read routines.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string FILE_IN_NAME, the name of the file.
        //
    {
        int xsize = 0;
        int ysize = 0;

        int[] b = null;
        //
        //  Read the data.
        //
        pbma_read(file_in_name, ref xsize, ref ysize, ref b);

        Console.WriteLine("");
        Console.WriteLine("  PBMA_READ was able to read \"" + file_in_name + "\".");
        //
        //  Check the data.
        //
        pbma_check_data(xsize, ysize, b);

        Console.WriteLine();
        Console.WriteLine("  PBMA_CHECK_DATA approved the data from the file.");

    }

    public static void pbma_write(string file_out_name, int xsize, int ysize, int[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PBMA_WRITE writes the header and data for an ASCII PBM file.
        //
        //  Example:
        //
        //    P1
        //    # feep.pbm
        //    24 7
        //    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        //    0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0
        //    0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 1 0
        //    0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 1 0
        //    0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0
        //    0 1 0 0 0 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 0 0 0 0
        //    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string FILE_OUT_NAME, the name of the file.
        //
        //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
        //
        //    Input, int *B, the array of XSIZE by YSIZE data values.
        //
    {
        List<string> file_out = new();
        //
        //  Write the header.
        //
        pbma_write_header(ref file_out, file_out_name, xsize, ysize);
        //
        //  Write the data.
        //
        pbma_write_data(ref file_out, xsize, ysize, b);

        try
        {
            File.WriteAllLines(file_out_name, file_out);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("PBMA_WRITE - Fatal error!");
            Console.WriteLine("  Cannot open the output file \"" + file_out_name + "\".");
        }
    }

    public static void pbma_write_data(ref List<string> file_out, int xsize, int ysize, int[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PBMA_WRITE_DATA writes the data for an ASCII PBM file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, ofstream &FILE_OUT, a pointer to the file.
        //
        //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
        //
        //    Input, int *B, the arrays of XSIZE by YSIZE data values.
        //
    {
        int j;

        int numval = 0;
        string cout = "";
        int bIndex = 0;

        for (j = 0; j < ysize; j++)
        {
            int i;
            for (i = 0; i < xsize; i++)
            {
                cout += b[bIndex] + " ";
                numval += 1;
                bIndex++;

                if (numval % 12 == 0 || i == xsize - 1 || numval == xsize * ysize)
                {
                    file_out.Add(cout);
                    cout = "";
                }
                else
                {
                    cout += " ";
                }

            }
        }
    }

    public static void pbma_write_header(ref List<string> file_out, string file_out_name, int xsize,
            int ysize)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PBMA_WRITE_HEADER writes the header of an ASCII PBM file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, ofstream &FILE_OUT, a pointer to the file.
        //
        //    Input, string FILE_OUT_NAME, the name of the file.
        //
        //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
        //
    {
        file_out.Add("P1");
        file_out.Add("# " + file_out_name + " created by PBMA_IO::PBMA_WRITE.");
        file_out.Add(xsize + "  " + ysize + "");
    }
        
    public static void pbma_write_test ( string file_out_name )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PBMA_WRITE_TEST tests the ASCII PBM routines.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string FILE_OUT_NAME, the name of the file.
        //
    {
        const int xsize = 200;
        const int ysize = 200;
        //
        //  Allocate memory.
        //
        int[] b = new int[xsize * ysize];
        //
        //  Set the data.
        //
        pbma_example ( xsize, ysize, ref b );
        //
        //  Write the data to the file.
        //
        pbma_write ( file_out_name, xsize, ysize, b );
    }
}