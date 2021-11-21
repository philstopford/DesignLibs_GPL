﻿using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt.IO;

public static class PPM_ASCII
{
    public static bool ppma_check_data(int xsize, int ysize, ref int rgb_max, int[] r,
            int[] g, int[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PPMA_CHECK_DATA checks the data for an ASCII portable pixel map file.
        //
        //  Discussion:
        //
        //    XSIZE and YSIZE must be positive, the pointers must not be null,
        //    and the data must be nonnegative and no greater than RGB_MAX.
        //
        //  Example:
        //
        //    P3
        //    # feep.ppm
        //    4 4
        //    15
        //     0  0  0    0  0  0    0  0  0   15  0 15
        //     0  0  0    0 15  7    0  0  0    0  0  0
        //     0  0  0    0  0  0    0 15  7    0  0  0
        //    15  0 15    0  0  0    0  0  0    0  0  0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 February 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
        //
        //    Input, int RGB_MAX, the maximum RGB value.
        //
        //    Input, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.
        //
        //    Output, bool PPMA_CHECK_DATA, is
        //    true, if an error was detected, or
        //    false, if the data was legal.
        //
    {
        char c = 'x';
        int index = 0;
        int k;

        switch (xsize)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("PPMA_CHECK_DATA: Fatal error!");
                Console.WriteLine("  xsize <= 0.");
                Console.WriteLine("  xsize = " + xsize + "");
                return true;
        }

        switch (ysize)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("PPMA_CHECK_DATA: Fatal error!");
                Console.WriteLine("  ysize <= 0.");
                Console.WriteLine("  ysize = " + ysize + "");
                return true;
        }

        switch (r)
        {
            case null:
                Console.WriteLine("");
                Console.WriteLine("PPMA_CHECK_DATA: Fatal error!");
                Console.WriteLine("  Null pointer to R.");
                return true;
        }

        switch (g)
        {
            case null:
                Console.WriteLine("");
                Console.WriteLine("PPMA_CHECK_DATA: Fatal error!");
                Console.WriteLine("  Null pointer to G.");
                return true;
        }

        switch (b)
        {
            case null:
                Console.WriteLine("");
                Console.WriteLine("PPMA_CHECK_DATA: Fatal error!");
                Console.WriteLine("  Null pointer to B.");
                return true;
        }

        for ( k = 0; k < 3; k++ )
        {
            float val = 0;

            switch (k)
            {
                case 0:
                    val = r[index];
                    c = 'R';
                    break;
                case 1:
                    val = g[index];
                    c = 'G';
                    break;
                case 2:
                    val = b[index];
                    c = 'B';
                    break;
            }

            int j;
            for ( j = 0; j < ysize; j++ )
            {
                int i;
                for ( i = 0; i < xsize; i++ )
                {
                    switch (val)
                    {
                        case < 0:
                            Console.WriteLine("");
                            Console.WriteLine("PPMA_CHECK_DATA - Fatal error!");
                            Console.WriteLine("  Negative data.");
                            Console.WriteLine("  " + c + "(" + i + "," + j + ")=" + val + "");
                            return true;
                    }

                    if ( rgb_max < val )
                    {
                        Console.WriteLine("");
                        Console.WriteLine("PPMA_CHECK_DATA - Fatal error!");
                        Console.WriteLine("  Data exceeds RGB_MAX = " + rgb_max + "");
                        Console.WriteLine("  " + c + "(" + i + "," + j + ")=" + val + "");
                        return true;
                    }

                    index += 1;
                }
            } 
            if (index >= r.Length)
            {
                break;
            }
        }

        return false;
    }

    public static bool ppma_example(int xsize, int ysize, ref int[] r, ref int[] g, ref int[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PPMA_EXAMPLE sets up some RGB data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 February 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
        //
        //    Output, int *R, *G, *B, the arrays of XSIZE by YSIZE RGB values.
        //
        //    Output, bool PPMA_EXAMPLE, is
        //    false, if no error occurred,
        //    true, if an error occurred.
        //
    {
        int i;

        int index = 0;
        for (i = 0; i < ysize; i++)
        {
            float y = (ysize + 1 - i) / (float) (ysize - 1);
            int j;
            for (j = 0; j < xsize; j++)
            {
                float x = j / (float) (xsize - 1);

                float f1 = 4.0f * (x - 0.5f) * (x - 0.5f);
                float f2 = (float) Math.Sin(3.14159265f * x);

                if (y <= f1)
                {
                    r[index] = (int) (255.0 * f1);
                }
                else
                {
                    r[index] = 50;
                }

                if (y <= f2)
                {
                    g[index] = (int) (255.0 * f2);
                }
                else
                {
                    g[index] = 150;
                }

                if (y <= x)
                {
                    b[index] = (int) (255.0 * x);
                }
                else
                {
                    b[index] = 250;
                }

                index++;
            }
        }

        return false;
    }

    public static bool ppma_read(string input_name, ref int xsize, ref int ysize, ref int rgb_max,
            ref int[] r, ref int[] g, ref int[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PPMA_READ reads the header and data from an ASCII portable pixel map file.
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
        //    Input, string INPUT_NAME, the name of the file containing the ASCII
        //    portable pixel map data.
        //
        //    Output, int &XSIZE, &YSIZE, the number of rows and columns of data.
        //
        //    Output, int &RGB_MAX, the maximum RGB value.
        //
        //    Output, int **R, **G, **B, the arrays of XSIZE by YSIZE data values.
        //
        //    Output, bool PPMA_READ, is
        //    true, if an error was detected, or
        //    false, if the file was read.
        //
    {
        string[] input;

        int index = 0;

        try
        {
            input = File.ReadAllLines(input_name);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("PPMA_READ - Fatal error!");
            Console.WriteLine("  Cannot open the input file \"" + input_name + "\".");
            return true;
        }

        //
        //  Read the header.
        //
        bool error = ppma_read_header(ref input, ref index, ref xsize, ref ysize, ref rgb_max);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("PPMA_READ - Fatal error!");
                Console.WriteLine("  PPMA_READ_HEADER failed.");
                return true;
        }

        //
        //  Allocate storage for the data.
        //
        int numbytes = xsize * ysize;

        r = new int[numbytes];
        g = new int[numbytes];
        b = new int[numbytes];
        //
        //  Read the data.
        //
        error = ppma_read_data(ref input, ref index, xsize, ysize, ref r, ref g, ref b);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("PPMA_READ - Fatal error!");
                Console.WriteLine("  PPMA_READ_DATA failed.");
                return true;
            default:
                return false;
        }
    }

    public static bool ppma_read_data(ref string[] input, ref int index, int xsize, int ysize, ref int[] r,
            ref int[] g, ref int[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PPMA_READ_DATA reads the data in an ASCII portable pixel map file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 August 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, ifstream &INPUT, a pointer to the file containing the ASCII
        //    portable pixel map data.
        //
        //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
        //
        //    Output, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.
        //
        //    Output, bool PPMA_READ_DATA, is
        //    true, if an error was detected, or
        //    false, if the data was read.
        //
    {
        int arrayIndex = 0;
        while (index < input.Length)
        {

            string[] tokens = Helpers.splitStringByWhitespace(input[index]);
            for (int i = 0; i < tokens.Length; i+=3)
            {
                r[arrayIndex] = Convert.ToInt32(tokens[i]);
                g[arrayIndex] = Convert.ToInt32(tokens[i+1]);
                b[arrayIndex] = Convert.ToInt32(tokens[i+2]);

                arrayIndex++;
            }

            index++;
        }

        return false;
    }

    public static bool ppma_read_header(ref string[] input, ref int index, ref int xsize, ref int ysize,
            ref int rgb_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PPMA_READ_HEADER reads the header of an ASCII portable pixel map file.
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
        //    Input, ifstream &INPUT, a pointer to the file containing the ASCII
        //    portable pixel map data.
        //
        //    Output, int &XSIZE, &YSIZE, the number of rows and columns of data.
        //
        //    Output, int &RGB_MAX, the maximum RGB value.
        //
        //    Output, bool PPMA_READ_HEADER, is
        //    true, if an error was detected, or
        //    false, if the header was read.
        //
    {
        string rest = "";
        string word = "";

        int step = 0;

        while (true)
        {
            string line;
            try
            {
                line = input[index];
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("PPMA_READ_HEADER - Fatal error!");
                Console.WriteLine("  End of file.");
                return true;
            }

            switch (line[0])
            {
                case '#':
                    index++;
                    continue;
            }

            switch (step)
            {
                case 0:
                {
                    typeMethods.s_word_extract_first(line, ref word, ref rest);

                    if (typeMethods.s_len_trim(word) <= 0)
                    {
                        continue;
                    }

                    if (!typeMethods.s_eqi(word, "P3"))
                    {
                        Console.WriteLine("");
                        Console.WriteLine("PPMA_READ_HEADER - Fatal error.");
                        Console.WriteLine("  Bad magic number = \"" + word + "\".");
                        return true;
                    }

                    line = rest;
                    step = 1;
                    index++;
                    continue;
                }
                case 1:
                {
                    typeMethods.s_word_extract_first(line, ref word, ref rest);

                    if (typeMethods.s_len_trim(word) <= 0)
                    {
                        continue;
                    }

                    xsize = Convert.ToInt32(word);
                    line = rest;
                    step = 2;
                    break;
                }
            }

            switch (step)
            {
                case 2:
                {
                    typeMethods.s_word_extract_first(line, ref word, ref rest);

                    if (typeMethods.s_len_trim(word) <= 0)
                    {
                        continue;
                    }

                    ysize = Convert.ToInt32(word);
                    line = rest;
                    step = 3;
                    index++;
                    continue;
                }
            }

            if (step == 3)
            {
                typeMethods.s_word_extract_first(line, ref word, ref rest);

                if (typeMethods.s_len_trim(word) <= 0)
                {
                    continue;
                }

                rgb_max = Convert.ToInt32(word);
                line = rest;
                index++;
                break;
            }

            index++;
        }

        return false;
    }

    public static bool ppma_read_test(string input_name)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PPMA_READ_TEST tests the ASCII portable pixel map read routines.
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
        //    Input, string INPUT_NAME, the name of the file containing the ASCII
        //    portable pixel map data.
        //
        //    Output, bool PPMA_READ_TEST, is
        //    true, if an error was detected, or
        //    false, if the test was carried out.
        //
    {
        int rgb_max = 0;
        int xsize = 0;
        int ysize = 0;

        int[] r = null;
        int[] g = null;
        int[] b = null;
        //
        //  Read the data.
        //
        bool error = ppma_read(input_name, ref xsize, ref ysize, ref rgb_max, ref r, ref g, ref b);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("PPMA_READ_TEST - Fatal error!");
                Console.WriteLine("  PPMA_READ failed.");

                return true;
        }

        Console.WriteLine("");
        Console.WriteLine("PPMA_READ_TEST:");
        Console.WriteLine("  PPMA_READ was able to read \"" + input_name + "\".");
        //
        //  Check the data.
        //
        error = ppma_check_data(xsize, ysize, ref rgb_max, r, g, b);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("PPMA_READ_TEST - Fatal error!");
                Console.WriteLine("  PPMA_CHECK_DATA reports bad data in the file.");
                return true;
        }

        Console.WriteLine("");
        Console.WriteLine("PPMA_READ_TEST:");
        Console.WriteLine("  PPMA_CHECK_DATA has approved the data from the file.");

        return false;
    }

    public static bool ppma_write(string file_out_name, int xsize, int ysize, int[] r,
            int[] g, int[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PPMA_WRITE writes the header and data for an ASCII portable pixel map file.
        // 
        //  Example:
        //
        //    P3
        //    # feep.ppm
        //    4 4
        //    15
        //     0  0  0    0  0  0    0  0  0   15  0 15
        //     0  0  0    0 15  7    0  0  0    0  0  0
        //     0  0  0    0  0  0    0 15  7    0  0  0
        //    15  0 15    0  0  0    0  0  0    0  0  0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        // 
        //    28 February 2003
        // 
        //  Author:
        // 
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string FILE_OUT_NAME, the name of the file to contain the ASCII
        //    portable pixel map data.
        //
        //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
        //
        //    Input, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.
        //
        //    Output, bool PPMA_WRITE, is
        //    true, if an error was detected, or
        //    false, if the file was written.
        //
    {
        List<string> file_out = new();
        int j;

        //
        //  Compute the maximum.
        //
        int rgb_max = 0;
        int index = 0;

        for (j = 0; j < ysize; j++)
        {
            int i;
            for (i = 0; i < xsize; i++)
            {
                if (rgb_max < r[index])
                {
                    rgb_max = r[index];
                }

                if (rgb_max < g[index])
                {
                    rgb_max = g[index];
                }

                if (rgb_max < b[index])
                {
                    rgb_max = b[index];
                }

                index++;
            }
        }

        //
        //  Write the header.
        //
        bool error = ppma_write_header(ref file_out, file_out_name, xsize, ysize, rgb_max);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("PPMA_WRITE - Fatal error!");
                Console.WriteLine("  PPMA_WRITE_HEADER failed.");
                return true;
        }

        //
        //  Write the data.
        //
        error = ppma_write_data(ref file_out, xsize, ysize, r, g, b);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("PPMA_WRITE - Fatal error!");
                Console.WriteLine("  PPMA_WRITE_DATA failed.");
                return true;
            default:
                try
                {
                    File.WriteAllLines(file_out_name, file_out);
                }
                catch
                {
                    Console.WriteLine("");
                    Console.WriteLine("PPMA_WRITE - Fatal error!");
                    Console.WriteLine("  Cannot open the output file \"" + file_out_name + "\".");
                    return true;
                }


                return false;
        }
    }

    public static bool ppma_write_data(ref List<string> file_out, int xsize, int ysize, int[] r,
            int[] g, int[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PPMA_WRITE_DATA writes the data for an ASCII portable pixel map file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 February 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, ofstream &FILE_OUT, a pointer to the file to contain the ASCII
        //    portable pixel map data.
        //
        //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
        //
        //    Input, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.
        //
        //    Output, bool PPMA_WRITE_DATA, is
        //    true, if an error was detected, or
        //    false, if the data was written.
        //
    {
        int j;

        int rgb_num = 0;

        int index = 0;
        string tmp = "";
        for (j = 0; j < ysize; j++)
        {
            int i;
            for (i = 0; i < xsize; i++)
            {
                tmp += r[index] + " " + g[index] + " " + b[index];
                rgb_num += 3;

                if (rgb_num % 12 == 0 || i == xsize - 1 || rgb_num == 3 * xsize * ysize)
                {
                    file_out.Add(tmp);
                    tmp = "";
                }
                else
                {
                    tmp += " ";
                }
                index++;
            }
        }

        return false;
    }

    public static bool ppma_write_header(ref List<string> file_out, string file_out_name, int xsize,
            int ysize, int rgb_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PPMA_WRITE_HEADER writes the header of an ASCII portable pixel map file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 February 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, ofstream &FILE_OUT, a pointer to the file to contain the ASCII
        //    portable pixel map data.
        //
        //    Input, string FILE_OUT_NAME, the name of the file.
        //
        //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
        //
        //    Input, int RGB_MAX, the maximum RGB value.
        //
        //    Output, bool PPMA_WRITE_HEADER, is
        //    true, if an error was detected, or
        //    false, if the header was written.
        //
    {
        file_out.Add("P3");
        file_out.Add("# " + file_out_name + " created by PPMA_WRITE.C.");
        file_out.Add(xsize + "  " + ysize + "");
        file_out.Add(rgb_max + "");

        return false;
    }
            
    public static bool ppma_write_test ( string output_name )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PPMA_WRITE_TEST tests the ASCII portable pixel map write routines.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 February 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string OUTPUT_NAME, the name of the file to contain the ASCII
        //    portable pixel map data.
        //
        //    Output, bool PPMA_WRITE_TEST, equals
        //    true, if the test could not be carried out,
        //    false, if the test was carried out.
        //
    {
        const int xsize = 300;
        const int ysize = 300;
        //
        //  Allocate memory.
        //
        int[] r = new int[xsize * ysize];
        int[] g = new int[xsize * ysize];
        int[] b = new int[xsize * ysize];
        //
        //  Set the data.
        //
        bool error = ppma_example ( xsize, ysize, ref r, ref g, ref b );

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("PPMA_WRITE_TEST - Fatal error!");
                Console.WriteLine("  PPMA_EXAMPLE failed.");
                return true;
        }
        //
        //  Write the data to the file.
        //
        error = ppma_write ( output_name, xsize, ysize, r, g, b );
 
        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("PPMA_WRITE_TEST - Fatal error!");
                Console.WriteLine("  PPMA_WRITE failed.");
                return true;
            default:
                return false;
        }
    }
}