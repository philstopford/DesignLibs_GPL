﻿using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt
{
    public static class PGMA
    {
        public static void pgma_check_data(int xsize, int ysize, int maxg, ref int[] g)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMA_CHECK_DATA checks the data for an ASCII PGM file.
            //
            //  Discussion:
            //
            //    XSIZE and YSIZE must be positive, the pointers must not be null,
            //    and the data must be nonnegative and no greater than MAXG.
            //
            //  Example:
            //
            //    P2
            //    # feep.pgm
            //    24 7
            //    15
            //    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
            //    0  3  3  3  3  0  0  7  7  7  7  0  0 11 11 11 11  0  0 15 15 15 15  0
            //    0  3  0  0  0  0  0  7  0  0  0  0  0 11  0  0  0  0  0 15  0  0 15  0
            //    0  3  3  3  0  0  0  7  7  7  0  0  0 11 11 11  0  0  0 15 15 15 15  0
            //    0  3  0  0  0  0  0  7  0  0  0  0  0 11  0  0  0  0  0 15  0  0  0  0
            //    0  3  0  0  0  0  0  7  7  7  7  0  0 11 11 11 11  0  0 15  0  0  0  0
            //    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
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
            //    Input, int MAXG, the maximum gray value.
            //
            //    Input, int *G, the array of XSIZE by YSIZE data values.
            //
        {
            int i;
            int[] index;
            int j;

            if (xsize <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("PGMA_CHECK_DATA: Error!");
                Console.WriteLine("  XSIZE <= 0.");
                Console.WriteLine("  XSIZE = " + xsize + "");
                return;
            }

            if (ysize <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("PGMA_CHECK_DATA: Error!");
                Console.WriteLine("  YSIZE <= 0.");
                Console.WriteLine("  YSIZE = " + ysize + "");
                return;
            }

            if (g == null)
            {
                Console.WriteLine("");
                Console.WriteLine("PGMA_CHECK_DATA: Error!");
                Console.WriteLine("  Null pointer to g.");
                return;
            }

            index = g;

            int iCount = 0;

            for (j = 0; j < ysize; j++)
            {
                for (i = 0; i < xsize; i++)
                {
                    if (iCount < 0)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("PGMA_CHECK_DATA - Fatal error!");
                        Console.WriteLine("  Negative data.");
                        Console.WriteLine("  G(" + i + "," + j + ")=" + index[iCount] + "");
                        return;
                    }
                    else if (maxg < iCount)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("PGMA_CHECK_DATA - Fatal error!");
                        Console.WriteLine("  Data exceeds MAXG = " + maxg + "");
                        Console.WriteLine("  G(" + i + "," + j + ")=" + index[iCount] + "");
                        return;
                    }

                    iCount = iCount + 1;
                }
            }
        }

        public static void pgma_example(int xsize, int ysize, ref int[] g)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMA_EXAMPLE sets up some data for an ASCII PGM file.
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
            //    Output, int *G, the array of XSIZE by YSIZE gray values.
            //
        {
            int i;
            int[] indexg;
            int j;
            int periods = 3;
            double PI = 3.14159265;
            double x;
            double y;

            indexg = g;
            int iCount = 0;

            for (i = 0; i < ysize; i++)
            {
                y = ((double) (2 * i)) / ((double) (ysize - 1)) - 1.0;

                for (j = 0; j < xsize; j++)
                {
                    x = (2.0 * PI * (double) (periods * j)) / ((double) (xsize - 1));

                    indexg[iCount] = (int) (20.0 * (Math.Sin(x) - y + 2));

                    iCount = iCount + 1;
                }
            }
        }

        public static void pgma_read(string input_name, ref int xsize, ref int ysize, ref int maxg, ref int[] g)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMA_READ reads the header and data from an ASCII PGM file.
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
            //    Input, string INPUT_NAME, the name of the file.
            //
            //    Output, int &XSIZE, &YSIZE, the number of rows and columns of data.
            //
            //    Output, int &MAXG, the maximum gray value.
            //
            //    Output, int **G, the array of XSIZE by YSIZE data values.
            //
        {
            string[] input;
            int numbytes;
            int inputIndex = 0;

            try
            {
                input = File.ReadAllLines(input_name);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("PGMA_READ - Fatal error!");
                Console.WriteLine("  Cannot open the input file \"" + input_name + "\".");
                return;
            }

            //
            //  Read the header.
            //
            pgma_read_header(ref input, ref inputIndex, ref xsize, ref ysize, ref maxg);
            //
            //  Allocate storage for the data.
            //
            numbytes = xsize * ysize * sizeof(int);

            g = new int[numbytes];
            //
            //  Read the data.
            //
            pgma_read_data(ref input, ref inputIndex, xsize, ysize, ref g);
        }

        public static void pgma_read_data(ref string[] input, ref int inputIndex, int xsize, int ysize, ref int[] g)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMA_READ_DATA reads the data in an ASCII PGM file.
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
            //    Input, ifstream &INPUT, a pointer to the file containing the data.
            //
            //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
            //
            //    Output, int *G, the array of XSIZE by YSIZE data values.
            //
        {
            int i;
            int j;

            int gIndex = 0;
            for (j = 0; j < ysize; j++)
            {
                for (i = 0; i < xsize; i++)
                {
                    g[gIndex] = Convert.ToInt32(input[inputIndex]);
                    gIndex++;
                    inputIndex++;
                }
            }
        }

        public static void pgma_read_header(ref string[] input, ref int inputIndex, ref int xsize, ref int ysize,
                ref int maxg)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMA_READ_HEADER reads the header of an ASCII PGM file.
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
            //    Input, ifstream &INPUT, a pointer to the file.
            //
            //    Output, int &XSIZE, &YSIZE, the number of rows and columns of data.
            //
            //    Output, int &MAXG, the maximum gray value.
            //
        {
            string line;
            string rest = "";
            int step;
            string word = "";

            step = 0;

            while (true)
            {
                try
                {
                    line = input[inputIndex];
                    inputIndex++;
                }
                catch (Exception e)
                {
                    Console.WriteLine("");
                    Console.WriteLine("PGMA_READ_HEADER - Fatal error!");
                    Console.WriteLine("  End of file.");
                    return;
                }

                if (line[0] == '#')
                {
                    continue;
                }

                if (step == 0)
                {
                    typeMethods.s_word_extract_first(line, ref word, ref rest);

                    if (typeMethods.s_len_trim(word) == 0)
                    {
                        continue;
                    }

                    line = rest;

                    if ((word[0] != 'P' && word[0] != 'p') ||
                        word[1] != '2')
                    {
                        Console.WriteLine("");
                        Console.WriteLine("PGMA_READ_HEADER - Fatal error.");
                        Console.WriteLine("  Bad magic number = \"" + word + "\".");
                        return;
                    }

                    step = 1;
                }

                if (step == 1)
                {
                    typeMethods.s_word_extract_first(line, ref word, ref rest);

                    if (typeMethods.s_len_trim(word) == 0)
                    {
                        continue;
                    }

                    xsize = Convert.ToInt32(word);
                    line = rest;
                    step = 2;
                }

                if (step == 2)
                {
                    typeMethods.s_word_extract_first(line, ref word, ref rest);

                    if (typeMethods.s_len_trim(word) == 0)
                    {
                        continue;
                    }

                    ysize = Convert.ToInt32(word);
                    line = rest;
                    step = 3;
                }

                if (step == 3)
                {
                    typeMethods.s_word_extract_first(line, ref word, ref rest);

                    if (typeMethods.s_len_trim(word) == 0)
                    {
                        continue;
                    }

                    maxg = Convert.ToInt32(word);
                    break;
                }

            }
        }

        public static void pgma_read_test(string input_name)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMA_READ_TEST tests the ASCII PGM read routines.
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
            //    Input, string INPUT_NAME, the name of the file.
            //
        {
            int[] g = new int[1];
            int maxg = 0;
            int xsize = 0;
            int ysize = 0;
            int inputIndex = 0;

            //
            //  Read the data.
            //
            pgma_read(input_name, ref xsize, ref ysize, ref maxg, ref g);

            Console.WriteLine("");
            Console.WriteLine("PGMA_READ_TEST:");
            Console.WriteLine("  PGMA_READ was able to read \"" + input_name + "\".");
            //
            //  Check the data.
            //
            pgma_check_data(xsize, ysize, maxg, ref g);

            Console.WriteLine("");
            Console.WriteLine("PGMA_READ_TEST:");
            Console.WriteLine("  PGMA_CHECK_DATA has approved the data from the file.");
        }

        public static void pgma_write(string output_name, int xsize, int ysize, int[] g)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMA_WRITE writes the header and data for an ASCII PGM file.
            // 
            //  Example:
            //
            //    P2
            //    # feep.pgm
            //    24 7
            //    15
            //    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
            //    0  3  3  3  3  0  0  7  7  7  7  0  0 11 11 11 11  0  0 15 15 15 15  0
            //    0  3  0  0  0  0  0  7  0  0  0  0  0 11  0  0  0  0  0 15  0  0 15  0
            //    0  3  3  3  0  0  0  7  7  7  0  0  0 11 11 11  0  0  0 15 15 15 15  0
            //    0  3  0  0  0  0  0  7  0  0  0  0  0 11  0  0  0  0  0 15  0  0  0  0
            //    0  3  0  0  0  0  0  7  7  7  7  0  0 11 11 11 11  0  0 15  0  0  0  0
            //    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
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
            //    Input, string OUTPUT_NAME, the name of the file.
            //
            //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
            //
            //    Input, int *G, the array of XSIZE by YSIZE data values.
            //
        {
            List<string> output = new List<string>();
            int i;
            int[] indexg;
            int j;
            int maxg;
            //
            //  Compute the maximum.
            //
            maxg = 0;
            indexg = g;
            int gIndex = 0;

            for (j = 0; j < ysize; j++)
            {
                for (i = 0; i < xsize; i++)
                {
                    if (maxg < indexg[gIndex])
                    {
                        maxg = indexg[gIndex];
                    }

                    gIndex++;

                }
            }

            //
            //  Write the header.
            //
            pgma_write_header(ref output, output_name, xsize, ysize, maxg);
            //
            //  Write the data.
            //
            pgma_write_data(ref output, xsize, ysize, g);
            //
            //  Open the output file.
            //
            try
            {
                File.WriteAllLines(output_name, output);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("PGMA_WRITE - Fatal error!");
                Console.WriteLine("  Cannot open the output file \"" + output_name + "\".");
            }
        }

        public static void pgma_write_data(ref List<string> output, int xsize, int ysize, int[] g)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMA_WRITE_DATA writes the data for an ASCII PGM file.
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
            //    Input, ofstream &OUTPUT, a pointer to the file.
            //
            //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
            //
            //    Input, int *G, the array of XSIZE by YSIZE data.
            //
        {
            int i;
            int[] indexg;
            int j;
            int numval;

            indexg = g;
            numval = 0;
            int gIndex = 0;

            string cout = "";
            for (j = 0; j < ysize; j++)
            {
                for (i = 0; i < xsize; i++)
                {
                    cout += indexg[gIndex].ToString();
                    numval++;
                    gIndex++;

                    if (numval % 12 == 0 || i == xsize - 1 || numval == xsize * ysize)
                    {
                        output.Add(cout);
                        cout = "";
                    }
                    else
                    {
                        cout += " ";
                    }

                }
            }
        }

        public static void pgma_write_header(ref List<string> output, string output_name, int xsize,
                int ysize, int maxg)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMA_WRITE_HEADER writes the header of an ASCII PGM file.
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
            //    Input, ofstream &OUTPUT, a pointer to the file.
            //
            //    Input, string OUTPUT_NAME, the name of the file.
            //
            //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
            //
            //    Input, int MAXG, the maximum gray value.
            //
        {
            output.Add("P2");
            output.Add("# " + output_name + " created by PGMA_IO::PGMA_WRITE.");
            output.Add(xsize + "  " + ysize + "");
            output.Add(maxg + "");
        }

        public static void pgma_write_test(string output_name)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMA_WRITE_TEST tests the ASCII PGM write routines.
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
            //    Input, string OUTPUT_NAME, the name of the file.
            //
        {
            int[] g;
            int xsize;
            int ysize;

            xsize = 300;
            ysize = 300;
            //
            //  Allocate memory.
            //
            g = new int[xsize * ysize];
            //
            //  Set the data.
            //
            pgma_example(xsize, ysize, ref g);
            //
            //  Write the data to the file.
            //
            pgma_write(output_name, xsize, ysize, g);
        }
    }
}