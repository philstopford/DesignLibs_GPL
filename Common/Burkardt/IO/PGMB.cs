using System;
using System.Collections.Generic;
using System.IO;
using System.Security.Principal;
using Burkardt.Types;

namespace Burkardt.IO
{
    public static class PGMB
    {
        public static bool pgmb_check_data(int xsize, int ysize, ref int maxg,
                ref int[] g)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMB_CHECK_DATA checks the data for a binary portable gray map file.
            //
            //  Discussion:
            //
            //    XSIZE and YSIZE must be positive, the pointers must not be null,
            //    and the data must be nonnegative and no greater than MAXG.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    15 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
            //
            //    Input, unsigned char MAXG, the maximum gray value.
            //
            //    Input, unsigned char *G, the array of XSIZE by YSIZE data values.
            //
            //    Output, bool PGMB_CHECK_DATA, is
            //    true, if an error was detected, or
            //    false, if the data was legal.
            //
        {
            int i;
            int j;

            if (xsize <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("PGMB_CHECK_DATA: Error!");
                Console.WriteLine("  xsize <= 0.");
                Console.WriteLine("  xsize = " + xsize + "");
                return true;
            }

            if (ysize <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("PGMB_CHECK_DATA: Error!");
                Console.WriteLine("  ysize <= 0.");
                Console.WriteLine("  ysize = " + ysize + "");
                return true;
            }

            if (g == null)
            {
                Console.WriteLine("");
                Console.WriteLine("PGMB_CHECK_DATA: Error!");
                Console.WriteLine("  Null pointer to g.");
                return true;
            }

            int index = 0;
            for (j = 0; j < ysize; j++)
            {
                for (i = 0; i < xsize; i++)
                {
                    if (maxg < g[index])
                    {
                        Console.WriteLine("");
                        Console.WriteLine("PGMB_CHECK_DATA - Fatal error!");
                        Console.WriteLine("  Data exceeds MAXG = " + (int)maxg + "");
                        Console.WriteLine("  G(" + i + "," + j + ")=" + (int)(g[index]) + "");
                        return true;
                    }

                    index++;
                }
            }

            return false;
        }

        public static bool pgmb_example(int xsize, int ysize, ref int[] g)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMB_EXAMPLE sets up some data for a binary portable gray map file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    15 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
            //
            //    Output, unsigned char *G, the array of XSIZE by YSIZE gray values.
            //
            //    Output, bool PGMB_EXAMPLE, is
            //    false, if no error occurred,
            //    true, if an error occurred.
            //
        {
            int i;
            int indexg;
            int j;
            int periods = 3;
            float x;
            float y;

            indexg = 0;

            for (i = 0; i < ysize; i++)
            {
                y = ((float)(2 * i)) / ((float)(ysize - 1)) - 1.0f;

                for (j = 0; j < xsize; j++)
                {
                    x = (float)(2.0 * Math.PI * (float)(periods * j)) / ((float)(xsize - 1));

                    g[indexg] = (int)(20.0 * (Math.Sin(x) - y + 2));

                    indexg = indexg + 1;
                }
            }

            return false;
        }

        public static bool pgmb_read(string input_name, ref int xsize, ref int ysize,
                ref int maxg, ref int[] g)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMB_READ reads the header and data from a binary portable gray map file.
            // 
            //  Discussion:
            //
            //    Thanks to Jonas Schwertfeger for pointing out that, especially on 
            //    Microsoft Windows systems, a binary file needs to be opened as a 
            //    binary file!
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
            //    Input, string INPUT_NAME, the name of the file containing the binary
            //    portable gray map data.
            //
            //    Output, int &XSIZE, &YSIZE, the number of rows and columns of data.
            //
            //    Output, unsigned char &MAXG, the maximum gray value.
            //
            //    Output, unsigned char **G, the array of XSIZE by YSIZE data values.
            //
            //    Output, bool PGMB_READ, is true if an error occurred.
            //
        {
            bool error;
            string[] input;
            int inputIndex = -1;

            try
            {
                input = File.ReadAllLines(input_name);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("PGMB_READ: Fatal error!");
                Console.WriteLine("  Cannot open the input file " + input_name + "");
                return true;
            }

            //
            //  Read the header.
            //
            error = pgmb_read_header(input, ref inputIndex, ref xsize, ref ysize, ref maxg);

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("PGMB_READ: Fatal error!");
                Console.WriteLine("  PGMB_READ_HEADER failed.");
                return true;
            }

            //
            //  Allocate storage for the data.
            //
            g = new int[xsize * ysize];
            //
            //  Read the data.
            //
            error = pgmb_read_data(input, ref inputIndex, xsize, ysize, ref g);

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("PGMB_READ: Fatal error!");
                Console.WriteLine("  PGMB_READ_DATA failed.");
                return true;
            }

            return false;
        }

        public static bool pgmb_read_data(string[] input, ref int inputIndex, int xsize, int ysize,
                ref int[] g)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMB_READ_DATA reads the data in a binary portable gray map file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 December 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, ifstream &INPUT, a pointer to the file containing the binary
            //    portable gray map data.
            //
            //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
            //
            //    Output, unsigned char *G, the array of XSIZE by YSIZE data values.
            //
            //    Output, bool PGMB_READ_DATA, is true if an error occurred.
            //
        {
            int c;
            int i;
            int j;

            c = 0;            
            
            for (j = 0; j < ysize; j++)
            {
                for (i = 0; i < xsize; i++)
                {
                    g[c] = Convert.ToInt32(input[inputIndex]);
                    inputIndex++;
                    c++;
                }
            }

            return false;
        }

        public static bool pgmb_read_header(string[] input, ref int inputIndex, ref int xsize, ref int ysize,
                ref int maxg)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMB_READ_HEADER reads the header of a binary portable gray map file.
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
            //    Input, ifstream &INPUT, a pointer to the file containing the binary
            //    portable gray map data.
            //
            //    Output, int &XSIZE, &YSIZE, the number of rows and columns of data.
            //
            //    Output, unsigned char &MAXG, the maximum gray value.
            //
            //    Output, bool PGMB_READ_HEADER, is true if an error occurred.
            //
        {
            int fred;
            string line = "";
            string rest = "";
            int step;
            string word = "";

            step = 0;
            inputIndex = 0;

            while (true)
            {
                line = input[inputIndex];
                inputIndex++;

                if (line[0] == '#')
                {
                    continue;
                }

                if (step == 0)
                {
                    typeMethods.s_word_extract_first(line, ref word, ref rest);

                    if (typeMethods.s_len_trim(word) <= 0)
                    {
                        continue;
                    }

                    if (!typeMethods.s_eqi(word, "P5"))
                    {
                        Console.WriteLine("");
                        Console.WriteLine("PGMB_READ_HEADER - Fatal error.");
                        Console.WriteLine("  Bad magic number = \"" + word + "\".");
                        return true;
                    }

                    line = rest;
                    step = 1;
                }

                if (step == 1)
                {
                    typeMethods.s_word_extract_first(line, ref word, ref rest);

                    if (typeMethods.s_len_trim(word) <= 0)
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

                    if (typeMethods.s_len_trim(word) <= 0)
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

                    if (typeMethods.s_len_trim(word) <= 0)
                    {
                        continue;
                    }

                    fred = Convert.ToInt32(word);
                    maxg = fred;
                    line = rest;
                    break;
                }
            }

            return false;
        }

        public static bool pgmb_read_test(string input_name)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMB_READ_TEST tests the binary portable gray map read routines.
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
            //    Input, string INPUT_NAME, the name of the file containing the binary
            //    portable gray map data.
            //
            //    Output, bool PGMB_READ_TEST, is true if an error occurred.
            //
        {
            bool error;
            int[] g = null;
            int maxg = 0;
            int xsize = 0;
            int ysize = 0;
            //
            //  Read the data.
            //
            error = pgmb_read(input_name, ref xsize, ref ysize, ref maxg, ref g);

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("PGMB_READ_TEST: Fatal error!");
                Console.WriteLine("  PGMB_READ failed.");
                return true;
            }

            //
            //  Check the data.
            //
            error = pgmb_check_data(xsize, ysize, ref maxg, ref g);

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("  PGMB_CHECK_DATA reports bad data from the file.");
                return true;
            }

            Console.WriteLine("");
            Console.WriteLine("  PGMB_CHECK_DATA passes the data from the file.");

            return false;
        }

        public static bool pgmb_write(string output_name, int xsize, int ysize, int[] g)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMB_WRITE writes the header and data for a binary portable gray map file.
            //
            //  Discussion:
            //
            //    Thanks to Jonas Schwertfeger for pointing out that, especially on Microsoft
            //    Windows systems, a binary file needs to be opened as a binary file!
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            // 
            //    02 April 2005
            // 
            //  Author:
            // 
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string OUTPUT_NAME, the name of the file to contain the binary
            //    portable gray map data.
            //
            //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
            //
            //    Input, int *G, the array of XSIZE by YSIZE data values.
            //
            //    Output, bool PGMB_WRITE, is true if an error occurred.
            //
        {
            bool error;
            List<string> output = new List<string>();
            int i;
            int indexg;
            int j;
            int maxg;
            //
            //  Determine the maximum gray value.
            //
            maxg = 0;
            indexg = 0;

            for (i = 0; i < xsize; i++)
            {
                for (j = 0; j < ysize; j++)
                {
                    if (maxg < g[indexg])
                    {
                        maxg = g[indexg];
                    }

                    indexg = indexg + 1;
                }
            }

            //
            //  Write the header.
            //
            error = pgmb_write_header(ref output, xsize, ysize, maxg);

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("PGMB_WRITE: Fatal error!");
                Console.WriteLine("  PGMB_WRITE_HEADER failed.");
                return true;
            }

            //
            //  Write the data.
            //
            error = pgmb_write_data(ref output, xsize, ysize, g);

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("PGMB_WRITE: Fatal error!");
                Console.WriteLine("  PGMB_WRITE_DATA failed.");
                return true;
            }

            //
            //  Open the file.
            //
            try
            {
                File.WriteAllLines(output_name, output);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("PGMB_WRITE: Fatal error!");
                Console.WriteLine("  Cannot open the output file " + output_name + "");
                return true;
            }


            return false;
        }

        public static bool pgmb_write_data(ref List<string> output, int xsize, int ysize,
                int[] g)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMB_WRITE_DATA writes the data for a binary portable gray map file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, ofstream &OUTPUT, a pointer to the file to contain the binary
            //    portable gray map data.
            //
            //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
            //
            //    Input, unsigned char *G, the array of XSIZE by YSIZE data values.
            //
            //    Output, bool PGMB_WRITE_DATA, is true if an error occurred.
            //
        {
            int i;
            int j;
            int index = 0;

            for (j = 0; j < ysize; j++)
            {
                for (i = 0; i < xsize; i++)
                {
                    output.Add(g[index].ToString());
                    index++;
                }
            }

            return false;
        }

        public static bool pgmb_write_header(ref List<string> output, int xsize, int ysize,
                int maxg)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMB_WRITE_HEADER writes the header of a binary portable gray map file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, ofstream &OUTPUT, a pointer to the file to contain the binary
            //    portable gray map data.
            //
            //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
            //
            //    Input, unsigned char MAXG, the maximum gray value.
            //
            //    Output, bool PGMB_WRITE_HEADER, is true if an error occurred.
            //
        {
            output.Add("P5" + " "
                            + xsize + " "
                            + ysize + " "
                            + (int)maxg + "");

            return false;
        }

        public static bool pgmb_write_test(string output_name)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PGMB_WRITE_TEST tests the binary portable gray map write routines.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string OUTPUT_NAME, the name of the file to contain the binary
            //    portable gray map data.
            //
            //    Output, bool PGMB_WRITE_TEST, is true if an error occurred.
            //
        {
            bool error;
            int[] g;
            int xsize;
            int ysize;
            //
            //  Set the data.
            //
            xsize = 300;
            ysize = 200;

            g = new int[xsize * ysize];

            error = pgmb_example(xsize, ysize, ref g);

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("PGMB_WRITE_TEST: Fatal error!");
                Console.WriteLine("  PGM_EXAMPLE failed.");
                return true;
            }

            error = pgmb_write(output_name, xsize, ysize, g);

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("PGMB_WRITE_TEST: Fatal error!");
                Console.WriteLine("  PGMB_WRITE failed.");
                return true;
            }

            return false;
        }
    }
}