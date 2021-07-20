using System;
using System.Collections.Generic;
using System.IO;

namespace Burkardt
{
    public static class PPM_ASCII
    {
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
            int[] b_index;
            bool error;
            List<string> file_out = new List<string>();
            int[] g_index;
            int i;
            int j;
            int[] r_index;
            int rgb_max;

            //
            //  Compute the maximum.
            //
            rgb_max = 0;
            r_index = r;
            g_index = g;
            b_index = b;

            int index = 0;
            for (j = 0; j < ysize; j++)
            {
                for (i = 0; i < xsize; i++)
                {
                    if (rgb_max < r_index[index])
                    {
                        rgb_max = r_index[index];
                    }

                    if (rgb_max < g_index[index])
                    {
                        rgb_max = g_index[index];
                    }

                    if (rgb_max < b_index[index])
                    {
                        rgb_max = b_index[index];
                    }

                    index++;
                }
            }

            //
            //  Write the header.
            //
            error = ppma_write_header(ref file_out, file_out_name, xsize, ysize, rgb_max);

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("PPMA_WRITE - Fatal error!");
                Console.WriteLine("  PPMA_WRITE_HEADER failed.");
                return true;
            }

            //
            //  Write the data.
            //
            error = ppma_write_data(ref file_out, xsize, ysize, r, g, b);

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("PPMA_WRITE - Fatal error!");
                Console.WriteLine("  PPMA_WRITE_DATA failed.");
                return true;
            }

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
            int[] b_index;
            int[] g_index;
            int i;
            int j;
            int[] r_index;
            int rgb_num;

            r_index = r;
            g_index = g;
            b_index = b;
            rgb_num = 0;

            int index = 0;
            for (j = 0; j < ysize; j++)
            {
                string tmp = "";
                for (i = 0; i < xsize; i++)
                {
                    tmp += r_index[index] + " " + g_index[index] + " " + b_index[index];
                    rgb_num = rgb_num + 3;

                    if (rgb_num % 12 == 0 || i == xsize - 1 || rgb_num == 3 * xsize * ysize)
                    {
                        file_out.Add(tmp);
                        tmp = "";
                    }
                    else
                    {
                        tmp += " ";
                    }
                }

                index++;
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
    }
}