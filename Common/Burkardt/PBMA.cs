using System;
using System.Collections.Generic;
using System.IO;

namespace Burkardt
{
    public static class PBMA
    {
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
            List<string> file_out = new List<string>();
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
            int i;
            int[] indexb;
            int j;
            int numval;

            indexb = b;
            numval = 0;
            string cout = "";
            int bIndex = 0;

            for (j = 0; j < ysize; j++)
            {
                for (i = 0; i < xsize; i++)
                {
                    cout += indexb[bIndex] + " ";
                    numval = numval + 1;
                    bIndex++;

                    if (numval % 12 == 0 || i == xsize - 1 || numval == xsize * ysize)
                    {
                        file_out.Add(cout);
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
    }
}