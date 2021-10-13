using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using Burkardt.Types;
using MiscUtil.Conversion;
using MiscUtil.IO;

namespace Burkardt.IO
{
    public static class PBMB
    {
        public static bool pbmb_check_data(int xsize, int ysize, int[] barray)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PBMB_CHECK_DATA checks the data for a binary portable bit map file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
            //
            //    Input, int[] BARRAY, the array of XSIZE by YSIZE bits.
            //
            //    Output, bool PBMB_CHECK_DATA, is true if an error occurred.
            //
        {
            int i;
            int indexb;
            int j;

            if (xsize <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("PBMB_CHECK_DATA - Fatal error!");
                Console.WriteLine("  XSIZE <= 0");
                Console.WriteLine("  XSIZE = " + xsize + "");
                return true;
            }

            if (ysize <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("PBMB_CHECK_DATA - Fatal error!");
                Console.WriteLine("  YSIZE <= 0");
                Console.WriteLine("  YSIZE = " + ysize + "");
                return true;
            }

            indexb = 0;

            for (j = 0; j < ysize; j++)
            {
                for (i = 0; i < xsize; i++)
                {
                    if (barray[indexb] != 0 && barray[indexb] != 1)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("PBMB_CHECK_DATA - Fatal error!");
                        Console.WriteLine("  b(" + i + "," + j + ") = "
                                          + barray[indexb] + ".");
                        return true;
                    }

                    indexb = indexb + 1;
                }
            }

            return false;
        }

        public static bool pbmb_example(int xsize, int ysize, ref int[] barray)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PBMB_EXAMPLE sets up some sample PBMB data.
            //
            //  Discussion:
            //
            //    The data represents an ellipse.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
            //    Values of 200 would be reasonable.
            //
            //    Output, int[] BARRAY, the array of XSIZE by YSIZE data values.
            //
            //    Output, bool PBMB_EXAMPLE, is true if an error occurred.
            //
        {
            int i;
            int indexb;
            int j;
            float r;
            float test;
            float x;
            float xc;
            float y;
            float yc;

            indexb = 0;
            if (xsize < ysize)
            {
                r = (float) xsize / 3.0f;
            }
            else
            {
                r = (float) ysize / 3.0f;
            }

            xc = (xsize) / 2.0f;
            yc = (ysize) / 2.0f;

            for (i = 0; i < ysize; i++)
            {
                y = (float) i;
                for (j = 0; j < xsize; j++)
                {
                    x = (float) j;
                    test = r - (float) Math.Sqrt((x - xc) * (x - xc)
                                                 + 0.75f * (y - yc) * (y - yc));
                    if (Math.Abs(test) <= 3.0)
                    {
                        barray[indexb] = 1;
                    }
                    else
                    {
                        barray[indexb] = 0;
                    }

                    indexb = indexb + 1;
                }
            }

            return false;
        }

        public static bool pbmb_read(string input_name, ref int xsize, ref int ysize, ref int[] barray)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PBMB_READ reads the header and data from a binary portable bit map file.
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
            //    22 July 2011
            // 
            //  Author:
            // 
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string INPUT_NAME, the name of the file containing the binary
            //    portable bit map data.
            //
            //    Output, int &XSIZE, &YSIZE, the number of rows and columns of data.
            //
            //    Output, int[] *BARRAY, the array of XSIZE by YSIZE data values.
            //
            //    Output, bool PBMB_READ, is true if an error occurred.
            //
        {
            bool error;
            string[] input;

            Stream file_in_s;
            EndianBinaryReader file_in;

            try
            {
                file_in_s = File.OpenRead(input_name);
                file_in = new EndianBinaryReader(EndianBitConverter.Little, file_in_s);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("PBMB_READ: Fatal error!");
                Console.WriteLine("  Cannot open the input file " + input_name + "");
                return true;
            }

            //
            //  Read the header.
            //
            error = pbmb_read_header(ref file_in, ref xsize, ref ysize);

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("PBMB_READ: Fatal error!");
                Console.WriteLine("  PBMB_READ_HEADER failed.");
                return true;
            }

            //
            //  Allocate storage for the data.
            //
            barray = new int [xsize * ysize];
            //
            //  Read the data.
            //
            error = pbmb_read_data(ref file_in, xsize, ysize, ref barray);

            file_in.Close();
            
            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("PBMB_READ: Fatal error!");
                Console.WriteLine("  PBMB_READ_DATA failed.");
                return true;
            }

            return false;
        }

        public static bool pbmb_read_data(ref EndianBinaryReader br, int xsize, int ysize, ref int[] barray)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PBMB_READ_DATA reads the data in a binary portable bit map file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 October 2017
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, ifstream &INPUT, a pointer to the file containing the binary
            //    portable bit map data.
            //
            //    Input, int &XSIZE, &YSIZE, the number of rows and columns of data.
            //
            //    Output, int[] BARRAY, the array of XSIZE by YSIZE data values.
            //
            //    Output, bool PBMB_READ_DATA, is true if an error occurred.
            //
        {
            int bit;
            short c = 0;
            ushort c2 = 0;
            int i;
            int j;
            int k;
            int numbyte;

            int indexb = 0;
            numbyte = 0;

            for ( j = 0; j < ysize; j++ )
            {
                for ( i = 0; i < xsize; i++ )
                {
                    if ( i%8 == 0 )
                    {
                        try
                        {
                            c = br.ReadByte();
                            c2 = ( ushort ) c;
                        }
                        catch (Exception e)
                        {
                            /*
                            Console.WriteLine();
                            Console.WriteLine("PBMB_CHECK_DATA - Fatal error!");
                            Console.WriteLine("  Failed reading byte " + numbyte);
                            */
                        }
                        numbyte = numbyte + 1;
                    }

                    k = 7 - ( i % 8 );
                    bit = ( c2 >> k ) & 1;

                    barray[indexb] = bit;
                    indexb++;
                }
            }
            return false;
        }

        public static bool pbmb_read_header(ref EndianBinaryReader br, ref int xsize, ref int ysize)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PBMB_READ_HEADER reads the header of a binary portable bit map file.
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
            //    portable bit map data.
            //
            //    Output, int &XSIZE, &YSIZE, the number of rows and columns of data.
            //
            //    Output, bool PBMB_READ_HEADER, is true if an error occurred.
            //
        {
            string line;
            string rest = "";
            int step;
            string word = "";

            step = 0;
            int index = 0;

            
            while (true)
            {
                try
                {
                    // Read null-terminated string.
                    List<byte> strBytes = new List<byte>();
                    int b;
                    while ((b = br.ReadByte()) != 0x00)
                    {
                        strBytes.Add((byte) b);
                    }

                    line = Encoding.ASCII.GetString(strBytes.ToArray());
                }
                catch
                {
                    Console.WriteLine("");
                    Console.WriteLine("PBMB_READ_HEADER - Fatal error!");
                    Console.WriteLine("  End of file.");
                    return true;
                }

                if (line.StartsWith('#'))
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

                    if (!typeMethods.s_eqi(word, "P4"))
                    {
                        Console.WriteLine("");
                        Console.WriteLine("PBMB_READ_HEADER - Fatal error.");
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
                    break;
                }

            }
            return false;
        }

        public static bool pbmb_read_test(string input_name)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PBMB_READ_TEST tests the binary portable bit map read routines.
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
            //    portable bit map data.
            //
            //    Output, bool PBMB_READ_TEST, is true if an error occurred.
            //
        {
            int[] barray = null;
            bool error;
            int xsize = 0;
            int ysize = 0;

            //
            //  Read the data.
            //
            error = pbmb_read(input_name, ref xsize, ref ysize, ref barray);

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("PBMB_READ_TEST: Fatal error!");
                Console.WriteLine("  PBMB_READ failed.");
                return true;
            }

            //
            //  Check the data.
            //
            error = pbmb_check_data(xsize, ysize, barray);


            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("  PBMB_CHECK_DATA reports bad data from the file.");
                return true;
            }

            Console.WriteLine("");
            Console.WriteLine("  PBMB_CHECK_DATA passes the data from the file.");

            return false;
        }

        public static bool pbmb_write(string output_name, int xsize, int ysize, int[] barray)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PBMB_WRITE writes the header and data for a binary portable bit map file.
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
            //    portable bit map data.
            //
            //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
            //
            //    Input, int[] BARRAY, the array of XSIZE by YSIZE data values.
            //
            //    Output, bool PBMB_WRITE, is true if an error occurred.
            //
        {
            bool error = false;
            Stream file_out_s;
            EndianBinaryWriter file_out = null;

            try
            {
                file_out_s = File.OpenWrite(output_name);
                file_out = new EndianBinaryWriter(EndianBitConverter.Little, file_out_s);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("PBMB_WRITE: Fatal error!");
                Console.WriteLine("  Cannot open the output file " + output_name + "");
                error = true;
            }

            if (!error)
            {
                //
                //  Write the header.
                //
                error = pbmb_write_header(ref file_out, xsize, ysize);

                if (error)
                {
                    Console.WriteLine("");
                    Console.WriteLine("PBMB_WRITE: Fatal error!");
                    Console.WriteLine("  PBMB_WRITE_HEADER failed.");
                }
                else
                {
                    //
                    //  Write the data.
                    //
                    error = pbmb_write_data(ref file_out, xsize, ysize, barray);

                    if (error)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("PBMB_WRITE: Fatal error!");
                        Console.WriteLine("  PBMB_WRITE_DATA failed.");
                    }
                }
                file_out.Close();
            }

            return error;
        }

        public static bool pbmb_write_data(ref EndianBinaryWriter file_out, int xsize, int ysize, int[] barray)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PBMB_WRITE_DATA writes the data for a binary portable bit map file.
            //
            //  Discussion:
            //
            //    Thanks to Andreas Wagner for correcting the computation of "bit".
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 October 2017
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, List<string> &OUTPUT, a pointer to the file to contain the binary
            //    portable bit map data.
            //
            //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
            //
            //    Input, int[] BARRAY, the array of XSIZE by YSIZE data values.
            //
            //    Output, bool PBMB_WRITE_DATA, is true if an error occurred.
            //
        {
            int bit;
            byte c;
            int i;
            int indexb;
            int j;
            int k;

            indexb = 0;
            c = 0;

            for ( j = 0; j < ysize; j++ )
            {
                for ( i = 0; i < xsize; i++ )
                {
                    k = 7 - ( i % 8 );
                    bit = ( barray[indexb] ) & 1;
                    c = (byte)(c | ( bit << k ));

                    indexb++;

                    if ( (i+1)%8 == 0 || i == ( xsize - 1 ) )
                    {
                        file_out.Write(c);
                        c = 0;
                    }
                }
            }
            return false;
        }

        public static bool pbmb_write_header(ref EndianBinaryWriter file_out, int xsize, int ysize)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PBMB_WRITE_HEADER writes the header of a binary portable bit map file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, List<string> &OUTPUT, a pointer to the file to contain the binary
            //    portable bit map data.
            //
            //    Input, int XSIZE, YSIZE, the number of rows and columns of data.
            //
            //    Output, bool PBMB_WRITE_HEADER, is true if an error occurred.
            //
        {
            file_out.Write(("P4" + " "
                            + xsize + " "
                            + ysize + "").ToCharArray());

            return false;
        }

        public static bool pbmb_write_test(string output_name)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PBMB_WRITE_TEST tests the binary portable bit map write routines.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string OUTPUT_NAME, the name of the file to contain the binary
            //    portable bit map data.
            //
            //    Output, bool PBMB_WRITE_TEST, is true if an error occurred.
            //
        {
            int[] barray;
            bool error;
            int xsize;
            int ysize;
            //
            //  Set the data.
            //
            xsize = 250;
            ysize = 150;

            barray = new int [xsize * ysize];

            error = pbmb_example(xsize, ysize, ref barray);

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("PBMB_WRITE_TEST: Fatal error!");
                Console.WriteLine("  PBM_EXAMPLE failed.");
                return true;
            }

            error = pbmb_write(output_name, xsize, ysize, barray);

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("PBMB_WRITE_TEST: Fatal error!");
                Console.WriteLine("  PBMB_WRITE failed.");
                return true;
            }

            return false;
        }
    }
}