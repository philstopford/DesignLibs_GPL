using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt.IO
{
    public static class XYZ
    {
        public static void xyz_data_print(int point_num, double[] xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZ_DATA_PRINT prints the data for an XYZ file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    03 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, double XY[3*POINT_NUM], the arrays of coordinate data.
            //
        {
            int j;

            Console.WriteLine("");
            for (j = 0; j < point_num; j++)
            {
                Console.WriteLine(xyz[0 + j * 3].ToString().PadLeft(10) + "  "
                                                                        + xyz[1 + j * 3].ToString().PadLeft(10) + "  "
                                                                        + xyz[2 + j * 3].ToString().PadLeft(10) + "");

            }
        }

        public static void xyz_data_read(string input_filename, int point_num, ref double[] xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZ_DATA_READ reads the data in an XYZ file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 October 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string INPUT_FILENAME, the name of the file.
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Output, double XYZ[3*POINT_NUM], the point coordinates.
            //
        {
            bool error;
            int i;
            string[] input;
            int j;
            string text;
            double[] temp = new double[3];

            try
            {
                input = File.ReadAllLines(input_filename);
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("XYZ_DATA_READ - Fatal error!");
                Console.WriteLine("  Cannot open the input file \"" + input_filename + "\".");
                return;
            }

            j = 0;

            int index = 0;

            while (j < point_num)
            {
                try
                {
                    text = input[index];
                    index++;
                }
                catch (Exception e)
                {
                    break;
                }

                if (text[0] == '#' || typeMethods.s_len_trim(text) == 0)
                {
                    continue;
                }

                //
                //  Extract two real numbers.
                //
                r8vec r = typeMethods.s_to_r8vec(text, 3);

                if (r.error)
                {
                    Console.WriteLine("");
                    Console.WriteLine("XYZ_DATA_READ - Fatal error!");
                    Console.WriteLine("  S_TO_R8VEC returned error flag.");
                    return;
                }

                xyz[0 + j * 3] = r.rvec[0];
                xyz[1 + j * 3] = r.rvec[1];
                xyz[2 + j * 3] = r.rvec[2];
                j = j + 1;
            }
        }

        public static void xyz_data_write(ref List<string> output_unit, int point_num, double[] xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZ_DATA_WRITE writes the data for an XYZ file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, ofstream &OUTPUT_UNIT, a pointer to the file.
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, double XY[3*POINT_NUM], the arrays of coordinate data.
            //
        {
            int j;

            for (j = 0; j < point_num; j++)
            {
                output_unit.Add(xyz[0 + j * 3].ToString().PadLeft(10) + "  "
                                                                      + xyz[1 + j * 3].ToString().PadLeft(10) + "  "
                                                                      + xyz[2 + j * 3].ToString().PadLeft(10) + "");

            }
        }

        public static void xyz_example(int point_num, ref double[] xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZ_EXAMPLE sets up data suitable for an XYZ file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Output, double XYZ[3*POINT_NUM], the point coordinates.
            //
        {
            double a;
            int i;
            double r;
            
            double theta;
            double theta1;
            double theta2;

            r = 1.0;
            theta1 = 0.0;
            theta2 = 3.0 * 2.0 * Math.PI;
            a = 2.0 / (theta2 - theta1);

            for (i = 0; i < point_num; i++)
            {
                if (point_num == 1)
                {
                    theta = 0.5 * (theta1 + theta2);
                }
                else
                {
                    theta = ((double) (point_num - i - 1) * theta1
                             + (double) (i) * theta2)
                            / (double) (point_num - 1);
                }

                xyz[0 + i * 3] = r * Math.Cos(theta);
                xyz[1 + i * 3] = r * Math.Sin(theta);
                xyz[2 + i * 3] = a * theta;
            }

        }

        public static int xyz_example_size()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZ_EXAMPLE_SIZE sizes an example XYZ dataset.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, integer XYZ_EXAMPLE_SIZE, the number of points.
            //
        {
            int point_num;

            point_num = 101;

            return point_num;
        }

        public static void xyz_header_print(int point_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZ_HEADER_PRINT prints the header of an XYZ file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    03 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int POINT_NUM, the number of points.
            //
        {
            Console.WriteLine("");
            Console.WriteLine("  Number of points = " + point_num + "");

            return;
        }

        public static void xyz_header_read(string input_filename, ref int point_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZ_HEADER_READ reads the header of an XYZ file.
            //
            //  Discussion:
            //
            //    All we do here is count the number of lines that are not comments 
            //    and not blank.  Each such line is assumed to represent a single point 
            //    coordinate record.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 October 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string INPUT_FILENAME, the name of the file.
            //
            //    Output, int *POINT_NUM, the number of points.
            //
        {
            string[] input;
            string text;

            point_num = 0;

            try
            {
                input = File.ReadAllLines(input_filename);
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("XYZ_HEADER_READ - Fatal error!");
                Console.WriteLine("  Cannot open the input file \"" + input_filename + "\".");
                return;
            }

            int index = 0;

            while (true)
            {

                try
                {
                    text = input[index];
                    index++;
                }
                catch (Exception e)
                {
                    break;
                }

                if (text[0] == '#' || typeMethods.s_len_trim(text) == 0)
                {
                    continue;
                }

                point_num = point_num + 1;
            }

        }

        public static void xyz_header_write(string output_filename, ref List<string> output_unit,
                int point_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZ_HEADER_WRITE writes the header of an XYZ file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string OUTPUT_FILENAME, the name of the file.
            //
            //    Input, ofstream &OUTPUT_UNIT, a pointer to the file to contain the data.
            //
            //    Input, int POINT_NUM, the number of points.
            //
        {
            output_unit.Add("#  " + output_filename + "");
            output_unit.Add("#  created by xyz_io::xyz_header_write.C");
            output_unit.Add("#");
            output_unit.Add("#  Number of points = " + point_num + "");
            output_unit.Add("#");
        }

        public static void xyz_read(string input_filename, ref int point_num, ref double[] xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZ_READ reads the header and data from an XYZ file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            // 
            //    05 January 2009
            // 
            //  Author:
            // 
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string INPUT_FILENAME, the name of the file.
            //
            //    Output, int *POINT_NUM, the number of points.
            //
            //    Output, double *XYZ[3*(*POINT_NUM)], the point coordinates.
            //
        {
            xyz_header_read(input_filename, ref point_num);

            xyz = new double[3 * (point_num)];

            xyz_data_read(input_filename, point_num, ref xyz);

        }

        public static void xyz_read_test(string input_filename)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZ_READ_TEST tests the XYZ file read routines.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string INPUT_FILENAME, the name of the file.
            //
        {
            int point_num = 0;
            double[] xyz = null;

            xyz_read(input_filename, ref point_num, ref xyz);

        }

        public static void xyz_write(string output_filename, int point_num, double[] xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZ_WRITE writes the header and data for an XYZ file.
            // 
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            // 
            //    01 January 2009
            // 
            //  Author:
            // 
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string OUTPUT_FILENAME, the name of the file.
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, double XYZ[3*POINT_NUM], the arrays of coordinate data.
            //
        {
            List<string> output_unit = new List<string>();

            //
            //  Write the header.
            //
            xyz_header_write(output_filename, ref output_unit, point_num);
            //
            //  Write the data.
            //
            xyz_data_write(ref output_unit, point_num, xyz);

            //
            //  Open the output file.
            //
            try
            {
                File.WriteAllLines(output_filename, output_unit);
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("XYZ_WRITE - Fatal error!");
                Console.WriteLine("  Cannot open the output file \"" + output_filename + "\".");
            }

        }

        public static void xyz_write_test(string output_filename)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZ_WRITE_TEST tests the XYZ write routines.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string OUTPUT_FILENAME, the name of the file to contain the data.
            //
        {
            int point_num;
            double[] xyz;
            //
            //  Set the data.
            //
            point_num = xyz_example_size();
            //
            //  Allocate memory.
            //
            xyz = new double[3 * point_num];
            //
            //  Set the data.
            //
            xyz_example(point_num, ref xyz);
            //
            //  Write the data to the file.
            //
            xyz_write(output_filename, point_num, xyz);

        }

        public static void xyzf_data_print(int point_num, int face_num,
                int face_data_num, int[] face_pointer, int[] face_data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZF_DATA_PRINT prints the data of an XYZF file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int FACE_NUM, the number of faces.
            //
            //    Input, int FACE_DATA_NUM, the number of face items.
            //
            //    Input, int FACE_POINTER[FACE_NUM+1], pointers to the
            //    first face item for each face.
            //
            //    Input, int FACE_DATA[FACE_DATA_NUM], indices
            //    of points that form faces.
            //
        {
            int i;
            int face;

            Console.WriteLine("");
            for (face = 0; face < face_num; face++)
            {
                Console.WriteLine("  " + face.ToString().PadLeft(4) + "  "
                                  + "  " + face_pointer[face].ToString().PadLeft(8)
                                  + "  " + (face_pointer[face + 1] - 1).ToString().PadLeft(8) + "");
            }

            Console.WriteLine("");
            for (face = 0; face < face_num; face++)
            {
                string cout = "  " + face.ToString().PadLeft(4) + "  ";
                for (i = face_pointer[face]; i < face_pointer[face + 1]; i++)
                {
                    cout += "  " + face_data[i].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }

        public static void xyzf_data_read(string input_filename, int face_num, int face_data_num,
                ref int[] face_pointer, ref int[] face_data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZF_DATA_READ reads the data in an XYZF file.
            //
            //  Discussion:
            //
            //    This routine assumes that the file contains exactly three kinds of
            //    records:
            //
            //    COMMENTS which begin with a '#' character in column 1;
            //    BLANKS which contain nothing but 'whitespace';
            //    FACE ITEMS, which are indices of points on a face.
            //
            //    The routine ignores comments and blanks and returns
            //    the number of face items.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 October 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string INPUT_FILENAME, the name of the input file.
            //
            //    Input, int FACE_NUM, the number of faces.
            //
            //    Input, int FACE_DATA_NUM, the number of face items.
            //
            //    Output, int FACE_POINTER[FACE_NUM+1], pointers to the
            //    first face item for each face.
            //
            //    Output, int FACE_DATA[FACE_DATA_NUM], the face items.
            //
        {
            int face;
            int ierror;
            int ilo;
            string[] input;
            int n;
            string text;

            try
            {
                input = File.ReadAllLines(input_filename);
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("XYZF_DATA_READ - Fatal error!");
                Console.WriteLine("  Cannot open the input file \"" + input_filename + "\".");
                return;
            }

            face = 0;
            face_pointer[0] = 0;

            int index = 0;

            while (face < face_num)
            {
                try
                {
                    text = input[index];
                    index++;
                }
                catch (Exception e)
                {
                    Console.WriteLine("");
                    Console.WriteLine("XYZF_DATA_READ - Fatal error!");
                    Console.WriteLine("  Unexpected end of information.");
                    return;
                }

                if (text[0] == '#' || typeMethods.s_len_trim(text) == 0)
                {
                    continue;
                }

                n = typeMethods.s_word_count(text);
                face_pointer[face + 1] = face_pointer[face] + n;

                ilo = face_pointer[face];

                i4vec r = typeMethods.s_to_i4vec(text, n);

                if (r.error)
                {
                    Console.WriteLine("");
                    Console.WriteLine("XYZF_DATA_READ - Fatal error!");
                    Console.WriteLine("  Error from S_TO_I4VEC.");
                    return;
                }

                for (int ix = 0; ix < r.ivec.Length; ix++)
                {
                    face_data[ilo + ix] = r.ivec[ix];
                }

                face = face + 1;
            }
        }

        public static void xyzf_data_write(ref List<string> output_unit, int point_num, int face_num,
                int face_data_num, int[] face_pointer, int[] face_data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZF_DATA_WRITE writes the data of an XYZF file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, ofstream &OUTPUT_UNIT, a pointer to the XY file.
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int FACE_NUM, the number of faces.
            //
            //    Input, int FACE_DATA_NUM, the number of face items.
            //
            //    Input, int FACE_POINTER[FACE_NUM+1], pointers to the
            //    first face item for each face.
            //
            //    Input, int FACE_DATA[FACE_DATA_NUM], indices
            //    of points that form faces.
            //
        {
            int i;
            int face;

            for (face = 0; face < face_num; face++)
            {
                string cout = "";
                for (i = face_pointer[face]; i < face_pointer[face + 1]; i++)
                {
                    cout += "  " + face_data[i];
                }

                output_unit.Add(cout);
            }
        }

        public static void xyzf_example(int point_num, int face_num, int face_data_num,
                ref double[] xyz, ref int[] face_pointer, ref int[] face_data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZF_EXAMPLE sets data suitable for a pair of XYZ and XYZF files.
            //
            //  Discussion:
            //
            //    There are 8 points.
            //    There are 6 faces.
            //    There are 24 face items.
            //
            //       8------7
            //      /|     /|
            //     / |    / |
            //    5------6  |
            //    |  4---|--3
            //    | /    | /
            //    |/     |/
            //    1------2
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int FACE_NUM, the number of faces.
            //
            //    Input, int FACE_DATA_NUM, the number of face items.
            //
            //    Output, double XYZ[3*POINT_NUM], the point coordinates.
            //
            //    Output, int FACE_POINTER[FACE_NUM+1], pointers to the
            //    first face item for each face.
            //
            //    Output, int FACE_DATA[FACE_DATA_NUM], indices
            //    of points that form faces.
            //
        {
            int FACE_DATA_NUM = 24;
            int FACE_NUM = 6;
            int POINT_NUM = 8;

            int[] FACE_DATA =
            {
                0, 3, 2, 1,
                1, 2, 6, 5,
                4, 5, 6, 7,
                4, 7, 3, 0,
                0, 1, 5, 4,
                2, 3, 7, 6
            };
            int[] FACE_POINTER = {0, 4, 8, 12, 16, 20, 24};
            double[] XYZ =
            {
                0.0, 0.0, 0.0,
                1.0, 0.0, 0.0,
                1.0, 1.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0,
                1.0, 0.0, 1.0,
                1.0, 1.0, 1.0,
                0.0, 1.0, 1.0
            };

            typeMethods.i4vec_copy(face_data_num, FACE_DATA, ref face_data);
            typeMethods.i4vec_copy(face_num + 1, FACE_POINTER, ref face_pointer);
            typeMethods.r8vec_copy(3 * point_num, XYZ, ref xyz);

        }

        public static void xyzf_example_size(ref int point_num, ref int face_num, ref int face_data_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZF_EXAMPLE_SIZE sizes the data to be created by XYZF_EXAMPLE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, int *POINT_NUM, the number of points.
            //
            //    Output, int *FACE_NUM, the number of faces.
            //
            //    Output, int *FACE_DATA_NUM, the number of face items.
            //
        {
            face_data_num = 24;
            face_num = 6;
            point_num = 8;
        }

        public static void xyzf_header_print(int point_num, int face_num, int face_data_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZF_HEADER_PRINT prints the header of an XYZF file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int FACE_NUM, the number of faces.
            //
            //    Input, int FACE_DATA_NUM, the number of face items.
            //
        {
            Console.WriteLine("");
            Console.WriteLine("  Number of points     = " + point_num + "");
            Console.WriteLine("  Number of faces      = " + face_num + "");
            Console.WriteLine("  Number of face items = " + face_data_num + "");

            return;
        }

        public static void xyzf_header_read(string input_filename, ref int face_num,
                ref int face_data_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZF_HEADER_READ determines the number of face items in an XYZF file.
            //
            //  Discussion:
            //
            //    This routine assumes that the file contains exactly three kinds of
            //    records:
            //
            //    COMMENTS which begin with a '#' character in column 1;
            //    BLANKS which contain nothing but 'whitespace';
            //    FACE ITEMS, which are indices of points on a face.
            //
            //    The routine ignores comments and blanks and returns
            //    the number of face items.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 October 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string INPUT_FILENAME, the name of the input file.
            //
            //    Output, int *FACE_NUM, the number of faces.
            //
            //    Output, int *FACE_DATA_NUM, the number of face items.
            //
        {
            int i;
            int i4_val;
            int ierror;
            string[] input;
            int length;
            int n;
            string text;

            face_data_num = 0;
            face_num = 0;

            try
            {
                input = File.ReadAllLines(input_filename);
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("XYZF_HEADER_READ - Fatal error!");
                Console.WriteLine("  Cannot open the input file \"" + input_filename + "\".");
                return;
            }

            int index = 0;

            for (;;)
            {
                try
                {
                    text = input[index];
                    index++;
                }
                catch (Exception e)
                {
                    break;
                }

                if (text[0] == '#' || typeMethods.s_len_trim(text) == 0)
                {
                    continue;
                }

                n = typeMethods.s_word_count(text);

                face_data_num = face_data_num + n;

                face_num = face_num + 1;
            }

        }

        public static void xyzf_header_write(string output_filename, ref List<string> output_unit,
                int point_num, int face_num, int face_data_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZF_HEADER_WRITE writes the header of an XYZF file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string OUTPUT_FILENAME, the name of the file.
            //
            //    Input, ofstream &OUTPUT_UNIT, a pointer to the file to contain the data.
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int FACE_NUM, the number of faces.
            //
            //    Input, int FACE_DATA_NUM, the number of face items.
            //
        {
            output_unit.Add("#  " + output_filename + "");
            output_unit.Add("#  created by xyz_io::xyzf_header_write.C");
            output_unit.Add("#");
            output_unit.Add("#  Number of points     = " + point_num + "");
            output_unit.Add("#  Number of faces      = " + face_num + "");
            output_unit.Add("#  Number of face items = " + face_data_num + "");
            output_unit.Add("#");
        }

        public static void xyzf_write(string output_filename, int point_num, int face_num,
                int face_data_num, int[] face_pointer, int[] face_data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZF_WRITE writes the header and data for an XYZF file.
            // 
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            // 
            //    07 January 2009
            // 
            //  Author:
            // 
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string OUTPUT_FILENAME, the name of the file.
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int FACE_NUM, the number of faces.
            //
            //    Input, int FACE_DATA_NUM, the number of face items.
            //
            //    Input, int FACE_POINTER[FACE_NUM+1], pointers to the
            //    first face item for each face.
            //
            //    Input, int FACE_DATA[FACE_DATA_NUM], indices
            //    of points that form faces.
            //
        {
            List<string> output_unit = new List<string>();

            //
            //  Write the header.
            //
            xyzf_header_write(output_filename, ref output_unit, point_num, face_num,
                face_data_num);
            //
            //  Write the data.
            //
            xyzf_data_write(ref output_unit, point_num, face_num, face_data_num, face_pointer,
                face_data);

            //
            //  Open the output file.
            //
            try
            {
                File.WriteAllLines(output_filename, output_unit);
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("XYZF_WRITE - Fatal error!");
                Console.WriteLine("  Cannot open the output file \"" + output_filename + "\".");
            }

        }

        public static void xyzl_data_print(int point_num, int line_num,
                int line_data_num, int[] line_pointer, int[] line_data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZL_DATA_PRINT prints the data of an XYZL file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int LINE_NUM, the number of lines.
            //
            //    Input, int LINE_DATA_NUM, the number of line items.
            //
            //    Input, int LINE_POINTER[LINE_NUM+1], pointers to the
            //    first line item for each line.
            //
            //    Input, int LINE_DATA[LINE_DATA_NUM], indices
            //    of points that form lines.
            //
        {
            int i;
            int line;

            Console.WriteLine("");
            for (line = 0; line < line_num; line++)
            {
                Console.WriteLine("  " + line.ToString().PadLeft(4) + "  "
                                  + "  " + line_pointer[line].ToString().PadLeft(8)
                                  + "  " + (line_pointer[line + 1] - 1).ToString().PadLeft(8) + "");
            }

            Console.WriteLine("");
            for (line = 0; line < line_num; line++)
            {
                string cout = "  " + line.ToString().PadLeft(4) + "  ";
                for (i = line_pointer[line]; i < line_pointer[line + 1]; i++)
                {
                    cout += "  " + line_data[i].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }

        public static void xyzl_data_read(string input_filename, int line_num, int line_data_num,
                ref int[] line_pointer, ref int[] line_data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZL_DATA_READ reads the data in an XYZL file.
            //
            //  Discussion:
            //
            //    This routine assumes that the file contains exactly three kinds of
            //    records:
            //
            //    COMMENTS which begin with a '#' character in column 1;
            //    BLANKS which contain nothing but 'whitespace';
            //    LINE ITEMS, which are indices of points on a line.
            //
            //    The routine ignores comments and blanks and returns
            //    the number of line items.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 October 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string INPUT_FILENAME, the name of the input file.
            //
            //    Input, int LINE_NUM, the number of lines.
            //
            //    Input, int LINE_DATA_NUM, the number of line items.
            //
            //    Output, int LINE_POINTER[LINE_NUM+1], pointers to the
            //    first line item for each line.
            //
            //    Output, int LINE_DATA[LINE_DATA_NUM], the line items.
            //
        {
            int ierror;
            int ilo;
            string[] input;
            int line;
            int n;
            string text;

            try
            {
                input = File.ReadAllLines(input_filename);
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("XYZL_DATA_READ - Fatal error!");
                Console.WriteLine("  Cannot open the input file \"" + input_filename + "\".");
                return;
            }

            line = 0;
            line_pointer[0] = 0;

            int index = 0;

            while (line < line_num)
            {
                try
                {
                    text = input[index];
                    index++;
                }
                catch (Exception e)
                {
                    Console.WriteLine("");
                    Console.WriteLine("XYZL_DATA_READ - Fatal error!");
                    Console.WriteLine("  Unexpected end of information.");
                    return;
                }

                if (text[0] == '#' || typeMethods.s_len_trim(text) == 0)
                {
                    continue;
                }

                n = typeMethods.s_word_count(text);
                line_pointer[line + 1] = line_pointer[line] + n;

                ilo = line_pointer[line];

                i4vec r = typeMethods.s_to_i4vec(text, n);

                if (r.error)
                {
                    Console.WriteLine("");
                    Console.WriteLine("XYZL_DATA_READ - Fatal error!");
                    Console.WriteLine("  Error from S_TO_I4VEC.");
                    return;
                }

                for (int ix = 0; ix < r.ivec.Length; ix++)
                {
                    line_data[ilo + ix] = r.ivec[ix];
                }

                line = line + 1;
            }

        }

        public static void xyzl_data_write(ref List<string> output_unit, int point_num, int line_num,
                int line_data_num, int[] line_pointer, int[] line_data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZL_DATA_WRITE writes the data of an XYZL file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, ofstream &OUTPUT_UNIT, a pointer to the XY file.
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int LINE_NUM, the number of lines.
            //
            //    Input, int LINE_DATA_NUM, the number of line items.
            //
            //    Input, int LINE_POINTER[LINE_NUM+1], pointers to the
            //    first line item for each line.
            //
            //    Input, int LINE_DATA[LINE_DATA_NUM], indices
            //    of points that form lines.
            //
        {
            int i;
            int line;

            for (line = 0; line < line_num; line++)
            {
                string cout = "";
                for (i = line_pointer[line]; i < line_pointer[line + 1]; i++)
                {
                    cout += "  " + line_data[i];
                }

                output_unit.Add(cout);
            }
        }

        public static void xyzl_example(int point_num, int line_num, int line_data_num,
                ref double[] xyz, ref int[] line_pointer, ref int[] line_data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZL_EXAMPLE sets data suitable for a pair of XYZ and XYZL files.
            //
            //  Discussion:
            //
            //    There are 8 points.
            //    There are 6 lines.
            //    There are 18 line items.
            //
            //       8------7
            //      /|     /|
            //     / |    / |
            //    5------6  |
            //    |  4---|--3
            //    | /    | /
            //    |/     |/
            //    1------2
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int LINE_NUM, the number of lines.
            //
            //    Input, int LINE_DATA_NUM, the number of line items.
            //
            //    Output, double XYZ[3*POINT_NUM], the point coordinates.
            //
            //    Output, int LINE_POINTER[LINE_NUM+1], pointers to the
            //    first line item for each line.
            //
            //    Output, int LINE_DATA[LINE_DATA_NUM], indices
            //    of points that form lines.
            //
        {
            int LINE_DATA_NUM = 18;
            int LINE_NUM = 6;
            int POINT_NUM = 8;

            int[] LINE_DATA =
            {
                0, 1, 2, 3, 0,
                4, 5, 6, 7, 4,
                0, 4,
                1, 5,
                2, 6,
                3, 7
            };
            int[] LINE_POINTER = {0, 5, 10, 12, 14, 16, 18};
            double[] XYZ =
            {
                0.0, 0.0, 0.0,
                1.0, 0.0, 0.0,
                1.0, 1.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0,
                1.0, 0.0, 1.0,
                1.0, 1.0, 1.0,
                0.0, 1.0, 1.0
            };

            typeMethods.i4vec_copy(line_data_num, LINE_DATA, ref line_data);
            typeMethods.i4vec_copy(line_num + 1, LINE_POINTER, ref line_pointer);
            typeMethods.r8vec_copy(3 * point_num, XYZ, ref xyz);

        }

        public static void xyzl_example_size(ref int point_num, ref int line_num, ref int line_data_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZL_EXAMPLE_SIZE sizes the data to be created by XYZL_EXAMPLE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, int *POINT_NUM, the number of points.
            //
            //    Output, int *LINE_NUM, the number of lines.
            //
            //    Output, int *LINE_DATA_NUM, the number of line items.
            //
        {
            line_data_num = 18;
            line_num = 6;
            point_num = 8;
        }

        public static void xyzl_header_print(int point_num, int line_num, int line_data_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZL_HEADER_PRINT prints the header of an XYZL file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    03 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int LINE_NUM, the number of lines.
            //
            //    Input, int LINE_DATA_NUM, the number of line items.
            //
        {
            Console.WriteLine("");
            Console.WriteLine("  Number of points     = " + point_num + "");
            Console.WriteLine("  Number of lines      = " + line_num + "");
            Console.WriteLine("  Number of line items = " + line_data_num + "");

            return;
        }

        public static void xyzl_header_read(string input_filename, ref int line_num,
                ref int line_data_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZL_HEADER_READ determines the number of line items in an XYZL file.
            //
            //  Discussion:
            //
            //    This routine assumes that the file contains exactly three kinds of
            //    records:
            //
            //    COMMENTS which begin with a '#' character in column 1;
            //    BLANKS which contain nothing but 'whitespace';
            //    LINE ITEMS, which are indices of points on a line.
            //
            //    The routine ignores comments and blanks and returns
            //    the number of line items.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 October 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string INPUT_FILENAME, the name of the input file.
            //
            //    Output, int *LINE_NUM, the number of lines.
            //
            //    Output, int *LINE_DATA_NUM, the number of line items.
            //
        {
            int i;
            int i4_val;
            int ierror;
            string[] input;
            int length;
            int n;
            string text;

            line_data_num = 0;
            line_num = 0;

            try
            {
                input = File.ReadAllLines(input_filename);
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("XYZL_HEADER_READ - Fatal error!");
                Console.WriteLine("  Cannot open the input file \"" + input_filename + "\".");
                return;
            }

            int index = 0;

            for (;;)
            {
                try
                {
                    text = input[index];
                    index++;
                }
                catch (Exception e)
                {
                    break;
                }

                if (text[0] == '#' || typeMethods.s_len_trim(text) == 0)
                {
                    continue;
                }

                n = typeMethods.s_word_count(text);

                line_data_num = line_data_num + n;

                line_num = line_num + 1;
            }

        }

        public static void xyzl_header_write(string output_filename, ref List<string> output_unit,
                int point_num, int line_num, int line_data_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZL_HEADER_WRITE writes the header of an XYZL file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string OUTPUT_FILENAME, the name of the file.
            //
            //    Input, ofstream &OUTPUT_UNIT, a pointer to the file to contain the data.
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int LINE_NUM, the number of lines.
            //
            //    Input, int LINE_DATA_NUM, the number of line items.
            //
        {
            output_unit.Add("#  " + output_filename + "");
            output_unit.Add("#  created by xyz_io::xyzl_header_write.C");
            output_unit.Add("#");
            output_unit.Add("#  Number of points     = " + point_num + "");
            output_unit.Add("#  Number of lines      = " + line_num + "");
            output_unit.Add("#  Number of line items = " + line_data_num + "");
            output_unit.Add("#");
        }

        public static void xyzl_write(string output_filename, int point_num, int line_num,
                int line_data_num, int[] line_pointer, int[] line_data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XYZL_WRITE writes the header and data for an XYZL file.
            // 
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            // 
            //    05 January 2009
            // 
            //  Author:
            // 
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string OUTPUT_FILENAME, the name of the file.
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int LINE_NUM, the number of lines.
            //
            //    Input, int LINE_DATA_NUM, the number of line items.
            //
            //    Input, int LINE_POINTER[LINE_NUM+1], pointers to the
            //    first line item for each line.
            //
            //    Input, int LINE_DATA[LINE_DATA_NUM], indices
            //    of points that form lines.
            //
        {
            List<string> output_unit = new List<string>();
            //
            //  Write the header.
            //
            xyzl_header_write(output_filename, ref output_unit, point_num, line_num,
                line_data_num);
            //
            //  Write the data.
            //
            xyzl_data_write(ref output_unit, point_num, line_num, line_data_num, line_pointer,
                line_data);

            //
            //  Open the output file.
            //
            try
            {
                File.WriteAllLines(output_filename, output_unit);
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("XYZL_WRITE - Fatal error!");
                Console.WriteLine("  Cannot open the output file \"" + output_filename + "\".");
            }

        }
    }
}