using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using Burkardt.Types;

namespace Burkardt.IO;

public static class XY
{
    public static void xy_data_print(int point_num, double[] xy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XY_DATA_PRINT prints the data for an XY file.
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
        //    Input, double XY[2*POINT_NUM], the arrays of coordinate data.
        //
    {
        int j;

        for (j = 0; j < point_num; j++)
        {
            Console.WriteLine(xy[0 + 2 * j].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                                                                                   + xy[1 + 2 * j].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

        }
    }

    public static void xy_data_read(string input_filename, int point_num, ref double[] xy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XY_DATA_READ reads the data in an XY file.
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
        //    Input, int POINT_NUM, the number of points.
        //
        //    Output, double XY[2*POINT_NUM], the point coordinates.
        //
    {
        string[] input;

        try
        {
            input = File.ReadAllLines(input_filename);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("XY_DATA_READ - Fatal error!");
            Console.WriteLine("  Cannot open the input file \"" + input_filename + "\".");
            return;
        }

        int j = 0;
        int index = 0;

        while (j < point_num)
        {
            string text;
            try
            {
                text = input[index];
                index++;
            }
            catch (Exception)
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
            r8vec r = typeMethods.s_to_r8vec(text, 2);

            switch (r.error)
            {
                case true:
                    Console.WriteLine("");
                    Console.WriteLine("XY_DATA_READ - Fatal error!");
                    Console.WriteLine("  S_TO_R8VEC returned error flag.");
                    return;
            }

            xy[0 + j * 2] = r.rvec[0];
            xy[1 + j * 2] = r.rvec[1];
            j += 1;
        }

    }

    public static void xy_data_write(ref List<string> output_unit, int point_num, double[] xy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XY_DATA_WRITE writes the data for an XY file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 December 2008
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
        //    Input, double XY[2*POINT_NUM], the arrays of coordinate data.
        //
    {
        int j;

        for (j = 0; j < point_num; j++)
        {
            output_unit.Add(xy[0 + 2 * j].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                                                                                 + xy[1 + 2 * j].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

        }
    }

    public static void xy_example(int point_num, ref double[] xy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XY_EXAMPLE sets up sample XY data suitable for an XY file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 December 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int POINT_NUM, the number of points.
        //
        //    Output, double XY[2*POINT_NUM], the arrays of coordinate data.
        //
    {
        int j;

        const int turns = 5;

        for (j = 0; j < point_num; j++)
        {
            double r = j / (double) (point_num - 1);
            double theta = turns * r * (2.0 * Math.PI);
            xy[0 + j * 2] = r * Math.Cos(theta);
            xy[1 + j * 2] = r * Math.Sin(theta);
        }

        for (j = 0; j < point_num; j++)
        {
            xy[0 + j * 2] = 0.5 * (1.0 + xy[0 + j * 2]);
            xy[1 + j * 2] = 0.5 * (1.0 + xy[1 + j * 2]);
        }

    }

    public static void xy_header_print(int point_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XY_HEADER_PRINT prints the header of an XY file.
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
        Console.WriteLine("#");
        Console.WriteLine("#  Number of points = " + point_num + "");
    }

    public static void xy_header_read(string input_filename, ref int point_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XY_HEADER_READ reads the header of an XY file.
        //
        //  Discussion:
        //
        //    All we do here is count the number of records that are not comments 
        //    and not blank.  Each such record is assumed to represent a single point 
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
        //    Input, string INPUT_FILENAME, the name of the input file.
        //
        //    Output, int *POINT_NUM, the number of points.
        //
    {
        string[] input;

        point_num = 0;

        try
        {
            input = File.ReadAllLines(input_filename);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("XY_HEADER_READ - Fatal error!");
            Console.WriteLine("  Cannot open the input file \"" + input_filename + "\".");
            return;
        }

        int index = 0;
        while (true)
        {
            string text;
            try
            {
                text = input[index];
                index++;
            }
            catch (Exception)
            {
                break;
            }

            if (text[0] == '#' || typeMethods.s_len_trim(text) == 0)
            {
                continue;
            }

            point_num += 1;
        }

    }

    public static void xy_header_write(string output_filename, ref List<string> output_unit,
            int point_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XY_HEADER_WRITE writes the header of an XY file.
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
        output_unit.Add("#  created by xy_io::xy_header_write.C");
        output_unit.Add("#");
        output_unit.Add("#  Number of points = " + point_num + "");
        output_unit.Add("#");
    }

    public static void xy_read(string input_filename, ref int point_num, ref double[] xy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XY_READ reads the header and data from an XY file.
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
        //    Input, string INPUT_FILENAME, the name of the XY file.
        //
        //    Output, int *POINT_NUM, the number of points.
        //
        //    Output, double *XY[2*(*POINT_NUM)], the point coordinates.
        //
    {
        xy_header_read(input_filename, ref point_num);

        xy = new double[2 * point_num];

        xy_data_read(input_filename, point_num, ref xy);

    }

    public static void xy_read_test(string input_filename)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XY_READ_TEST tests the XY file read routines.
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
        //    Input, string INPUT_FILENAME, the name of the XY file.
        //
    {
        int point_num = 0;
        double[] xy = null;

        xy_read(input_filename, ref point_num, ref xy);

    }

    public static void xy_write(string output_filename, int point_num, double[] xy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XY_WRITE writes the header and data for an XY file.
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
        //    Input, double XY[2*POINT_NUM], the arrays of coordinate data.
        //
    {
        List<string> output_unit = new();

        //
        //  Write the header.
        //
        xy_header_write(output_filename, ref output_unit, point_num);
        //
        //  Write the data.
        //
        xy_data_write(ref output_unit, point_num, xy);

        //
        //  Open the output file.
        //
        try
        {
            File.WriteAllLines(output_filename, output_unit);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("XY_WRITE - Fatal error!");
            Console.WriteLine("  Cannot open the output file \"" + output_filename + "\".");
        }
    }

    public static void xy_write_test(string output_filename)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XY_WRITE_TEST tests the XY write routines.
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
        const int point_num = 100;
        //
        //  Allocate memory.
        //
        double[] xy = new double[2 * point_num];
        //
        //  Set the data.
        //
        xy_example(point_num, ref xy);
        //
        //  Write the data to the file.
        //
        xy_write(output_filename, point_num, xy);

    }

    public static void xyf_data_print(int point_num, int face_num,
            int face_data_num, int[] face_pointer, int[] face_data)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYF_DATA_PRINT prints the data of an XYF file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 January 2009
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
        int face;

        for (face = 0; face < face_num; face++)
        {
            Console.WriteLine("  " + face.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + face_pointer[face].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + (face_pointer[face + 1] - 1).ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }

        Console.WriteLine("");

        for (face = 0; face < face_num; face++)
        {
            string cout = "";
            int i;
            for (i = face_pointer[face]; i < face_pointer[face + 1]; i++)
            {
                cout += "  " + face_data[i];
            }

            Console.WriteLine(cout);
        }
    }

    public static void xyf_data_read(string input_filename, int face_num, int face_data_num,
            ref int[] face_pointer, ref int[] face_data)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYF_DATA_READ reads the data in an XYF file.
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
        //    The routine ignores comments and blank faces and returns
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
        string[] input;

        try
        {
            input = File.ReadAllLines(input_filename);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("XYF_DATA_READ - Fatal error!");
            Console.WriteLine("  Cannot open the input file \"" + input_filename + "\".");
            return;
        }

        int face = 0;
        face_pointer[0] = 0;
        int index = 0;

        while (face < face_num)
        {
            string text;
            try
            {
                text = input[index];
                index++;

            }
            catch (Exception)
            {
                Console.WriteLine("");
                Console.WriteLine("XYF_DATA_READ - Fatal error!");
                Console.WriteLine("  Unexpected end of information.");
                return;
            }

            if (text[0] == '#' || typeMethods.s_len_trim(text) == 0)
            {
                continue;
            }

            int n = typeMethods.s_word_count(text);
            face_pointer[face + 1] = face_pointer[face] + n;

            int ilo = face_pointer[face];

            i4vec r = typeMethods.s_to_i4vec(text, n);

            switch (r.error)
            {
                case true:
                    Console.WriteLine("");
                    Console.WriteLine("XYF_DATA_READ - Fatal error!");
                    Console.WriteLine("  Error from S_TO_I4VEC.");
                    return;
            }

            for (int ix = 0; ix < r.ivec.Length; ix++)
            {
                face_data[ilo + ix] = r.ivec[ix];
            }

            face += 1;
        }
    }

    public static void xyf_data_write(ref List<string> output_unit, int point_num, int face_num,
            int face_data_num, int[] face_pointer, int[] face_data)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYF_DATA_WRITE writes the data of an XYF file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 January 2009
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
        int face;

        for (face = 0; face < face_num; face++)
        {
            int i;
            for (i = face_pointer[face]; i < face_pointer[face + 1]; i++)
            {
                output_unit.Add("  " + face_data[i]);
            }

            output_unit.Add("");
        }
    }

    public static void xyf_example(int point_num, int face_num, int face_data_num, ref double[] xy,
            ref int[] face_pointer, ref int[] face_data)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYF_EXAMPLE sets data suitable for a pair of XY and XYF files.
        //
        //  Discussion:
        //
        //    There are 65 points.
        //    There are 48 faces.
        //    There are 48*4=192 face items.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 January 2009
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
        //    Output, double XY[2*POINT_NUM], the point coordinates.
        //
        //    Output, int FACE_POINTER[FACE_NUM+1], pointers to the
        //    first face item for each face.
        //
        //    Output, int FACE_DATA[FACE_DATA_NUM], indices
        //    of points that form faces.
        //
    {
        int i;
        int j;
        const int n_t = 13;
        const int n_r = 5;
        const double r_min = 1.0;
        const double r_max = 3.0;
        const double t_min = 3.141592653589793;
        const double t_max = 0.0;

        int k = 0;
        for (j = 1; j <= n_r; j++)
        {
            double r = ((n_r - j) * r_min
                        + (j - 1) * r_max)
                       / (n_r - 1);

            for (i = 1; i <= n_t; i++)
            {
                double t = ((n_t - i) * t_min
                            + (i - 1) * t_max)
                           / (n_t - 1);

                xy[0 + k * 2] = r * Math.Cos(t);
                xy[1 + k * 2] = r * Math.Sin(t);
                k += 1;
            }
        }

        int face = 0;
        k = 0;
        face_pointer[face] = k;

        for (j = 1; j < n_r; j++)
        {
            for (i = 1; i < n_t; i++)
            {
                face += 1;

                face_data[k] = i + (j - 1) * n_t - 1;
                k += 1;
                face_data[k] = i + 1 + (j - 1) * n_t - 1;
                k += 1;
                face_data[k] = i + 1 + j * n_t - 1;
                k += 1;
                face_data[k] = i + j * n_t - 1;
                k += 1;
                face_pointer[face] = k;
            }
        }

    }

    public static void xyf_example_size(ref int point_num, ref int face_num, ref int face_data_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYF_EXAMPLE_SIZE sizes the data to be created by XYF_EXAMPLE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 January 2009
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
        const int n_t = 13;
        const int n_r = 5;

        face_data_num = 4 * (n_t - 1) * (n_r - 1);
        face_num = (n_t - 1) * (n_r - 1);
        point_num = n_t * n_r;

    }

    public static void xyf_header_print(int point_num, int face_num, int face_data_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYF_HEADER_PRINT prints the header of an XYF file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 January 2009
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
    }

    public static void xyf_header_read(string input_filename, ref int face_num,
            ref int face_data_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYF_HEADER_READ determines the number of face items in an XYF file.
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
        string[] input;

        face_data_num = 0;
        face_num = 0;

        try
        {
            input = File.ReadAllLines(input_filename);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("XYF_HEADER_READ - Fatal error!");
            Console.WriteLine("  Cannot open the input file \"" + input_filename + "\".");
            return;
        }

        int index = 0;

        for (;;)
        {
            string text;
            try
            {
                text = input[index];
                index++;
            }
            catch (Exception)
            {
                break;
            }

            if (text[0] == '#' || typeMethods.s_len_trim(text) == 0)
            {
                continue;
            }

            int n = typeMethods.s_word_count(text);

            face_data_num += n;

            face_num += 1;
        }
    }

    public static void xyf_header_write(string output_filename, ref List<string> output_unit,
            int point_num, int face_num, int face_data_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYF_HEADER_WRITE writes the header of an XYF file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 January 2009
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
        output_unit.Add("#  created by xy_io::xyf_header_write.C");
        output_unit.Add("#");
        output_unit.Add("#  Number of points     = " + point_num + "");
        output_unit.Add("#  Number of faces      = " + face_num + "");
        output_unit.Add("#  Number of face items = " + face_data_num + "");
        output_unit.Add("#");

    }

    public static void xyf_write(string output_filename, int point_num, int face_num,
            int face_data_num, int[] face_pointer, int[] face_data)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYF_WRITE writes the header and data for an XYF file.
        // 
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        // 
        //    06 January 2009
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
        List<string> output_unit = new();

        //
        //  Write the header.
        //
        xyf_header_write(output_filename, ref output_unit, point_num, face_num,
            face_data_num);
        //
        //  Write the data.
        //
        xyf_data_write(ref output_unit, point_num, face_num, face_data_num, face_pointer,
            face_data);

        try
        {
            File.WriteAllLines(output_filename, output_unit);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("XYF_WRITE - Fatal error!");
            Console.WriteLine("  Cannot open the output file \"" + output_filename + "\".");
        }
    }

    public static void xyl_data_print(int point_num, int line_num,
            int line_data_num, int[] line_pointer, int[] line_data)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYL_DATA_PRINT prints the data of an XYL file.
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
        int line;

        for (line = 0; line < line_num; line++)
        {
            Console.WriteLine("  " + line.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + line_pointer[line].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + (line_pointer[line + 1] - 1).ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }

        Console.WriteLine("");

        for (line = 0; line < line_num; line++)
        {
            string cout = "";
            int i;
            for (i = line_pointer[line]; i < line_pointer[line + 1]; i++)
            {
                cout += "  " + line_data[i];
            }

            Console.WriteLine(cout);
        }
    }

    public static void xyl_data_read(string input_filename, int line_num, int line_data_num,
            ref int[] line_pointer, ref int[] line_data)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYL_DATA_READ reads the data in an XYL file.
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
        string[] input;

        try
        {
            input = File.ReadAllLines(input_filename);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("XYL_DATA_READ - Fatal error!");
            Console.WriteLine("  Cannot open the input file \"" + input_filename + "\".");
            return;
        }

        int line = 0;
        line_pointer[0] = 0;

        int index = 0;

        while (line < line_num)
        {
            string text;
            try
            {
                text = input[index];
                index++;
            }
            catch (Exception)
            {
                Console.WriteLine("");
                Console.WriteLine("XYL_DATA_READ - Fatal error!");
                Console.WriteLine("  Unexpected end of information.");
                return;
            }

            if (text[0] == '#' || typeMethods.s_len_trim(text) == 0)
            {
                continue;
            }

            int n = typeMethods.s_word_count(text);
            line_pointer[line + 1] = line_pointer[line] + n;

            int ilo = line_pointer[line];

            i4vec r = typeMethods.s_to_i4vec(text, n);

            switch (r.error)
            {
                case true:
                    Console.WriteLine("");
                    Console.WriteLine("XYL_DATA_READ - Fatal error!");
                    Console.WriteLine("  Error from S_TO_I4VEC.");
                    return;
            }

            for (int ix = 0; ix < r.ivec.Length; ix++)
            {
                line_data[ilo + ix] = r.ivec[ix];
            }

            line += 1;
        }

    }

    public static void xyl_data_write(ref List<string> output_unit, int point_num, int line_num,
            int line_data_num, int[] line_pointer, int[] line_data)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYL_DATA_WRITE writes the data of an XYL file.
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
        int line;

        for (line = 0; line < line_num; line++)
        {
            string cout = "";
            int i;
            for (i = line_pointer[line]; i < line_pointer[line + 1]; i++)
            {
                cout += "  " + line_data[i];
            }

            output_unit.Add(cout);
        }
    }

    public static void xyl_example(int point_num, int line_num, int line_data_num, ref double[] xy,
            ref int[] line_pointer, ref int[] line_data)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYL_EXAMPLE sets data suitable for a pair of XY and XYL files.
        //
        //  Discussion:
        //
        //    There are 13 points.
        //    There are 3 lines.
        //    There are 15 line data items.
        //
        //         4 12-11
        //         .. | |
        //        .  .| |
        //       .   13 |
        //      .      .10
        //     .        .
        //    5          3
        //    |          |
        //    |     9--8 |
        //    |     |  | |
        //    |     |  | |
        //    |     6--7 |
        //    |          |
        //    1----------2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 December 2008
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
        //    Output, double XY[2*POINT_NUM], the point coordinates.
        //
        //    Output, int LINE_POINTER[LINE_NUM+1], pointers to the
        //    first line item for each line.
        //
        //    Output, int LINE_DATA[LINE_DATA_NUM], indices
        //    of points that form lines.
        //
    {
        int[] LINE_DATA =
        {
            0, 1, 2, 3, 4, 0,
            5, 6, 7, 8, 5,
            9, 10, 11, 12
        };
        int[] LINE_POINTER = {0, 6, 11, 15};
        double[] XY =
        {
            0.0, 0.0,
            6.0, 0.0,
            6.0, 7.0,
            3.0, 10.0,
            0.0, 7.0,
            4.0, 1.0,
            5.0, 1.0,
            5.0, 4.0,
            4.0, 4.0,
            5.0, 8.0,
            5.0, 11.0,
            4.0, 11.0,
            4.0, 9.0
        };

        typeMethods.i4vec_copy(line_data_num, LINE_DATA, ref line_data);
        typeMethods.i4vec_copy(line_num + 1, LINE_POINTER, ref line_pointer);
        typeMethods.r8vec_copy(2 * point_num, XY, ref xy);

    }

    public static void xyl_example_size(ref int point_num, ref int line_num, ref int line_data_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYL_EXAMPLE_SIZE sizes the data to be created by XYL_EXAMPLE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 December 2008
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
        line_data_num = 15;
        line_num = 3;
        point_num = 13;
    }

    public static void xyl_header_print(int point_num, int line_num, int line_data_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYL_HEADER_PRINT prints the header of an XYL file.
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

    }

    public static void xyl_header_read(string input_filename, ref int line_num,
            ref int line_data_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYL_HEADER_READ determines the number of line items in an XYL file.
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
        //    08 January 2009
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
        string[] input;

        line_data_num = 0;
        line_num = 0;

        try
        {
            input = File.ReadAllLines(input_filename);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("XYL_HEADER_READ - Fatal error!");
            Console.WriteLine("  Cannot open the input file \"" + input_filename + "\".");
            return;
        }

        int index = 0;

        for (;;)
        {
            string text;
            try
            {
                text = input[index];
                index++;
            }
            catch (Exception)
            {
                break;
            }

            if (text[0] == '#' || typeMethods.s_len_trim(text) == 0)
            {
                continue;
            }

            int n = typeMethods.s_word_count(text);

            line_data_num += n;

            line_num += 1;
        }
    }

    public static void xyl_header_write(string output_filename, ref List<string> output_unit,
            int point_num, int line_num, int line_data_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYL_HEADER_WRITE writes the header of an XYL file.
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
        output_unit.Add("#  created by xy_io::xyl_header_write.C");
        output_unit.Add("#");
        output_unit.Add("#  Number of points     = " + point_num + "");
        output_unit.Add("#  Number of lines      = " + line_num + "");
        output_unit.Add("#  Number of line items = " + line_data_num + "");
        output_unit.Add("#");

    }

    public static void xyl_write(string output_filename, int point_num, int line_num,
            int line_data_num, int[] line_pointer, int[] line_data)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYL_WRITE writes the header and data for an XYL file.
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
        List<string> output_unit = new();

        //
        //  Write the header.
        //
        xyl_header_write(output_filename, ref output_unit, point_num, line_num,
            line_data_num);
        //
        //  Write the data.
        //
        xyl_data_write(ref output_unit, point_num, line_num, line_data_num, line_pointer,
            line_data);

        //
        //  Open the output file.
        //
        try
        {
            File.WriteAllLines(output_filename, output_unit);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("XYL_WRITE - Fatal error!");
            Console.WriteLine("  Cannot open the output file \"" + output_filename + "\".");
        }
    }
}