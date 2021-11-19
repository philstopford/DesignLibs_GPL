using System;
using Burkardt.IO;

namespace XYTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for XY_IO_TEST.
        //
        //  Discussion:
        //
        //    XY_IO_TEST tests the XY_IO library.
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
    {
        Console.WriteLine("");
        Console.WriteLine("XY_IO_TEST:");
            
        Console.WriteLine("  Test the XY_IO library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();

        Console.WriteLine("");
        Console.WriteLine("XY_IO_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");

    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests XY_EXAMPLE, XY_WRITE.
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
    {
        int POINT_NUM = 300;

        string file_name = "xy_io_prb_01.xy";
        int point_num = POINT_NUM;
        double[] xy = new double[2 * POINT_NUM];

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  XY_EXAMPLE sets up sample XY data.");
        Console.WriteLine("  XY_WRITE writes an XY file.");

        XY.xy_example(point_num, ref xy);

        Console.WriteLine("");
        Console.WriteLine("  XY_EXAMPLE has created the data.");

        XY.xy_write(file_name, point_num, xy);

        Console.WriteLine("");
        Console.WriteLine("  XY_WRITE wrote the header and data for \"" + file_name + "\"");
        Console.WriteLine("  Number of points = " + point_num + "");
        Console.WriteLine();
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests XY_READ.
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
    {
        string file_name = "xy_io_prb_02.xy";
        int i;
        int k;
        int point_num = 0;
        double[] xy = null;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  XY_READ reads the header and data of an XY file.");

        XY.xy_write_test(file_name);

        Console.WriteLine("");
        Console.WriteLine("  XY_WRITE_TEST created data and wrote it to \""
                          + file_name + "\"");
        //
        //  Now have XY_READ try to read it.
        //
        XY.xy_read(file_name, ref point_num, ref xy);

        Console.WriteLine("");
        Console.WriteLine("  XY_READ read the test data successfully.");
        Console.WriteLine("  Number of points = " + point_num + "");
        Console.WriteLine("");
        Console.WriteLine("  Sample data:");
        Console.WriteLine("");
        for (k = 0; k <= 9; k++)
        {
            i = ((9 - k) * 0 + k * (point_num - 1)) / 9;
            Console.WriteLine(i.ToString().PadLeft(4) + "  "
                                                      + xy[0 + i * 2].ToString().PadLeft(10) + "  "
                                                      + xy[1 + i * 2].ToString().PadLeft(10) + "");
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests XYL_EXAMPLE.
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
    {
        int[] line_data;
        int[] line_pointer;
        int line_data_num = 0;
        int line_num = 0;
        int point_num = 0;
        double[] xy;
        string xy_filename = "house.xy";
        string xyl_filename = "house.xyl";

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  XYL_EXAMPLE sets up XY and XYL data.");

        XY.xyl_example_size(ref point_num, ref line_num, ref line_data_num);

        Console.WriteLine("");
        Console.WriteLine("  Example has:");
        Console.WriteLine("");
        Console.WriteLine("  Number of points     = " + point_num + "");
        Console.WriteLine("  Number of lines      = " + line_num + "");
        Console.WriteLine("  Number of line items = " + line_data_num + "");

        line_data = new int[line_data_num];
        line_pointer = new int[line_num + 1];
        xy = new double[2 * point_num];

        XY.xyl_example(point_num, line_num, line_data_num, ref xy, ref line_pointer, ref line_data);

        XY.xy_write(xy_filename, point_num, xy);

        XY.xyl_write(xyl_filename, point_num, line_num, line_data_num,
            line_pointer, line_data);

        Console.WriteLine("");
        Console.WriteLine("  Wrote the XY file \"" + xy_filename + "\",");
        Console.WriteLine("  and the XYL file \"" + xyl_filename + "\".");

    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests XYL_READ.
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
    {
        int i;
        int j;
        int line;
        int[] line_data;
        int[] line_pointer;
        int line_data_num = 0;
        int line_num = 0;
        int point_num = 0;
        double[] xy;
        string xy_filename = "house.xy";
        string xyl_filename = "house.xyl";

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  XY_HEADER_READ  reads the header of an XY  file.");
        Console.WriteLine("  XY_DATA_READ    reads the data   of an XY  file.");
        Console.WriteLine("  XYL_HEADER_READ reads the header of an XYL file.");
        Console.WriteLine("  XYL_DATA_READ   reads the data   of an XYL file.");

        Console.WriteLine("");
        Console.WriteLine("  Examine XY file \"" + xy_filename + "\".");

        XY.xy_header_read(xy_filename, ref point_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of points     = " + point_num + "");

        xy = new double[2 * point_num];

        XY.xy_data_read(xy_filename, point_num, ref xy);

        Console.WriteLine("");
        Console.WriteLine("  Point data:");
        Console.WriteLine("");

        for (i = 0; i < point_num; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + xy[0 + 2 * i].ToString().PadLeft(10)
                                   + "  " + xy[1 + 2 * i].ToString().PadLeft(10) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Examine XYL file \"" + xyl_filename + "\".");

        XY.xyl_header_read(xyl_filename, ref line_num, ref line_data_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of lines      = " + line_num + "");
        Console.WriteLine("  Number of line items = " + line_data_num + "");

        line_data = new int[line_data_num];
        line_pointer = new int[line_num + 1];

        XY.xyl_data_read(xyl_filename, line_num, line_data_num, ref line_pointer, ref line_data);

        Console.WriteLine("");
        Console.WriteLine("  Line pointers:");
        Console.WriteLine("");

        for (line = 0; line < line_num; line++)
        {
            Console.WriteLine("  " + line.ToString().PadLeft(4)
                                   + "  " + line_pointer[line].ToString().PadLeft(8)
                                   + "  " + (line_pointer[line + 1] - 1).ToString().PadLeft(8) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Line data:");
        Console.WriteLine("");

        for (line = 0; line < line_num; line++)
        {
            string cout = "  " + line.ToString().PadLeft(4) + "  ";
            for (j = line_pointer[line]; j < line_pointer[line + 1]; j++)
            {
                cout += "  " + line_data[j].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests XYF_EXAMPLE.
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
    {
        int[] face_data;
        int[] face_pointer;
        int face_data_num = 0;
        int face_num = 0;
        int point_num = 0;
        double[] xy;
        string xy_filename = "annulus.xy";
        string xyf_filename = "annulus.xyf";

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  XYF_EXAMPLE sets up XY and XYF data.");

        XY.xyf_example_size(ref point_num, ref face_num, ref face_data_num);

        Console.WriteLine("");
        Console.WriteLine("  Example has:");
        Console.WriteLine("");
        Console.WriteLine("  Number of points     = " + point_num + "");
        Console.WriteLine("  Number of faces      = " + face_num + "");
        Console.WriteLine("  Number of face items = " + face_data_num + "");

        face_data = new int[face_data_num];
        face_pointer = new int[face_num + 1];
        xy = new double[2 * point_num];

        XY.xyf_example(point_num, face_num, face_data_num, ref xy, ref face_pointer, ref face_data);

        XY.xy_write(xy_filename, point_num, xy);

        XY.xyf_write(xyf_filename, point_num, face_num, face_data_num,
            face_pointer, face_data);

        Console.WriteLine("");
        Console.WriteLine("  Wrote the XY file \"" + xy_filename + "\",");
        Console.WriteLine("  and the XYF file \"" + xyf_filename + "\".");
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests XYF_READ.
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
    {
        int face;
        int[] face_data;
        int[] face_pointer;
        int face_data_num = 0;
        int face_num = 0;
        int i;
        int j;
        int point_num = 0;
        double[] xy;
        string xy_filename = "annulus.xy";
        string xyf_filename = "annulus.xyf";

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  XY_HEADER_READ  reads the header of an XY  file.");
        Console.WriteLine("  XY_DATA_READ    reads the data   of an XY  file.");
        Console.WriteLine("  XYF_HEADER_READ reads the header of an XYF file.");
        Console.WriteLine("  XYF_DATA_READ   reads the data   of an XYF file.");

        Console.WriteLine("");
        Console.WriteLine("  Examine XY file \"" + xy_filename + "\".");

        XY.xy_header_read(xy_filename, ref point_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of points     = " + point_num + "");

        xy = new double[2 * point_num];

        XY.xy_data_read(xy_filename, point_num, ref xy);

        Console.WriteLine("");
        Console.WriteLine("  Point data:");
        Console.WriteLine("");

        for (i = 0; i < point_num; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + xy[0 + 2 * i].ToString().PadLeft(10)
                                   + "  " + xy[1 + 2 * i].ToString().PadLeft(10) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Examine XYF file \"" + xyf_filename + "\".");

        XY.xyf_header_read(xyf_filename, ref face_num, ref face_data_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of faces      = " + face_num + "");
        Console.WriteLine("  Number of face items = " + face_data_num + "");

        face_data = new int[face_data_num];
        face_pointer = new int[face_num + 1];

        XY.xyf_data_read(xyf_filename, face_num, face_data_num, ref face_pointer, ref face_data);

        Console.WriteLine("");
        Console.WriteLine("  Face pointers:");
        Console.WriteLine("");

        for (face = 0; face < face_num; face++)
        {
            Console.WriteLine("  " + face.ToString().PadLeft(4)
                                   + "  " + face_pointer[face].ToString().PadLeft(8)
                                   + "  " + (face_pointer[face + 1] - 1).ToString().PadLeft(8) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Face data:");
        Console.WriteLine("");

        for (face = 0; face < face_num; face++)
        {
            string cout = "  " + face.ToString().PadLeft(4) + "  ";
            for (j = face_pointer[face]; j < face_pointer[face + 1]; j++)
            {
                cout += "  " + face_data[j].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
        }
    }
}