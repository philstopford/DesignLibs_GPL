using System;
using System.Globalization;
using Burkardt.IO;

namespace XYZTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for XYZ_IO_TEST.
        //
        //  Discussion:
        //
        //    XYZ_IO_TEST tests the XYZ_IO library.
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
        Console.WriteLine("XYZ_IO_TEST:");
        Console.WriteLine("  Test the XYZ_IO library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();
        Console.WriteLine("");
        Console.WriteLine("XYZ_IO_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");

    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //   TEST01 tests XYZ_EXAMPLE, XYZ_WRITE.
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
        const string file_name = "helix.xyz";

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  XYZ_EXAMPLE sets up sample XYZ data.");
        Console.WriteLine("  XYZ_WRITE writes an XYZ file.");

        int point_num = XYZ.xyz_example_size();

        Console.WriteLine("  Example dataset size is " + point_num + "");

        double[] xyz = new double[3 * point_num];

        XYZ.xyz_example(point_num, ref xyz);

        XYZ.xyz_write(file_name, point_num, xyz);

        Console.WriteLine("");
        Console.WriteLine("  XYZ_WRITE wrote the header and data for \"" + file_name + "\".");
        Console.WriteLine("  Number of points = " + point_num + "");

    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests XYZ_HEADER_READ, XYZ_DATA_READ.
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
        const string file_name = "xyz_io_prb_02.xyz";
        int k;
        int point_num = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  XYZ_HEADER_READ reads the header of an XYZ file.");
        Console.WriteLine("  XYZ_DATA_READ reads the data of an XYZ file.");

        XYZ.xyz_write_test(file_name);

        Console.WriteLine("");
        Console.WriteLine("  XYZ_WRITE_TEST created some data.");

        XYZ.xyz_header_read(file_name, ref point_num);

        Console.WriteLine("");
        Console.WriteLine("  XYZ_HEADER_READ has read the header.");
        Console.WriteLine("");
        Console.WriteLine("  Number of points = " + point_num + "");

        double[] xyz = new double[3 * point_num];

        XYZ.xyz_data_read(file_name, point_num, ref xyz);

        Console.WriteLine("");
        Console.WriteLine("  XYZ_DATA_READ has read the data.");

        Console.WriteLine("");
        Console.WriteLine("  Sample data:");
        Console.WriteLine("");

        for (k = 1; k <= 11; k++)
        {
            int i = ((11 - k) * 1
                     + (k - 1) * point_num)
                / (11 - 1) - 1;
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + xyz[0 + i * 3].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + xyz[1 + i * 3].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + xyz[2 + i * 3].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests XYZL_EXAMPLE.
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
        int line_data_num = 0;
        int line_num = 0;
        int point_num = 0;
        const string xyz_filename = "cube.xyz";
        const string xyzl_filename = "cube.xyzl";

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  XYZL_EXAMPLE sets up XYZ and XYZL data.");

        XYZ.xyzl_example_size(ref point_num, ref line_num, ref line_data_num);

        Console.WriteLine("");
        Console.WriteLine("  Example has:");
        Console.WriteLine("");
        Console.WriteLine("  Number of points     = " + point_num + "");
        Console.WriteLine("  Number of lines      = " + line_num + "");
        Console.WriteLine("  Number of line items = " + line_data_num + "");

        int[] line_data = new int[line_data_num];
        int[] line_pointer = new int[line_num + 1];
        double[] xyz = new double[3 * point_num];

        XYZ.xyzl_example(point_num, line_num, line_data_num, ref xyz, ref line_pointer,
            ref line_data);

        XYZ.xyz_write(xyz_filename, point_num, xyz);

        XYZ.xyzl_write(xyzl_filename, point_num, line_num, line_data_num,
            line_pointer, line_data);

        Console.WriteLine("");
        Console.WriteLine("  Wrote the XYZ file \"" + xyz_filename + "\".");
        Console.WriteLine("  and the XYZL file \"" + xyzl_filename + "\".");

    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests XYZL_READ.
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
        int line;
        int line_data_num = 0;
        int line_num = 0;
        int point_num = 0;
        const string xyz_filename = "cube.xyz";
        const string xyzl_filename = "cube.xyzl";

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  XYZ_HEADER_READ  reads the header of an XYZ  file.");
        Console.WriteLine("  XYZ_DATA_READ    reads the data   of an XYZ  file.");
        Console.WriteLine("  XYZL_HEADER_READ reads the header of an XYZL file.");
        Console.WriteLine("  XYZL_DATA_READ   reads the data   of an XYZL file.");

        Console.WriteLine("");
        Console.WriteLine("  Examine XYZ file \"" + xyz_filename + "\".");

        XYZ.xyz_header_read(xyz_filename, ref point_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of points = " + point_num + "");

        double[] xyz = new double[3 * point_num];

        XYZ.xyz_data_read(xyz_filename, point_num, ref xyz);

        Console.WriteLine("");
        Console.WriteLine("  Point data:");
        Console.WriteLine("");

        for (i = 0; i < point_num; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + xyz[0 + i * 3].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + xyz[1 + i * 3].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + xyz[2 + i * 3].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Examine XYZL file \"" + xyzl_filename + "\".");

        XYZ.xyzl_header_read(xyzl_filename, ref line_num, ref line_data_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of lines      = " + line_num + "");
        Console.WriteLine("  Number of line items = " + line_data_num + "");

        int[] line_data = new int[line_data_num];
        int[] line_pointer = new int[line_num + 1];

        XYZ.xyzl_data_read(xyzl_filename, line_num, line_data_num, ref line_pointer,
            ref line_data);

        Console.WriteLine("");
        Console.WriteLine("  Line pointers:");
        Console.WriteLine("");

        for (line = 0; line < line_num; line++)
        {
            Console.WriteLine("  " + line.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + line_pointer[line].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + (line_pointer[line + 1] - 1).ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Line data:");
        Console.WriteLine("");

        for (line = 0; line < line_num; line++)
        {
            string cout = "  " + line.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "    ";
            int j;
            for (j = line_pointer[line]; j <= line_pointer[line + 1] - 1; j++)
            {
                cout += "  " + line_data[j].ToString(CultureInfo.InvariantCulture).PadLeft(8);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests XYZF_EXAMPLE.
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
    {
        int face_data_num = 0;
        int face_num = 0;
        int point_num = 0;
        const string xyz_filename = "cube.xyz";
        const string xyzf_filename = "cube.xyzf";

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  XYZF_EXAMPLE sets up XYZ and XYZF data.");

        XYZ.xyzf_example_size(ref point_num, ref face_num, ref face_data_num);

        Console.WriteLine("");
        Console.WriteLine("  Example has:");
        Console.WriteLine("");
        Console.WriteLine("  Number of points     = " + point_num + "");
        Console.WriteLine("  Number of faces      = " + face_num + "");
        Console.WriteLine("  Number of face items = " + face_data_num + "");

        int[] face_data = new int[face_data_num];
        int[] face_pointer = new int[face_num + 1];
        double[] xyz = new double[3 * point_num];

        XYZ.xyzf_example(point_num, face_num, face_data_num, ref xyz, ref face_pointer,
            ref face_data);

        XYZ.xyz_write(xyz_filename, point_num, xyz);

        XYZ.xyzf_write(xyzf_filename, point_num, face_num, face_data_num,
            face_pointer, face_data);

        Console.WriteLine("");
        Console.WriteLine("  Wrote the XYZ file \"" + xyz_filename + "\".");
        Console.WriteLine("  and the XYZF file \"" + xyzf_filename + "\".");

    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests XYZF_READ.
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
    {
        int face;
        int face_data_num = 0;
        int face_num = 0;
        int i;
        int point_num = 0;
        const string xyz_filename = "cube.xyz";
        const string xyzf_filename = "cube.xyzf";

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  XYZ_HEADER_READ  reads the header of an XYZ  file.");
        Console.WriteLine("  XYZ_DATA_READ    reads the data   of an XYZ  file.");
        Console.WriteLine("  XYZF_HEADER_READ reads the header of an XYZF file.");
        Console.WriteLine("  XYZF_DATA_READ   reads the data   of an XYZF file.");

        Console.WriteLine("");
        Console.WriteLine("  Examine XYZ file \"" + xyz_filename + "\".");

        XYZ.xyz_header_read(xyz_filename, ref point_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of points = " + point_num + "");

        double[] xyz = new double[3 * point_num];

        XYZ.xyz_data_read(xyz_filename, point_num, ref xyz);

        Console.WriteLine("");
        Console.WriteLine("  Point data:");
        Console.WriteLine("");

        for (i = 0; i < point_num; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + xyz[0 + i * 3].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + xyz[1 + i * 3].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + xyz[2 + i * 3].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Examine XYZF file \"" + xyzf_filename + "\".");

        XYZ.xyzf_header_read(xyzf_filename, ref face_num, ref face_data_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of faces      = " + face_num + "");
        Console.WriteLine("  Number of face items = " + face_data_num + "");

        int[] face_data = new int[face_data_num];
        int[] face_pointer = new int[face_num + 1];

        XYZ.xyzf_data_read(xyzf_filename, face_num, face_data_num, ref face_pointer,
            ref face_data);

        Console.WriteLine("");
        Console.WriteLine("  Face pointers:");
        Console.WriteLine("");

        for (face = 0; face < face_num; face++)
        {
            Console.WriteLine("  " + face.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + face_pointer[face].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + (face_pointer[face + 1] - 1).ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Face data:");
        Console.WriteLine("");

        for (face = 0; face < face_num; face++)
        {
            string cout = "  " + face.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "    ";
            int j;
            for (j = face_pointer[face]; j <= face_pointer[face + 1] - 1; j++)
            {
                cout += "  " + face_data[j].ToString(CultureInfo.InvariantCulture).PadLeft(8);
            }

            Console.WriteLine(cout);
        }
    }
}