using System;
using Burkardt.SphereNS;
using Burkardt.Types;

namespace SphereLatLongTriGridTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SPHERE_LLT_GRID_TEST.
        //
        //  Discussion:
        //
        //    SPHERE_LLT_GRID_TEST tests the SPHERE_LLT_GRID library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SPHERE_LLT_GRID_TEST");
        Console.WriteLine("  Test the SPHERE_LLT_GRID library.");

        sphere_llt_grid_point_count_test();
        sphere_llt_grid_points_test();
        sphere_llt_grid_line_count_test();
        sphere_llt_grid_lines_test();
        sphere_llt_grid_display_test();

        Console.WriteLine("");
        Console.WriteLine("SPHERE_LLT_GRID_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void sphere_llt_grid_point_count_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_LLT_GRID_POINT_COUNT_TEST tests SPHERE_LLT_GRID_POINT_COUNT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int lat_num;
        int long_log;
        int long_num;
        int point_num;

        Console.WriteLine("");
        Console.WriteLine("SPHERE_LLT_GRID_POINT_COUNT_TEST");
        Console.WriteLine("  SPHERE_LLT_GRID_POINT_COUNT counts the points used for a");
        Console.WriteLine("  grid based on triangles defined by latitude and longitude");
        Console.WriteLine("  lines on a sphere in 3D.");
        Console.WriteLine("");
        Console.WriteLine("     LAT_NUM    LONG_NUM   POINT_NUM");

        for (lat_num = 1; lat_num <= 17; lat_num += 2)
        {
            Console.WriteLine("");
            long_num = 1;
            for (long_log = 1; long_log <= 4; long_log++)
            {
                long_num *= 2;
                point_num = Grid_LatLong.sphere_llt_grid_point_count(lat_num, long_num);
                Console.WriteLine("  " + lat_num.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + long_num.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + point_num.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
            }
        }
    }

    private static void sphere_llt_grid_points_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_LLT_GRID_POINTS_TEST tests SPHERE_LLT_GRID_POINTS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int k;
        int lat_num;
        int long_num;
        int node_num;
        double[] node_xyz;
        double[] pc = { 0.0, 0.0, 0.0 };
        double r;

        lat_num = 3;
        long_num = 4;

        r = 10.0;

        Console.WriteLine("");
        Console.WriteLine("SPHERE_LLT_GRID_POINTS_TEST");
        Console.WriteLine("  SPHERE_LLT_POINTS produces latitude/longitude");
        Console.WriteLine("  points on a sphere in 3D.");

        Console.WriteLine("");
        Console.WriteLine("  Radius = " + r + "");

        typeMethods.r8vec_print(3, pc, "  Center:");

        Console.WriteLine("");
        Console.WriteLine("  Number of latitudes is  " + lat_num + "");
        Console.WriteLine("  Number of longitudes is " + long_num + "");

        node_num = Grid_LatLong.sphere_llt_grid_point_count(lat_num, long_num);

        Console.WriteLine("");
        Console.WriteLine("  The number of grid points is " + node_num + "");

        node_xyz = Grid_LatLong.sphere_llt_grid_points(r, pc, lat_num, long_num, node_num);

        Console.WriteLine("");

        k = 0;
        Console.WriteLine("  " + k.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + node_xyz[0 + k * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + node_xyz[1 + k * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + node_xyz[2 + k * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        k += 1;

        Console.WriteLine("");

        for (i = 1; i <= lat_num; i++)
        {
            Console.WriteLine("");
            for (j = 0; j < long_num; j++)
            {
                Console.WriteLine("  " + k.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + node_xyz[0 + k * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + node_xyz[1 + k * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + node_xyz[2 + k * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                k += 1;
                Console.WriteLine("");
            }
        }

        Console.WriteLine("");

        Console.WriteLine("  " + k.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + node_xyz[0 + k * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + node_xyz[1 + k * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + node_xyz[2 + k * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        k += 1;
        Console.WriteLine("");
    }

    private static void sphere_llt_grid_line_count_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_LLT_GRID_LINE_COUNT_TEST tests SPHERE_LLT_GRID_LINE_COUNT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int lat_num;
        int line_num;
        int long_log;
        int long_num;

        lat_num = 3;
        long_num = 4;

        Console.WriteLine("");
        Console.WriteLine("SPHERE_LLT_GRID_LINE_COUNT_TEST");
        Console.WriteLine("  SPHERE_LLT_GRID_LINE_COUNT counts the lines used for a");
        Console.WriteLine("  grid based on triangles defined by latitude and longitude");
        Console.WriteLine("  lines on a sphere in 3D.");
        Console.WriteLine("");
        Console.WriteLine("     LAT_NUM    LONG_NUM   LINE_NUM");

        for (lat_num = 1; lat_num <= 17; lat_num += 2)
        {
            Console.WriteLine("");
            long_num = 1;
            for (long_log = 1; long_log <= 4; long_log++)
            {
                long_num *= 2;
                line_num = Grid_LatLong.sphere_llt_grid_line_count(lat_num, long_num);
                Console.WriteLine("  " + lat_num.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + long_num.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + line_num.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
            }
        }
    }

    private static void sphere_llt_grid_lines_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_LLT_GRID_LINES_TEST tests SPHERE_LLT_GRID_LINES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int lat_num;
        int[] line_data;
        int line_num;
        int long_num;

        lat_num = 3;
        long_num = 4;

        Console.WriteLine("");
        Console.WriteLine("SPHERE_LLT_GRID_LINES_TEST");
        Console.WriteLine("  SPHERE_LLT_GRID_LINES computes grid lines");
        Console.WriteLine("  on a sphere in 3D.");
        Console.WriteLine("");
        Console.WriteLine("  Number of latitudes is  " + lat_num + "");
        Console.WriteLine("  Number of longitudes is " + long_num + "");

        line_num = Grid_LatLong.sphere_llt_grid_line_count(lat_num, long_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of line segments is " + line_num + "");

        line_data = Grid_LatLong.sphere_llt_grid_lines(lat_num, long_num, line_num);

        typeMethods.i4mat_transpose_print(2, line_num, line_data,
            "  Grid line vertices:");
    }

    private static void sphere_llt_grid_display_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_LLT_GRID_DISPLAY_TEST tests SPHERE_LLT_GRID_DISPLAY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int lat_num;
        int[] line_data;
        int line_num;
        int long_num;
        int node_num;
        double[] node_xyz;
        double[] pc = { 0.0, 0.0, 0.0 };
        string prefix;
        double r;

        lat_num = 10;
        long_num = 12;

        r = 10.0;

        Console.WriteLine("");
        Console.WriteLine("SPHERE_LLT_GRID_DISPLAY_TEST");
        Console.WriteLine("  SPHERE_LLT_GRID_DISPLAY displays an LLT grid on a sphere.");
        Console.WriteLine("");
        Console.WriteLine("  Number of latitudes is  " + lat_num + "");
        Console.WriteLine("  Number of longitudes is " + long_num + "");
        //
        //  Get points.
        //
        node_num = Grid_LatLong.sphere_llt_grid_point_count(lat_num, long_num);

        Console.WriteLine("");
        Console.WriteLine("  The number of grid points is " + node_num + "");

        node_xyz = Grid_LatLong.sphere_llt_grid_points(r, pc, lat_num, long_num, node_num);
        //
        //  Get lines.
        //
        line_num = Grid_LatLong.sphere_llt_grid_line_count(lat_num, long_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of line segments is " + line_num + "");

        line_data = Grid_LatLong.sphere_llt_grid_lines(lat_num, long_num, line_num);

        prefix = "sphere_llt_grid";

        Grid_LatLong.sphere_llt_grid_display(node_num, node_xyz, line_num, line_data, prefix);

    }
}