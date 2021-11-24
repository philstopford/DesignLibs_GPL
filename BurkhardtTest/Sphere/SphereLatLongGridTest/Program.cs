using System;
using System.Globalization;
using Burkardt.SphereNS;
using Burkardt.Types;

namespace SphereLatLongGridTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SPHERE_LLQ_GRID_TEST.
        //
        //  Discussion:
        //
        //    SPHERE_LLQ_GRID_TEST tests the SPHERE_LLQ_GRID library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SPHERE_LLQ_GRID_TEST");
        Console.WriteLine("  Test the SPHERE_LLQ_GRID library.");

        sphere_llq_grid_point_count_test();
        sphere_llq_grid_points_test();
        sphere_llq_grid_line_count_test();
        sphere_llq_grid_lines_test();
        sphere_llq_grid_display_test();

        Console.WriteLine("");
        Console.WriteLine("SPHERE_LLQ_GRID_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void sphere_llq_grid_point_count_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_LLQ_GRID_POINT_COUNT_TEST tests SPHERE_LLQ_GRID_POINT_COUNT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int lat_num;

        Console.WriteLine("");
        Console.WriteLine("SPHERE_LLQ_GRID_POINT_COUNT_TEST");
        Console.WriteLine("  SPHERE_LLQ_GRID_POINT_COUNT counts the points used for a");
        Console.WriteLine("  grid based on quadrilaterals defined by latitude and longitude");
        Console.WriteLine("  lines on a sphere in 3D.");
        Console.WriteLine("");
        Console.WriteLine("     LAT_NUM    LONG_NUM   POINT_NUM");

        for (lat_num = 1; lat_num <= 17; lat_num += 2)
        {
            Console.WriteLine("");
            int long_num = 1;
            int long_log;
            for (long_log = 1; long_log <= 4; long_log++)
            {
                long_num *= 2;
                int point_num = Grid_LatLong.sphere_llq_grid_point_count(lat_num, long_num);
                Console.WriteLine("  " + lat_num.ToString().PadLeft(8)
                                       + "  " + long_num.ToString().PadLeft(8)
                                       + "  " + point_num.ToString().PadLeft(8) + "");
            }
        }
    }

    private static void sphere_llq_grid_points_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_LLQ_GRID_POINTS_TEST tests SPHERE_LLQ_GRID_POINTS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        double[] pc = { 0.0, 0.0, 0.0 };

        const int lat_num = 3;
        const int long_num = 4;

        const double r = 10.0;

        Console.WriteLine("");
        Console.WriteLine("SPHERE_LLQ_GRID_POINTS_TEST");
        Console.WriteLine("  SPHERE_LLQ_POINTS produces latitude/longitude");
        Console.WriteLine("  points on a sphere in 3D.");

        Console.WriteLine("");
        Console.WriteLine("  Radius = " + r + "");

        typeMethods.r8vec_print(3, pc, "  Center:");

        Console.WriteLine("");
        Console.WriteLine("  Number of latitudes is  " + lat_num + "");
        Console.WriteLine("  Number of longitudes is " + long_num + "");

        int node_num = Grid_LatLong.sphere_llq_grid_point_count(lat_num, long_num);

        Console.WriteLine("");
        Console.WriteLine("  The number of grid points is " + node_num + "");

        double[] node_xyz = Grid_LatLong.sphere_llq_grid_points(r, pc, lat_num, long_num, node_num);

        Console.WriteLine("");

        int k = 0;
        Console.WriteLine("  " + k.ToString().PadLeft(8)
                               + "  " + node_xyz[0].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + node_xyz[1].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + node_xyz[2].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        k += 1;

        Console.WriteLine("");

        for (i = 1; i <= lat_num; i++)
        {
            Console.WriteLine("");
            int j;
            for (j = 0; j < long_num; j++)
            {
                Console.WriteLine("  " + k.ToString().PadLeft(8)
                                       + "  " + node_xyz[0 + k * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + node_xyz[1 + k * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + node_xyz[2 + k * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                k += 1;
                Console.WriteLine("");
            }
        }

        Console.WriteLine("");

        Console.WriteLine("  " + k.ToString().PadLeft(8)
                               + "  " + node_xyz[0 + k * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + node_xyz[1 + k * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + node_xyz[2 + k * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        k += 1;
        Console.WriteLine("");
    }

    private static void sphere_llq_grid_line_count_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_LLQ_GRID_LINE_COUNT_TEST tests SPHERE_LLQ_GRID_LINE_COUNT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int lat_num;

        Console.WriteLine("");
        Console.WriteLine("SPHERE_LLQ_GRID_LINE_COUNT_TEST");
        Console.WriteLine("  SPHERE_LLQ_GRID_LINE_COUNT counts the lines used for a");
        Console.WriteLine("  grid based on quadrilaterals defined by latitude and longitude");
        Console.WriteLine("  lines on a sphere in 3D.");
        Console.WriteLine("");
        Console.WriteLine("     LAT_NUM    LONG_NUM   LINE_NUM");

        for (lat_num = 1; lat_num <= 17; lat_num += 2)
        {
            Console.WriteLine("");
            int long_num = 1;
            int long_log;
            for (long_log = 1; long_log <= 4; long_log++)
            {
                long_num *= 2;
                int line_num = Grid_LatLong.sphere_llq_grid_line_count(lat_num, long_num);
                Console.WriteLine("  " + lat_num.ToString().PadLeft(8)
                                       + "  " + long_num.ToString().PadLeft(8)
                                       + "  " + line_num.ToString().PadLeft(8) + "");
            }
        }
    }

    private static void sphere_llq_grid_lines_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_LLQ_GRID_LINES_TEST tests SPHERE_LLQ_GRID_LINES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int lat_num = 3;
        const int long_num = 4;

        Console.WriteLine("");
        Console.WriteLine("SPHERE_LLQ_GRID_LINES_TEST");
        Console.WriteLine("  SPHERE_LLQ_GRID_LINES computes grid lines");
        Console.WriteLine("  on a sphere in 3D.");
        Console.WriteLine("");
        Console.WriteLine("  Number of latitudes is  " + lat_num + "");
        Console.WriteLine("  Number of longitudes is " + long_num + "");

        int line_num = Grid_LatLong.sphere_llq_grid_line_count(lat_num, long_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of line segments is " + line_num + "");

        int[] line_data = Grid_LatLong.sphere_llq_grid_lines(lat_num, long_num, line_num);

        typeMethods.i4mat_transpose_print(2, line_num, line_data,
            "  Grid line vertices:");

    }

    private static void sphere_llq_grid_display_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_LLQ_GRID_DISPLAY_TEST tests SPHERE_LLQ_GRID_DISPLAY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] pc = { 0.0, 0.0, 0.0 };

        const int lat_num = 10;
        const int long_num = 12;

        const double r = 10.0;

        Console.WriteLine("");
        Console.WriteLine("SPHERE_LLQ_GRID_DISPLAY_TEST");
        Console.WriteLine("  SPHERE_LLQ_GRID_DISPLAY displays an LLQ grid on a sphere.");
        Console.WriteLine("");
        Console.WriteLine("  Number of latitudes is  " + lat_num + "");
        Console.WriteLine("  Number of longitudes is " + long_num + "");
        //
        //  Get points.
        //
        int node_num = Grid_LatLong.sphere_llq_grid_point_count(lat_num, long_num);

        Console.WriteLine("");
        Console.WriteLine("  The number of grid points is " + node_num + "");

        double[] node_xyz = Grid_LatLong.sphere_llq_grid_points(r, pc, lat_num, long_num, node_num);
        //
        //  Get lines.
        //
        int line_num = Grid_LatLong.sphere_llq_grid_line_count(lat_num, long_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of line segments is " + line_num + "");

        int[] line_data = Grid_LatLong.sphere_llq_grid_lines(lat_num, long_num, line_num);

        string prefix = "sphere_llq_grid";

        Grid_LatLong.sphere_llq_grid_display(node_num, node_xyz, line_num, line_data, prefix);

    }
}