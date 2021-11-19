using System;
using System.Collections.Generic;
using Burkardt.SphereNS;
using Burkardt.Types;

namespace SphereGridTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SPHERE_GRID_TEST.
        //
        //  Discussion:
        //
        //    SPHERE_GRID_TEST tests the SPHERE_GRID library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SPHERE_GRID_TEST");
        Console.WriteLine("  Test the SPHERE_GRID library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();
        test07();
        test08();
        test09();

        test10();
        test11();
        test12();
        test13();

        Console.WriteLine("");
        Console.WriteLine("SPHERE_GRID_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests SPHERE_ICOS_POINT_NUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int edge_num;
        int factor;
        int factor_log;
        int node_num;
        int triangle_num;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  SPHERE_ICOS_POINT_NUM determines the size");
        Console.WriteLine("  (number of vertices, edges and faces) in a grid");
        Console.WriteLine("  on a sphere, made by subdividing an initial");
        Console.WriteLine("  projected icosahedron.");
        Console.WriteLine("");
        Console.WriteLine("  N determines the number of subdivisions of each");
        Console.WriteLine("  edge of the icosahedral faces.");
        Console.WriteLine("");
        Console.WriteLine("         N         V         E         F");
        Console.WriteLine("  --------  --------  --------  --------");
        Console.WriteLine("");

        for (factor = 1; factor <= 20; factor++)
        {
            node_num = Icosphere.sphere_icos_point_num(factor);
            edge_num = Icosphere.sphere_icos_edge_num(factor);
            triangle_num = Icosphere.sphere_icos_face_num(factor);
            Console.WriteLine("  " + factor.ToString().PadLeft(8)
                                   + "  " + node_num.ToString().PadLeft(8)
                                   + "  " + edge_num.ToString().PadLeft(8)
                                   + "  " + triangle_num.ToString().PadLeft(8) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Repeat, but using N constrained by doubling:");
        Console.WriteLine("");
        Console.WriteLine("         N         V         E         F");
        Console.WriteLine("  --------  --------  --------  --------");
        Console.WriteLine("");

        factor = 1;
        for (factor_log = 0; factor_log <= 10; factor_log++)
        {
            node_num = Icosphere.sphere_icos_point_num(factor);
            edge_num = Icosphere.sphere_icos_edge_num(factor);
            triangle_num = Icosphere.sphere_icos_face_num(factor);
            Console.WriteLine("  " + factor.ToString().PadLeft(8)
                                   + "  " + node_num.ToString().PadLeft(8)
                                   + "  " + edge_num.ToString().PadLeft(8)
                                   + "  " + triangle_num.ToString().PadLeft(8) + "");
            factor *= 2;
        }

    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests SPHERE_ICOS1_POINTS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int factor;
        string filename;
        int node_num;
        double[] node_xyz;
        List<string> output = new();

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  SPHERE_ICOS_POINT_NUM \"sizes\" a grid generated");
        Console.WriteLine("  on an icosahedron and projected to a sphere.");
        Console.WriteLine("  SPHERE_ICOS1_POINTS creates the grid points.");

        factor = 3;

        Console.WriteLine("");
        Console.WriteLine("  Sizing factor =       " + factor + "");

        node_num = Icosphere.sphere_icos_point_num(factor);

        Console.WriteLine("");
        Console.WriteLine("  Number of vertices =  " + node_num + "");

        node_xyz = Icosphere.sphere_icos1_points(factor, node_num);

        typeMethods.r8mat_transpose_print_some(3, node_num, node_xyz, 1, 1, 3, 20,
            "  Initial part of NODE_XYZ array:");
        //
        //  Write the nodes to a file.
        //
        filename = "sphere_grid_icos1_points_f"
                   + (factor, "%d") + ".xyz";

        typeMethods.r8mat_write(filename, 3, node_num, node_xyz);

        Console.WriteLine("");
        Console.WriteLine("  Wrote data to \"" + filename + "\"");
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests SPHERE_ICOS2_POINTS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int factor;
        string filename;
        int node_num;
        double[] node_xyz;
        List<string> output = new();

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  SPHERE_ICOS_POINT_NUM \"sizes\" a grid generated");
        Console.WriteLine("  on an icosahedron and projected to a sphere.");
        Console.WriteLine("  SPHERE_ICOS2_POINTS creates the grid.");

        factor = 3;

        Console.WriteLine("");
        Console.WriteLine("  Sizing factor FACTOR = " + factor + "");

        node_num = Icosphere.sphere_icos_point_num(factor);

        Console.WriteLine("");
        Console.WriteLine("  Number of nodes =     " + node_num + "");

        node_xyz = Icosphere.sphere_icos2_points(factor, node_num);

        typeMethods.r8mat_transpose_print_some(3, node_num, node_xyz, 1, 1, 3, 20,
            "  Initial part of NODE_XYZ array:");
        //
        //  Write the nodes to a file.
        //
        filename = "sphere_grid_icos2_points_f" + (factor, "%d")
                                                + ".xyz";

        typeMethods.r8mat_write(filename, 3, node_num, node_xyz);

        Console.WriteLine("");
        Console.WriteLine("  Wrote data to \"" + filename + "\"");

    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests SPHERE_LL_POINTS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int lat_num = 3;
        int lon_num = 4;

        double[] pc = { 0.0, 0.0, 0.0 };
        int i;
        int j;
        int k;
        int node_num;
        double[] node_xyz;
        double r = 10.0;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  SPHERE_LL_POINTS produces latitude/longitude");
        Console.WriteLine("  points on a sphere in 3D.");

        Console.WriteLine("");
        Console.WriteLine("  Radius = " + r + "");

        typeMethods.r8vec_print(3, pc, "  Center:");

        Console.WriteLine("");
        Console.WriteLine("  The number of latitudes =  " + lat_num + "");
        Console.WriteLine("  The number of longitudes = " + lon_num + "");

        node_num = Grid_LatLong.sphere_ll_point_num(lat_num, lon_num);
        Console.WriteLine("");
        Console.WriteLine("  The number of grid points is " + node_num + "");

        node_xyz = Grid_LatLong.sphere_ll_points(r, pc, lat_num, lon_num, node_num);

        Console.WriteLine("");

        k = 0;
        Console.WriteLine("  " + k.ToString().PadLeft(8)
                               + "  " + node_xyz[0 + k * 3].ToString().PadLeft(12)
                               + "  " + node_xyz[1 + k * 3].ToString().PadLeft(12)
                               + "  " + node_xyz[2 + k * 3].ToString().PadLeft(12) + "");

        for (i = 1; i <= lat_num; i++)
        {
            Console.WriteLine("");
            for (j = 0; j < lon_num; j++)
            {
                k += 1;
                Console.WriteLine("  " + k.ToString().PadLeft(8)
                                       + "  " + node_xyz[0 + k * 3].ToString().PadLeft(12)
                                       + "  " + node_xyz[1 + k * 3].ToString().PadLeft(12)
                                       + "  " + node_xyz[2 + k * 3].ToString().PadLeft(12) + "");
            }
        }

        Console.WriteLine("");
        k += 1;
        Console.WriteLine("  " + k.ToString().PadLeft(8)
                               + "  " + node_xyz[0 + k * 3].ToString().PadLeft(12)
                               + "  " + node_xyz[1 + k * 3].ToString().PadLeft(12)
                               + "  " + node_xyz[2 + k * 3].ToString().PadLeft(12) + "");
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests SPHERE_SPIRALPOINTS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        string filename;
        double[] center_xyz = { 0.0, 0.0, 0.0 };
        int node_num = 500;
        double[] node_xyz;
        double r = 1.0;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  SPHERE_SPIRALPOINTS produces a spiral of");
        Console.WriteLine("  points on an implicit sphere in 3D.");
        Console.WriteLine("");
        Console.WriteLine("  Radius = " + r + "");

        typeMethods.r8vec_print(3, center_xyz, "  Center:");

        Console.WriteLine("");
        Console.WriteLine("  The number of spiral points is " + node_num + "");

        node_xyz = Sphere.sphere_spiralpoints(r, center_xyz, node_num);

        typeMethods.r8mat_transpose_print_some(3, node_num, node_xyz, 1, 1, 3, 10,
            "  The spiral points:");
        //
        //  Write the nodes to a file.
        //
        filename = "sphere_grid_spiral_n" + (node_num, "%d") + ".xyz";

        typeMethods.r8mat_write(filename, 3, node_num, node_xyz);

        Console.WriteLine("");
        Console.WriteLine("  Wrote data to \"" + filename + "\"");

    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests SPHERE_LL_LINES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int lat_num = 3;
        int[] line;
        int line_num;
        int long_num = 4;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  SPHERE_LL_LINES computes latitude/longitude");
        Console.WriteLine("  lines on a sphere in 3D.");
        Console.WriteLine("");
        Console.WriteLine("  Number of latitudes is  " + lat_num + "");
        Console.WriteLine("  Number of longitudes is " + long_num + "");

        line_num = Grid_LatLong.sphere_ll_line_num(lat_num, long_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of line segments is " + line_num + "");

        line = Grid_LatLong.sphere_ll_lines(lat_num, long_num, line_num);

        typeMethods.i4mat_transpose_print(2, line_num, line, "  Grid line vertices:");

    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests SPHERE_GRID_Q4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int lat_num = 3;
        int long_num = 4;
        int rectangle_num = lat_num * long_num;
        int[] rectangle_node;

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  SPHERE_GRID_Q4 computes a grid");
        Console.WriteLine("  of Q4 rectangular elements on a sphere in 3D.");
        Console.WriteLine("");
        Console.WriteLine("  Number of latitudes is      " + lat_num + "");
        Console.WriteLine("  Number of longitudes is     " + long_num + "");
        Console.WriteLine("  The number of rectangles is " + rectangle_num + "");

        rectangle_node = Grid.sphere_grid_q4(lat_num, long_num);

        typeMethods.i4mat_transpose_print(4, rectangle_num, rectangle_node,
            "  Rectangle vertices:");

    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests SPHERE_GRID_T3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int lat_num = 3;
        int lon_num = 4;

        int triangle_num;
        int[] triangle_node;

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  SPHERE_GRID_T3 computes a grid");
        Console.WriteLine("  of T3 triangular elements on a sphere in 3D.");
        Console.WriteLine("");
        Console.WriteLine("  Number of latitudes is  " + lat_num + "");
        Console.WriteLine("  Number of longitudes is " + lon_num + "");

        triangle_node = Grid.sphere_grid_t3(lat_num, lon_num);

        triangle_num = 2 * (lat_num + 1) * lon_num;

        Console.WriteLine("");
        Console.WriteLine("  The number of triangles is " + triangle_num + "");

        typeMethods.i4mat_transpose_print(3, triangle_num, triangle_node,
            "  Triangle vertices:");

    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests SPHERE_UNIT_SAMPLE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        string filename;
        int node_num;
        double[] node_xyz;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  For the unit sphere in 3 dimensions:");
        Console.WriteLine("  SPHERE_UNIT_SAMPLE does a random sampling.");

        node_num = 1000;

        node_xyz = Sphere.sphere_unit_sample(node_num, ref seed);

        typeMethods.r8mat_transpose_print_some(3, node_num, node_xyz, 1, 1, 3, 10,
            "  First 10 values:");
        //
        //  Write the nodes to a file.
        //
        filename = "sphere_grid_sample_n" + (node_num, "%d") + ".xyz";

        typeMethods.r8mat_write(filename, 3, node_num, node_xyz);

        Console.WriteLine("");
        Console.WriteLine("  Wrote data to \"" + filename + "\"");

    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tests SPHERE_CUBED_POINTS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        string filename;
        int n;
        int ns;
        double[] xyz;

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  SPHERE_CUBED_POINTS computes points on a cubed sphere grid.");

        n = 10;
        Console.WriteLine("");
        Console.WriteLine("  Number of divisions on each face = " + n + "");

        ns = Cubed.sphere_cubed_point_num(n);
        Console.WriteLine("  Total number of points = " + ns + "");

        xyz = Cubed.sphere_cubed_points(n, ns);

        typeMethods.r8mat_transpose_print_some(3, ns, xyz, 1, 1, 3, 20, "  Initial part of XYZ array:");
        //
        //  Write the nodes to a file.
        //
        filename = "sphere_grid_cubed_f" + (n, "%d") + ".xyz";

        typeMethods.r8mat_write(filename, 3, n, xyz);

        Console.WriteLine("");
        Console.WriteLine("  Wrote data to \"" + filename + "\"");

    }

    private static void test11()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST11 is a dummy file.  See the MATLAB source code for details.
        //
    {
    }

    private static void test12()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST12 is a dummy file.  See the MATLAB source code for details.
        //
    {
    }

    private static void test13()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST13 is a dummy file.  See the MATLAB source code for details.
        //
    {
    }
}