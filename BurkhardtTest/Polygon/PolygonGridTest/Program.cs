using System;
using Burkardt.Types;

namespace PolygonGridTest
{
    using Grid = Burkardt.Polygon.Grid;

    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for POLYGON_GRID_TEST.
            //
            //  Discussion:
            //
            //    POLYGON_GRID_TEST tests POLYGON_GRID.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("POLYGON_GRID_TEST:");
            Console.WriteLine("  Test the POLYGON_GRID library.");

            polygon_grid_count_test();

            polygon_grid_points_test01();
            polygon_grid_points_test02();
            polygon_grid_points_test03();

            Console.WriteLine("");
            Console.WriteLine("POLYGON_GRID_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void polygon_grid_count_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_GRID_COUNT_TEST tests POLYGON_GRID_COUNT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n;
            int ng;
            int nv;

            Console.WriteLine("");
            Console.WriteLine("POLYGON_GRID_COUNT_TEST:");
            Console.WriteLine("  POLYGON_GRID_COUNT counts NG, the number of points in");
            Console.WriteLine("  a grid defined with N+1 points on each side of a");
            Console.WriteLine("  polygon of NV vertices.");

            for (nv = 3; nv <= 5; nv++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Polygonal vertex count NV = " + nv + "");
                Console.WriteLine("");
                Console.WriteLine("   N     NG");
                Console.WriteLine("");
                for (n = 0; n <= 5; n++)
                {
                    ng = Grid.polygon_grid_count(n, nv);
                    Console.WriteLine("  " + n.ToString().PadLeft(2)
                                           + "  " + ng.ToString().PadLeft(5) + "");
                }
            }
        }

        static void polygon_grid_points_test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_GRID_POINTS_TEST01 tests POLYGON_GRID_POINTS
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string filename;
            int n;
            int ng;
            int nv = 3;
            string prefix;
            double[] v =
            {
                0.0, 0.0,
                1.0, 0.0,
                0.5, 0.86602540378443860
            };
            double[] xg;

            Console.WriteLine("");
            Console.WriteLine("POLYGON_GRID_POINTS_TEST01:");
            Console.WriteLine("  POLYGON_GRID_POINTS returns grid points for a polygon");
            Console.WriteLine("  of NV vertices, with N+1 points on a side");
            Console.WriteLine("");
            Console.WriteLine("  For this test, the polygon is a triangle.");

            typeMethods.r8mat_transpose_print(2, nv, v, "  Polygon vertices:");
            //
            //  Count the grid points.
            //
            n = 5;
            ng = Grid.polygon_grid_count(n, nv);

            Console.WriteLine("");
            Console.WriteLine("  N = " + n + "");
            Console.WriteLine("  Number of grid points will be NG = " + ng + "");
            //
            //  Compute the grid points.
            //
            xg = Grid.polygon_grid_points(n, nv, v, ng);

            typeMethods.r8mat_transpose_print(2, ng, xg, "  The grid point array:");
            //
            //  Display the points.
            //
            prefix = "triangle";

            Grid.polygon_grid_display(n, nv, v, ng, xg, prefix);
            //
            //  Write the points to a file.
            //
            filename = "triangle.xy";

            typeMethods.r8mat_write(filename, 2, ng, xg);

            Console.WriteLine("");
            Console.WriteLine("  Data written to the file '" + filename + "'");
        }

        static void polygon_grid_points_test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_GRID_POINTS_TEST02 tests POLYGON_GRID_POINTS
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string filename;
            int n;
            int ng;
            int nv = 4;
            string prefix;
            double[] v =
            {
                1.0, 1.0,
                2.0, 0.0,
                4.0, 3.0,
                0.0, 5.0
            };
            double[] xg;

            Console.WriteLine("");
            Console.WriteLine("POLYGON_GRID_POINTS_TEST02:");
            Console.WriteLine("  POLYGON_GRID_POINTS returns grid points for a polygon");
            Console.WriteLine("  of NV vertices, with N+1 points on a side");
            Console.WriteLine("");
            Console.WriteLine("  For this test, the polygon is a convex quadrilateral");
            Console.WriteLine("  with sides of varying length.");
            //
            //  Define the polygon.
            //
            typeMethods.r8mat_transpose_print(2, nv, v, "  Polygon vertices:");
            //
            //  Count the grid points.
            //
            n = 7;
            ng = Grid.polygon_grid_count(n, nv);

            Console.WriteLine("");
            Console.WriteLine("  N = " + n + "");
            Console.WriteLine("  Number of grid points will be NG = " + ng + "");
            //
            //  Compute the grid points.
            //
            xg = Grid.polygon_grid_points(n, nv, v, ng);

            typeMethods.r8mat_transpose_print(2, ng, xg, "  The grid point array:");
            //
            //  Display the points.
            //
            prefix = "quad";

            Grid.polygon_grid_display(n, nv, v, ng, xg, prefix);
            //
            //  Write the points to a file.
            //
            filename = "quad.xy";

            typeMethods.r8mat_write(filename, 2, ng, xg);

            Console.WriteLine("");
            Console.WriteLine("  Data written to the file '" + filename + "'");

        }

        static void polygon_grid_points_test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_GRID_POINTS_TEST03 tests POLYGON_GRID_POINTS
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string filename;
            int n;
            int ng;
            int nv = 6;
            string prefix;
            double[] v =
            {
                0.0, 0.0,
                2.0, 0.0,
                2.0, 1.0,
                1.0, 1.0,
                1.0, 2.0,
                0.0, 2.0
            };
            double[] xg;

            Console.WriteLine("");
            Console.WriteLine("POLYGON_GRID_POINTS_TEST03:");
            Console.WriteLine("  POLYGON_GRID_POINTS returns grid points for a polygon");
            Console.WriteLine("  of NV vertices, with N+1 points on a side");
            Console.WriteLine("");
            Console.WriteLine("  For this test, the polygon is nonconvex and six sided.");
            Console.WriteLine("  Two degenerate triangles are created, and some grid points");
            Console.WriteLine("  are generated several times.");
            //
            //  Define the polygon.
            //
            typeMethods.r8mat_transpose_print(2, nv, v, "  Polygon vertices:");
            //
            //  Count the grid points.
            //
            n = 5;
            ng = Grid.polygon_grid_count(n, nv);

            Console.WriteLine("");
            Console.WriteLine("  N = " + n + "");
            Console.WriteLine("  Number of grid points will be NG = " + ng + "");
            //
            //  Compute the grid points.
            //
            xg = Grid.polygon_grid_points(n, nv, v, ng);

            typeMethods.r8mat_transpose_print(2, ng, xg, "  The grid point array:");
            //
            //  Display the points.
            //
            prefix = "ell";

            Grid.polygon_grid_display(n, nv, v, ng, xg, prefix);
            //
            //  Write the points to a file.
            //
            filename = "ell.xy";

            typeMethods.r8mat_write(filename, 2, ng, xg);

            Console.WriteLine("");
            Console.WriteLine("  Data written to the file '" + filename + "'");
        }
    }
}