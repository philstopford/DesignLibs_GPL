using System;
using Burkardt.Types;

namespace EllipsoidGridTest
{
    using Grid = Burkardt.Ellipsoid.Grid;
    
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for ELLIPSOID_GRID_TEST.
            //
            //  Discussion:
            //
            //    ELLIPSOID_GRID_TEST tests the ELLIPSOID_GRID library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 November 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("ELLIPSOID_GRID_TEST:");
            Console.WriteLine("  Test the ELLIPSOID_GRID library.");

            ellipsoid_grid_test01();

            Console.WriteLine("");
            Console.WriteLine("ELLIPSOID_GRID_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void ellipsoid_grid_test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPSOID_GRID_TEST01 tests ELLIPSOID_GRID.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 November 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] c = new double[3];
            string filename = "ellipsoid_grid_test01.xyz";
            int n;
            int ng;
            double[] r = new double[3];
            double[] xyz;

            Console.WriteLine("");
            Console.WriteLine("TEST01:");
            Console.WriteLine("  ELLIPSOID_GRID can define a grid of points");
            Console.WriteLine("  with N+1 points on the minor half axis,");
            Console.WriteLine("  based on any ellipsoid.");

            n = 4;
            r[0] = 2.0;
            r[1] = 1.0;
            r[2] = 1.5;
            c[0] = 1.0;
            c[1] = 2.0;
            c[2] = 1.5;

            Console.WriteLine("");
            Console.WriteLine("  We use N = " + n + "");
            Console.WriteLine("  Radius R = (" + r[0] + "," + r[1] + "," + r[2] + ")");
            Console.WriteLine("  Center C = (" + c[0] + "," + c[1] + "," + c[2] + ")");

            ng = Grid.ellipsoid_grid_count(n, r, c);

            Console.WriteLine("");
            Console.WriteLine("  Number of grid points will be " + ng + "");

            xyz = Grid.ellipsoid_grid(n, r, c, ng);

            typeMethods.r83vec_print_part(ng, xyz, 20, "  Part of the grid point array:");

            typeMethods.r8mat_write(filename, 3, ng, xyz);

            Console.WriteLine("");
            Console.WriteLine("  Data written to the file \"" + filename + "\".");
        }
    }
}