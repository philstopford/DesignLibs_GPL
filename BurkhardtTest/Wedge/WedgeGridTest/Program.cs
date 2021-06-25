using System;
using System.Collections.Generic;
using System.IO;

using Grid = Burkardt.Wedge.Grid;

namespace WedgeGridTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for WEDGE_GRID_TEST.
            //
            //  Discussion:
            //
            //    WEDGE_GRID_TEST tests the WEDGE_GRID library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("WEDGE_GRID_TEST:");
            Console.WriteLine("  Test the WEDGE_GRID library.");

            test01();
            test02();

            Console.WriteLine("");
            Console.WriteLine("WEDGE_GRID_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests WEDGE_GRID.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] g;
            int j;
            int n = 5;
            int ng;
            string output_filename;
            List<string> output_unit = new List<string>();

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  WEDGE_GRID can define a grid of points");
            Console.WriteLine("  with N+1 points on a side");
            Console.WriteLine("  over the interior of the unit wedge in 3D.");

            Console.WriteLine("");
            Console.WriteLine("  Grid order N = " + n + "");

            ng = Grid.wedge_grid_size(n);

            Console.WriteLine("  Grid count NG = " + ng + "");

            g = Grid.wedge_grid(n, ng);

            Console.WriteLine("");
            Console.WriteLine("     J      X                Y               Z");
            Console.WriteLine("");
            for (j = 0; j < ng; j++)
            {
                Console.WriteLine(j.ToString().PadLeft(6) + "  "
                    + g[0 + j * 3].ToString().PadLeft(14) + "  "
                    + g[1 + j * 3].ToString().PadLeft(14) + "  "
                    + g[2 + j * 3].ToString().PadLeft(14) + "");
            }

            output_filename = "wedge_grid.xy";

            for (j = 0; j < ng; j++)
            {
                output_unit.Add(g[0 + j * 3] + "  "
                    + g[1 + j * 3] + "  "
                    + g[2 + j * 3] + "");
            }
            
            File.WriteAllLines(output_filename, output_unit);

            Console.WriteLine("");
            Console.WriteLine("  Data written to '" + output_filename + "'");
        }

        static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests WEDGE_GRID_PLOT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] g;
            string header;
            int n = 5;
            int ng;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  WEDGE_GRID_PLOT can create GNUPLOT data files");
            Console.WriteLine("  for displaying a wedge grid.");

            Console.WriteLine("");
            Console.WriteLine("  Grid order N = " + n + "");

            ng = Grid.wedge_grid_size(n);

            Console.WriteLine("  Grid count NG = " + ng + "");

            g = Grid.wedge_grid(n, ng);

            header = "wedge";

            Grid.wedge_grid_plot(n, ng, g, header);
        }
    }
}