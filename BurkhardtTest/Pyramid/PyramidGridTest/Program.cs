﻿using System;
using Burkardt.PyramidNS;
using Burkardt.Types;

namespace PyramidGridTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for PYRAMID_GRID_TEST.
            //
            //  Discussion:
            //
            //    PYRAMID_GRID_TEST tests the PYRAMID_GRID library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("PYRAMID_GRID_TEST:");
            Console.WriteLine("  Test the PYRAMID_GRID library.");

            test01();
            test02();
            test03();

            Console.WriteLine("");
            Console.WriteLine("PYRAMID_GRID_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests PYRAMID_GRID_SIZE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n;
            int ng;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  PYRAMID_GRID_SIZE determines the size of a");
            Console.WriteLine("  pyramid grid with N+1 points along each edge.");

            Console.WriteLine("");
            Console.WriteLine("   N    Size");
            Console.WriteLine("");
            for (n = 0; n <= 10; n++)
            {
                ng = Grid.pyramid_grid_size(n);
                Console.WriteLine(n.ToString().PadLeft(4) + "  "
                                                          + ng.ToString().PadLeft(6) + "");
            }
        }

        static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests PYRAMID_UNIT_GRID.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n;
            int ng;
            double[] pg;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  PYRAMID_UNIT_GRID determines a unit pyramid");
            Console.WriteLine("  grid with N+1 points along each edge.");

            n = 4;
            typeMethods.r8_print(n, "  Grid parameter N:");

            ng = Grid.pyramid_grid_size(n);
            typeMethods.r8_print(ng, "  Grid size NG:");

            pg = Grid.pyramid_unit_grid(n, ng);

            typeMethods.r8mat_transpose_print(3, ng, pg, "  Pyramid grid points:");
        }

        static void test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 tests PYRAMID_UNIT_GRID_PLOT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string header;
            int n;
            int ng;
            double[] pg;

            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  PYRAMID_UNIT_GRID_PLOT plots a unit pyramid");
            Console.WriteLine("  grid with N+1 points along each edge.");

            n = 5;
            typeMethods.r8_print(n, "  Grid parameter N:");

            ng = Grid.pyramid_grid_size(n);
            typeMethods.r8_print(ng, "  Grid size NG:");

            pg = Grid.pyramid_unit_grid(n, ng);

            header = "pyramid_unit";
            Grid.pyramid_unit_grid_plot(n, ng, pg, header);
        }
    }
}