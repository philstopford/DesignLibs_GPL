﻿using System;
using System.Collections.Generic;
using System.IO;
using Grid = Burkardt.Triangle.Grid;

namespace TriangleGridTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TRIANGLE_GRID_TEST.
            //
            //  Discussion:
            //
            //    TRIANGLE_GRID_TEST tests the TRIANGLE_GRID library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_GRID_TEST:");
            Console.WriteLine("  C++ version");
            Console.WriteLine("  Test the TRIANGLE_GRID library.");

            test01();
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_GRID_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests TRIANGLE_GRID.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n = 10;
            int ng = ((n + 1) * (n + 2)) / 2;

            string filename;
            int j;
            List<string> output = new List<string>();
            double[] t =  {
                0.0, 0.0,
                1.0, 0.0,
                0.5, 0.86602540378443860
            }
            ;
            double[] tg;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  TRIANGLE_GRID can define a triangular grid of points");
            Console.WriteLine("  with N+1 points on a side, based on any triangle.");

            Console.WriteLine("");
            Console.WriteLine("  Defining triangle:");
            Console.WriteLine("     J      X      Y");
            Console.WriteLine("");
            for (j = 0; j < 3; j++)
            {
                Console.WriteLine("  " + j.ToString().PadLeft(4)
                    + "  " + t[0 + j * 2].ToString().PadLeft(12)
                    + "  " + t[1 + j * 2].ToString().PadLeft(12) + "");
            }

            tg = Grid.triangle_grid(n, t);

            Console.WriteLine("");
            Console.WriteLine("     J      X      Y");
            Console.WriteLine("");
            for (j = 0; j < ng; j++)
            {
                Console.WriteLine("  " + j.ToString().PadLeft(4)
                    + "  " + tg[0 + j * 2].ToString().PadLeft(12)
                    + "  " + tg[1 + j * 2].ToString().PadLeft(12) + "");
            }

            filename = "triangle_grid_test01.xy";

            for (j = 0; j < ng; j++)
            {
                output.Add("  " + tg[0 + j * 2].ToString().PadLeft(12)
                    + "  " + tg[1 + j * 2].ToString().PadLeft(12) + "");
            }

            File.WriteAllLines(filename, output);

            Console.WriteLine("");
            Console.WriteLine("  Data written to \"" + filename + "\"");
        }
    }
}