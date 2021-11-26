﻿using System;
using System.Globalization;
using Burkardt.Composition;
using Burkardt.Table;
using Burkardt.TriangleNS;
using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace TriangleExactnessTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGLE_EXACTNESS.
        //
        //  Discussion:
        //
        //    This program investigates the polynomial exactness of a quadrature
        //    rule for the triangle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int degree;
        int degree_max;
        bool error = false;
        int last = 0;
        int point;
        string quad_filename;

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_EXACTNESS");
        Console.WriteLine("");
        Console.WriteLine("  Investigate the polynomial exactness of a quadrature");
        Console.WriteLine("  rule for the triangle by integrating all monomials");
        Console.WriteLine("  of a given degree.");
        Console.WriteLine("");
        Console.WriteLine("  The rule will be adjusted to the unit triangle.");
        //
        //  Get the quadrature file root name:
        //
        try
        {
            quad_filename = args[0];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_EXACTNESS:");
            Console.WriteLine("  Enter the \"root\" name of the quadrature files.");

            quad_filename = Console.ReadLine();
        }

        //
        //  Create the names of:
        //    the quadrature X file;
        //    the quadrature W file;
        //    the quadrature R file;
        //
        string quad_r_filename = quad_filename + "_r.txt";
        string quad_w_filename = quad_filename + "_w.txt";
        string quad_x_filename = quad_filename + "_x.txt";
        //
        //  The second command line argument is the maximum degree.
        //
        try
        {
            degree_max = typeMethods.s_to_i4(args[1], ref last, ref error);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_EXACTNESS:");
            Console.WriteLine("  Please enter the maximum total degree to check.");

            degree_max = Convert.ToInt32(Console.ReadLine());
        }

        //
        //  Summarize the input.
        //
        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_EXACTNESS: User input:");
        Console.WriteLine("  Quadrature rule X file = \"" + quad_x_filename
                                                          + "\".");
        Console.WriteLine("  Quadrature rule W file = \"" + quad_w_filename
                                                          + "\".");
        Console.WriteLine("  Quadrature rule R file = \"" + quad_r_filename
                                                          + "\".");
        Console.WriteLine("  Maximum total degree to check = " + degree_max + "");
        //
        //  Read the X file.
        //
        TableHeader hdr = typeMethods.r8mat_header_read(quad_x_filename);
        int dim_num = hdr.m;
        int point_num = hdr.n;

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points  = " + point_num + "");

        if (dim_num != 2)
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_EXACTNESS - Fatal error!");
            Console.WriteLine("  The quadrature abscissas must be two dimensional.");
            return;
        }

        double[] x = typeMethods.r8mat_data_read(quad_x_filename, dim_num, point_num);
        //
        //  Read the W file.
        //
        hdr = typeMethods.r8mat_header_read(quad_w_filename);
        int dim_num2 = hdr.m;
        int point_num2 = hdr.n;

        if (dim_num2 != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_EXACTNESS - Fatal error!");
            Console.WriteLine("  The quadrature weight file should have exactly");
            Console.WriteLine("  one value on each line.");
            return;
        }

        if (point_num2 != point_num)
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_EXACTNESS - Fatal error!");
            Console.WriteLine("  The quadrature weight file should have exactly");
            Console.WriteLine("  the same number of lines as the abscissa file.");
            return;
        }

        double[] w = typeMethods.r8mat_data_read(quad_w_filename, 1, point_num);
        //
        //  Read the R file.
        //
        hdr = typeMethods.r8mat_header_read(quad_r_filename);
        int dim_num3 = hdr.m;
        int point_num3 = hdr.n;

        if (dim_num3 != dim_num)
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_EXACTNESS - Fatal error!");
            Console.WriteLine("  The quadrature region file should have the same");
            Console.WriteLine("  number of values on each line as the abscissa file");
            Console.WriteLine("  does.");
            return;
        }

        if (point_num3 != 3)
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_EXACTNESS - Fatal error!");
            Console.WriteLine("  The quadrature region file should have 3 lines.");
            return;
        }

        double[] r = typeMethods.r8mat_data_read(quad_r_filename, dim_num, 3);
        //
        //  Rescale the weights.
        //
        double area = Integrals.triangle_area(r);

        for (point = 0; point < point_num; point++)
        {
            w[point] = 0.5 * w[point] / area;
        }

        //
        //  Translate the abscissas.
        //
        double[] x_ref = new double[2 * point_num];

        Triangulation.triangle_order3_physical_to_reference(r, point_num, x, ref x_ref);
        //
        //  Explore the monomials.
        //
        int[] expon = new int[dim_num];

        Console.WriteLine("");
        Console.WriteLine("      Error    Degree  Exponents");

        for (degree = 0; degree <= degree_max; degree++)
        {
            Console.WriteLine("");
            bool more = false;
            int h = 0;
            int t = 0;

            for (;;)
            {
                Comp.comp_next(degree, dim_num, ref expon, ref more, ref h, ref t);

                double quad_error = Integrals.triangle01_monomial_quadrature(dim_num, expon, point_num,
                    x_ref, w);

                string cout = "  " + quad_error.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "     " + degree.ToString().PadLeft(2)
                                   + "  ";

                int dim;
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += expon[dim].ToString().PadLeft(3);
                }

                Console.WriteLine(cout);

                if (!more)
                {
                    break;
                }
            }
        }

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_EXACTNESS:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}