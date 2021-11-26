using System;
using System.Globalization;
using Burkardt.Composition;
using Burkardt.MonomialNS;
using Burkardt.Table;
using Burkardt.Types;
using Burkardt.Wedge;

namespace WedgeExactnessTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for WEDGE_EXACTNESS.
        //
        //  Discussion:
        //
        //    This program investigates the polynomial exactness of a quadrature
        //    rule for the unit wedge.
        //
        //    The interior of the unit wedge in 3D is defined by the constraints:
        //      0 <= X
        //      0 <= Y
        //           X + Y <= 1
        //     -1 <= Z <= +1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int degree;
        int degree_max;
        string quad_filename;

        Console.WriteLine("");
        Console.WriteLine("WEDGE_EXACTNESS");
        Console.WriteLine("  Investigate the polynomial exactness of a quadrature");
        Console.WriteLine("  rule for the unit wedge.");
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
            Console.WriteLine("WEDGE_EXACTNESS:");
            Console.WriteLine("  Enter the \"root\" name of the quadrature files.");

            quad_filename = Console.ReadLine();
        }

        //
        //  Create the names of:
        //    the quadrature X file;
        //    the quadrature W file;
        //
        string quad_w_filename = quad_filename + "_w.txt";
        string quad_x_filename = quad_filename + "_x.txt";
        //
        //  The second command line argument is the maximum degree.
        //
        try
        {
            degree_max = Convert.ToInt32(args[1]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("WEDGE_EXACTNESS:");
            Console.WriteLine("  Please enter the maximum total degree to check.");

            degree_max = Convert.ToInt32(Console.ReadLine());
        }

        //
        //  Summarize the input.
        //
        Console.WriteLine("");
        Console.WriteLine("WEDGE_EXACTNESS: User input:");
        Console.WriteLine("  Quadrature rule X file = \"" + quad_x_filename + "\".");
        Console.WriteLine("  Quadrature rule W file = \"" + quad_w_filename + "\".");
        Console.WriteLine("  Maximum total degree to check = " + degree_max + "");
        //
        //  Read the X file.
        //
        TableHeader hd = typeMethods.r8mat_header_read(quad_x_filename);
        int dim_num = hd.m;
        int order = hd.n;

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points  = " + order + "");

        if (dim_num != 3)
        {
            Console.WriteLine("");
            Console.WriteLine("WEDGE_EXACTNESS - Fatal error!");
            Console.WriteLine("  The quadrature abscissas must be 3 dimensional.");
            return;
        }

        double[] x = typeMethods.r8mat_data_read(quad_x_filename, dim_num, order);
        //
        //  Read the W file.
        //
        hd = typeMethods.r8mat_header_read(quad_w_filename);
        int dim_num2 = hd.m;
        int order2 = hd.n;

        if (dim_num2 != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("WEDGE_EXACTNESS - Fatal error!");
            Console.WriteLine("  The quadrature weight file should have exactly");
            Console.WriteLine("  one value on each line.");
            return;
        }

        if (order2 != order)
        {
            Console.WriteLine("");
            Console.WriteLine("WEDGE_EXACTNESS - Fatal error!");
            Console.WriteLine("  The quadrature weight file should have exactly");
            Console.WriteLine("  the same number of lines as the abscissa file.");
            return;
        }

        double[] w = typeMethods.r8mat_data_read(quad_w_filename, 1, order);
        //
        //  Explore the monomials.
        //
        int[] expon = new int[dim_num];

        Console.WriteLine("");
        Console.WriteLine("      Error    Degree  Exponents");
        Console.WriteLine("");

        for (degree = 0; degree <= degree_max; degree++)
        {
            bool more = false;
            int h = 0;
            int t = 0;

            for (;;)
            {
                Comp.comp_next(degree, dim_num, ref expon, ref more, ref h, ref t);

                double[] v = Monomial.monomial_value(dim_num, order, expon, x);

                double quad = Integrals.wedge01_volume() * typeMethods.r8vec_dot_product(order, w, v);

                double exact = Integrals.wedge01_integral(expon);

                double quad_error = Math.Abs(quad - exact);

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

            Console.WriteLine("");
        }

        Console.WriteLine("");
        Console.WriteLine("WEDGE_EXACTNESS:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

}