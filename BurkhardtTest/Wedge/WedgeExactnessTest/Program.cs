using System;
using Burkardt.Composition;
using Burkardt.MonomialNS;
using Burkardt.Table;
using Burkardt.Types;
using Burkardt.Wedge;

namespace WedgeExactnessTest
{
    class Program
    {
        static void Main(string[] args)
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
            int dim;
            int dim_num;
            int dim_num2;
            bool error;
            double exact;
            int[] expon;
            int h;
            int last;
            bool more;
            int order;
            int order2;
            double quad;
            double quad_error;
            string quad_filename;
            string quad_w_filename;
            string quad_x_filename;
            int t;
            double[] v;
            double[] w;
            double[] x;

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
            quad_w_filename = quad_filename + "_w.txt";
            quad_x_filename = quad_filename + "_x.txt";
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
            dim_num = hd.m;
            order = hd.n;

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

            x = typeMethods.r8mat_data_read(quad_x_filename, dim_num, order);
            //
            //  Read the W file.
            //
            hd = typeMethods.r8mat_header_read(quad_w_filename);
            dim_num2 = hd.m;
            order2 = hd.n;

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

            w = typeMethods.r8mat_data_read(quad_w_filename, 1, order);
            //
            //  Explore the monomials.
            //
            expon = new int[dim_num];

            Console.WriteLine("");
            Console.WriteLine("      Error    Degree  Exponents");
            Console.WriteLine("");

            for (degree = 0; degree <= degree_max; degree++)
            {
                more = false;
                h = 0;
                t = 0;

                for (;;)
                {
                    Comp.comp_next(degree, dim_num, ref expon, ref more, ref h, ref t);

                    v = Monomial.monomial_value(dim_num, order, expon, x);

                    quad = Integrals.wedge01_volume() * typeMethods.r8vec_dot_product(order, w, v);

                    exact = Integrals.wedge01_integral(expon);

                    quad_error = Math.Abs(quad - exact);

                    string cout = "  " + quad_error.ToString().PadLeft(12)
                                       + "     " + degree.ToString().PadLeft(2)
                                       + "  ";

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
}