using System;
using Burkardt.Composition;
using Burkardt.Quadrature;
using Burkardt.Table;
using Burkardt.Types;

namespace nIntExactnessTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for NINT_EXACTNESS_MIXED.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] alpha;
        double[] beta;
        int degree;
        int degree_max;
        int dim;
        int dim_num = 0;
        int dim_num2;
        bool error = false;
        int[] expon;
        int h;
        int last = 0;
        bool more;
        int point_num = 0;
        int point_num2;
        double quad_error;
        string quad_filename;
        string quad_a_filename;
        string quad_b_filename;
        string quad_r_filename;
        string quad_w_filename;
        string quad_x_filename;
        int[] rule;
        int t;
        double[] weight;
        double[] x;
        double[] x_range;

        Console.WriteLine("");
        Console.WriteLine("NINT_EXACTNESS_MIXED");
        Console.WriteLine("");
        Console.WriteLine("  Investigate the polynomial exactness of");
        Console.WriteLine("  a multidimensional quadrature rule");
        Console.WriteLine("  for a region R = R1 x R2 x ... x RM.");
        Console.WriteLine("");
        Console.WriteLine("  Individual rules may be for:");
        Console.WriteLine("");
        Console.WriteLine("  Legendre:");
        Console.WriteLine("  region: [-1,+1]");
        Console.WriteLine("  weight: w(x)=1");
        Console.WriteLine("  rules: Gauss-Legendre, Clenshaw-Curtis, Fejer2, Gauss-Patterson");
        Console.WriteLine("");
        Console.WriteLine("  Jacobi:");
        Console.WriteLine("  region: [-1,+1]");
        Console.WriteLine("  weight: w(x)=(1-x)^alpha (1+x)^beta");
        Console.WriteLine("  rules: Gauss-Jacobi");
        Console.WriteLine("");
        Console.WriteLine("  Laguerre:");
        Console.WriteLine("  region: [0,+oo)");
        Console.WriteLine("  weight: w(x)=exp(-x)");
        Console.WriteLine("  rules: Gauss-Laguerre");
        Console.WriteLine("");
        Console.WriteLine("  Generalized Laguerre:");
        Console.WriteLine("  region: [0,+oo)");
        Console.WriteLine("  weight: w(x)=x^alpha exp(-x)");
        Console.WriteLine("  rules: Generalized Gauss-Laguerre");
        Console.WriteLine("");
        Console.WriteLine("  Hermite:");
        Console.WriteLine("  region: (-oo,+o)");
        Console.WriteLine("  weight: w(x)=exp(-x*x)");
        Console.WriteLine("  rules: Gauss-Hermite");
        Console.WriteLine("");
        Console.WriteLine("  Generalized Hermite:");
        Console.WriteLine("  region: (-oo,+oo)");
        Console.WriteLine("  weight: w(x)=|x|^alpha exp(-x*x)");
        Console.WriteLine("  rules: generalized Gauss-Hermite");
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
            Console.WriteLine("NINT_EXACTNESS_MIXED:");
            Console.WriteLine("  Enter the \"root\" name of the quadrature files.");

            quad_filename = Console.ReadLine();
        }

        //
        //  Create the names of:
        //    the quadrature X file;
        //    the quadrature W file;
        //    the quadrature R file;
        //    the output "exactness" file.
        //
        quad_a_filename = quad_filename + "_a.txt";
        quad_b_filename = quad_filename + "_b.txt";
        quad_r_filename = quad_filename + "_r.txt";
        quad_w_filename = quad_filename + "_w.txt";
        quad_x_filename = quad_filename + "_x.txt";
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
            Console.WriteLine("NINT_EXACTNESS_MIXED:");
            Console.WriteLine("  Please enter the maximum total degree to check.");

            degree_max = Convert.ToInt32(Console.ReadLine());
        }

        //
        //  Summarize the input.
        //
        Console.WriteLine("");
        Console.WriteLine("NINT_EXACTNESS: User input:");
        Console.WriteLine("  Quadrature rule A file = \"" + quad_a_filename + "\".");
        Console.WriteLine("  Quadrature rule B file = \"" + quad_b_filename + "\".");
        Console.WriteLine("  Quadrature rule R file = \"" + quad_r_filename + "\".");
        Console.WriteLine("  Quadrature rule W file = \"" + quad_w_filename + "\".");
        Console.WriteLine("  Quadrature rule X file = \"" + quad_x_filename + "\".");
        Console.WriteLine("  Maximum total degree to check = " + degree_max + "");
        //
        //  Read the X file.
        //
        TableHeader th = typeMethods.r8mat_header_read(quad_x_filename);
        dim_num = th.m;
        point_num = th.n;

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points  = " + point_num + "");

        x = typeMethods.r8mat_data_read(quad_x_filename, dim_num, point_num);
        //
        //  Read the W file.
        //
        th = typeMethods.r8mat_header_read(quad_w_filename);
        dim_num2 = th.m;
        point_num2 = th.n;

        if (dim_num2 != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("NINT_EXACTNESS_MIXED - Fatal error!");
            Console.WriteLine("  The quadrature weight file should have exactly");
            Console.WriteLine("  one value on each line.");
            return;
        }

        if (point_num2 != point_num)
        {
            Console.WriteLine("");
            Console.WriteLine("NINT_EXACTNESS_MIXED - Fatal error!");
            Console.WriteLine("  The quadrature weight file should have exactly");
            Console.WriteLine("  the same number of lines as the abscissa file.");
            return;
        }

        weight = typeMethods.r8mat_data_read(quad_w_filename, 1, point_num);
        //
        //  Read the R file.
        //
        th = typeMethods.r8mat_header_read(quad_r_filename);
        dim_num2 = th.m;
        point_num2 = th.n;

        if (dim_num2 != dim_num)
        {
            Console.WriteLine("");
            Console.WriteLine("NINT_EXACTNESS_MIXED - Fatal error!");
            Console.WriteLine("  The quadrature region file should have the same");
            Console.WriteLine("  number of values on each line as the abscissa file.");
            return;
        }

        if (point_num2 != 2)
        {
            Console.WriteLine("");
            Console.WriteLine("NINT_EXACTNESS_MIXED - Fatal error!");
            Console.WriteLine("  The quadrature region file should have two lines.");
            return;
        }

        x_range = beta = typeMethods.r8mat_data_read(quad_r_filename, dim_num, 2);
        //
        //  Read the A file.
        //
        th = typeMethods.r8mat_header_read(quad_a_filename);
        dim_num2 = th.m;
        point_num2 = th.n;

        if (dim_num2 != dim_num)
        {
            Console.WriteLine("");
            Console.WriteLine("NINT_EXACTNESS_MIXED - Fatal error!");
            Console.WriteLine("  The quadrature A file should have the same");
            Console.WriteLine("  number of values on each line as the abscissa file.");
            return;
        }

        if (point_num2 != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("NINT_EXACTNESS_MIXED - Fatal error!");
            Console.WriteLine("  The quadrature A file should have 1 line.");
            return;
        }

        alpha = beta = typeMethods.r8mat_data_read(quad_a_filename, dim_num, 1);
        //
        //  Read the B file.
        //
        th = typeMethods.r8mat_header_read(quad_b_filename);
        dim_num2 = th.m;
        point_num2 = th.n;

        if (dim_num2 != dim_num)
        {
            Console.WriteLine("");
            Console.WriteLine("NINT_EXACTNESS_MIXED - Fatal error!");
            Console.WriteLine("  The quadrature B file should have the same");
            Console.WriteLine("  number of values on each line as the abscissa file,");
            return;
        }

        if (point_num2 != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("NINT_EXACTNESS_MIXED - Fatal error!");
            Console.WriteLine("  The quadrature B file should have 1 line.");
            return;
        }

        beta = typeMethods.r8mat_data_read(quad_b_filename, dim_num, 1);
        //
        //  Try to determine the rule types.
        //
        rule = new int[dim_num];

        Console.WriteLine("");
        Console.WriteLine("  Analysis of integration region:");
        Console.WriteLine("");

        for (dim = 0; dim < dim_num; dim++)
        {
            switch (x_range[dim + 0 * dim_num])
            {
                case -1.0 when x_range[dim + 1 * dim_num] == +1.0:
                {
                    switch (alpha[dim])
                    {
                        case 0.0 when beta[dim] == 0.0:
                            rule[dim] = 1;
                            Console.WriteLine("  " + dim.ToString().PadLeft(4) + "  Gauss Legendre.");
                            break;
                        default:
                            rule[dim] = 2;
                            Console.WriteLine("  " + dim.ToString().PadLeft(4) + "  Gauss Jacobi"
                                              + ", ALPHA = " + alpha[dim]
                                              + ", BETA = " + beta[dim] + "");
                            break;
                    }

                    break;
                }
                case 0.0 when x_range[dim + 1 * dim_num] == typeMethods.r8_huge():
                {
                    switch (alpha[dim])
                    {
                        case 0.0:
                            rule[dim] = 3;
                            Console.WriteLine("  " + dim.ToString().PadLeft(4) + "  Gauss Laguerre.");
                            break;
                        default:
                            rule[dim] = 4;
                            Console.WriteLine("  " + dim.ToString().PadLeft(4) + "  Generalized Gauss Laguerre"
                                              + ", ALPHA = " + alpha[dim] + "");
                            break;
                    }

                    break;
                }
                default:
                {
                    if (x_range[dim + 0 * dim_num] == -typeMethods.r8_huge() &&
                        x_range[dim + 1 * dim_num] == +typeMethods.r8_huge())
                    {
                        switch (alpha[dim])
                        {
                            case 0.0:
                                rule[dim] = 5;
                                Console.WriteLine("  " + dim.ToString().PadLeft(4) + "  Gauss Hermite dimension.");
                                break;
                            default:
                                rule[dim] = 6;
                                Console.WriteLine("  " + dim.ToString().PadLeft(4) + "  Generalized Gauss Hermite"
                                                  + ", ALPHA = " + alpha[dim] + "");
                                break;
                        }
                    }
                    else
                    {
                        Console.WriteLine("");
                        Console.WriteLine("NINT_EXACTNESS_MIXED - Fatal error!");
                        Console.WriteLine("  Did not recognize region component.");
                        Console.WriteLine("  Dimension DIM = " + dim + "");
                        Console.WriteLine("  A = " + x_range[dim + 0 * dim_num] + "");
                        Console.WriteLine("  B = " + x_range[dim + 1 * dim_num] + "");
                        return;
                    }

                    break;
                }
            }
        }

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

                quad_error = MonomialQuadrature.monomial_quadrature(dim_num, point_num, rule, alpha,
                    beta, expon, weight, x);

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

            Console.WriteLine("");
            Console.WriteLine("NINT_EXACTNESS_MIXED:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}