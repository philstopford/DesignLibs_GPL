using System;
using Burkardt.Quadrature;
using Burkardt.Table;
using Burkardt.Types;

namespace PowerRuleTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for POWER_RULE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim_num;
        bool error = false;
        int last = 0;
        string quad_1d_filename;

        Console.WriteLine("");
        Console.WriteLine("POWER_RULE");
        Console.WriteLine("");
        Console.WriteLine("  Create a multidimensional power rule");
        Console.WriteLine("  as a product of identical 1D integration rules.");
        //
        //  Get the quadrature file root name:
        //
        try
        {
            quad_1d_filename = args[0];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("POWER_RULE:");
            Console.WriteLine("  Enter the \"root\" name of the 1D quadrature files.");

            quad_1d_filename = Console.ReadLine();
        }

        //
        //  Create the names of:
        //    the quadrature X file;
        //    the quadrature W file;
        //    the quadrature R file;
        //
        string quad_x_1d_filename = quad_1d_filename + "_x.txt";
        string quad_w_1d_filename = quad_1d_filename + "_w.txt";
        string quad_r_1d_filename = quad_1d_filename + "_r.txt";
        //
        //  The second command line argument is the spatial dimension.
        //
        try
        {
            dim_num = typeMethods.s_to_i4(args[1], ref last, ref error);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("POWER_RULE:");
            Console.WriteLine("  Please enter the desired spatial dimension of the rule.");

            dim_num = Convert.ToInt32(Console.ReadLine());
        }

        //
        //  Summarize the input.
        //
        Console.WriteLine("");
        Console.WriteLine("POWER_RULE: User input:");
        Console.WriteLine("  Quadrature rule X file = \"" + quad_x_1d_filename
                                                          + "\".");
        Console.WriteLine("  Quadrature rule W file = \"" + quad_w_1d_filename
                                                          + "\".");
        Console.WriteLine("  Quadrature rule R file = \"" + quad_r_1d_filename
                                                          + "\".");
        Console.WriteLine("  Spatial dimension = " + dim_num + "");
        //
        //  Read the X file.
        //
        TableHeader h = typeMethods.r8mat_header_read(quad_x_1d_filename);
        int dim_num_1d = h.m;
        int point_num_1d = h.n;

        if (dim_num_1d != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("POWER_RULE - Fatal error!");
            Console.WriteLine("  The 1D quadrature abscissa file should have exactly");
            Console.WriteLine("  one value on each line.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("  Number of points in 1D rule = " + point_num_1d + "");

        double[] x_1d = typeMethods.r8mat_data_read(quad_x_1d_filename, dim_num_1d, point_num_1d);
        //
        //  Read the W file.
        //
        h = typeMethods.r8mat_header_read(quad_w_1d_filename);

        dim_num_1d = h.m;
        int point_num_1d2 = h.n;

        if (dim_num_1d != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("POWER_RULE - Fatal error!");
            Console.WriteLine("  The 1D quadrature weight file should have exactly");
            Console.WriteLine("  one value on each line.");
            return;
        }

        if (point_num_1d2 != point_num_1d)
        {
            Console.WriteLine("");
            Console.WriteLine("POWER_RULE - Fatal error!");
            Console.WriteLine("  The 1D quadrature weight file should have exactly");
            Console.WriteLine("  the same number of lines as the abscissa file.");
            return;
        }

        double[] w_1d = typeMethods.r8mat_data_read(quad_w_1d_filename, dim_num_1d, point_num_1d);
        //
        //  Read the R file.
        //
        h = typeMethods.r8mat_header_read(quad_r_1d_filename);

        dim_num_1d = h.m;
        point_num_1d2 = h.n;

        if (dim_num_1d != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("POWER_RULE - Fatal error!");
            Console.WriteLine("  The 1D quadrature region file should have exactly");
            Console.WriteLine("  one value on each line.");
            return;
        }

        if (point_num_1d2 != 2)
        {
            Console.WriteLine("");
            Console.WriteLine("POWER_RULE - Fatal error!");
            Console.WriteLine("  The 1D quadrature region file should have two lines.");
            return;
        }

        double[] r_1d = typeMethods.r8mat_data_read(quad_r_1d_filename, 1, 2);
        //
        //  Determine size of the rule.
        //
        int point_num = PowerQuadrature.power_rule_size(point_num_1d, dim_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of points in rule = " + point_num + "");
        //
        //  Compute the rule.
        //
        double[] w = new double [point_num];
        double[] x = new double [dim_num * point_num];
        double[] r = new double [dim_num * 2];

        PowerQuadrature.power_rule_set(point_num_1d, x_1d, w_1d, r_1d, dim_num, point_num,
            ref x, ref w, ref r);
        //
        //  Write rule to files.
        //
        string quad_x_filename = "power_x.txt";
        string quad_w_filename = "power_w.txt";
        string quad_r_filename = "power_r.txt";

        Console.WriteLine("");
        Console.WriteLine("  Creating quadrature rule X file = \""
                          + quad_x_filename + "\".");

        typeMethods.r8mat_write(quad_x_filename, dim_num, point_num, x);

        Console.WriteLine("  Creating quadrature rule W file = \""
                          + quad_w_filename + "\".");

        typeMethods.r8mat_write(quad_w_filename, 1, point_num, w);

        Console.WriteLine("  Creating quadrature rule R file = \""
                          + quad_r_filename + "\".");

        typeMethods.r8mat_write(quad_r_filename, dim_num, 2, r);

        Console.WriteLine("");
        Console.WriteLine("POWER_RULE:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}