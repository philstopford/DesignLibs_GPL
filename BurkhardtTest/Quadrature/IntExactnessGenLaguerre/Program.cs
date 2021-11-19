using System;
using Burkardt.Laguerre;
using Burkardt.Table;
using Burkardt.Types;

namespace IntExactnessGenLaguerre;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for INT_EXACTNESS_GEN_LAGUERRE.
        //
        //  Discussion:
        //
        //    This program investigates a generalized Gauss-Laguerre quadrature rule
        //    by using it to integrate monomials over [0,+oo), and comparing the
        //    approximate result to the known exact value.
        //
        //    The user specifies:
        //    * the "root" name of the R, W and X files that specify the rule;
        //    * DEGREE_MAX, the maximum monomial degree to be checked.
        //    * ALPHA, the power of X in the weighting function.
        //    * OPTION, whether the rule is for x^alpha*exp(-x)*f(x) or f(x).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double alpha;
        int degree;
        int degree_max;
        int dim_num;
        int dim_num2;
        int i;
        int option;
        int order;
        int point_num;
        int point_num2;
        double quad_error;
        string quad_filename;
        string quad_r_filename;
        string quad_w_filename;
        string quad_x_filename;
        double[] r;
        double volume;
        double[] w;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("INT_EXACTNESS_GEN_LAGUERRE");
        Console.WriteLine("");
        Console.WriteLine("  Investigate the polynomial exactness of a generalized Gauss-Laguerre");
        Console.WriteLine("  quadrature rule by integrating exponentially weighted");
        Console.WriteLine("  monomials up to a given degree over the [0,+oo) interval.");
        Console.WriteLine("");
        Console.WriteLine("  The rule may be defined on another interval, [A,+oo)");
        Console.WriteLine("  in which case it is adjusted to the [0,+oo) interval.");
        //
        //  Get the quadrature file rootname.
        //
        try
        {
            quad_filename = args[0];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the quadrature file rootname:");
            quad_filename = Console.ReadLine();
        }

        Console.WriteLine("");
        Console.WriteLine("  The quadrature file rootname is \"" + quad_filename + "\".");
        //
        //  Create the names of:
        //    the quadrature X file;
        //    the quadrature W file;
        //    the quadrature R file;
        //
        quad_w_filename = quad_filename + "_w.txt";
        quad_x_filename = quad_filename + "_x.txt";
        quad_r_filename = quad_filename + "_r.txt";
        //
        //  Get the maximum degree:
        //
        try
        {
            degree_max = Convert.ToInt32(args[1]);

        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter DEGREE_MAX, the maximum monomial degree to check.");
            degree_max = Convert.ToInt32(Console.ReadLine());
        }

        Console.WriteLine("");
        Console.WriteLine("  The requested maximum monomial degree is = " + degree_max + "");
        //
        //  Get the exponent ALPHA:
        //
        try
        {
            alpha = Convert.ToDouble(args[2]);

        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  ALPHA is the power of X in the weighting function.");
            Console.WriteLine("");
            Console.WriteLine("  ALPHA is a real number greater than -1.0;");
            Console.WriteLine("");
            Console.WriteLine("  Enter ALPHA.");
            alpha = Convert.ToDouble(Console.ReadLine());
        }

        Console.WriteLine("");
        Console.WriteLine("  The requested value of ALPHA = " + alpha + "");
        //
        //  The fourth command line argument is OPTION.
        //  0 for the standard rule for integrating x^alpha*exp(-x)*f(x),
        //  1 for a rule for integrating f(x).
        //
        try
        {
            option = Convert.ToInt32(args[3]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("OPTION chooses the standard or modified rule.");
            Console.WriteLine("0: standard rule for integrating x^allpha*exp(-x)*f(x)");
            Console.WriteLine("1: modified rule for integrating                  f(x)");
            option = Convert.ToInt32(Console.ReadLine());
        }

        //
        //  Summarize the input.
        //
        Console.WriteLine("");
        Console.WriteLine("INT_EXACTNESS_GEN_LAGUERRE: User input:");
        Console.WriteLine("  Quadrature rule X file = \"" + quad_x_filename + "\".");
        Console.WriteLine("  Quadrature rule W file = \"" + quad_w_filename + "\".");
        Console.WriteLine("  Quadrature rule R file = \"" + quad_r_filename + "\".");
        Console.WriteLine("  Maximum degree to check = " + degree_max + "");
        Console.WriteLine("  Weighting exponent ALPHA = " + alpha + "");
        switch (option)
        {
            case 0:
                Console.WriteLine("  OPTION = 0, integrate x^alpha*exp(-x)*f(x)");
                break;
            default:
                Console.WriteLine("  OPTION = 1, integrate                 f(x)");
                break;
        }

        //
        //  Read the X file.
        //
        TableHeader h = typeMethods.r8mat_header_read(quad_x_filename);
        dim_num = h.m;
        order = h.n;

        if (dim_num != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("INT_EXACTNESS_GEN_LAGUERRE - Fatal error!");
            Console.WriteLine("  The spatial dimension of X should be 1.");
            Console.WriteLine(" The implicit input dimension was DIM_NUM = " + dim_num + "");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points  = " + order + "");

        x = typeMethods.r8mat_data_read(quad_x_filename, dim_num, order);
        //
        //  Read the W file.
        //
        h = typeMethods.r8mat_header_read(quad_w_filename);
        dim_num2 = h.m;
        point_num = h.n;

        if (dim_num2 != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("INT_EXACTNESS_GEN_LAGUERRE - Fatal error!");
            Console.WriteLine("  The quadrature weight file should have exactly");
            Console.WriteLine("  one value on each line.");
            return;
        }

        if (point_num != order)
        {
            Console.WriteLine("");
            Console.WriteLine("INT_EXACTNESS_GEN_LAGUERRE - Fatal error!");
            Console.WriteLine("  The quadrature weight file should have exactly");
            Console.WriteLine("  the same number of lines as the abscissa file.");
            return;
        }

        w = typeMethods.r8mat_data_read(quad_w_filename, dim_num, order);
        //
        //  Read the R file.
        //
        h = typeMethods.r8mat_header_read(quad_r_filename);
        dim_num2 = h.m;
        point_num2 = h.n;

        if (dim_num2 != dim_num)
        {
            Console.WriteLine("");
            Console.WriteLine("INT_EXACTNESS_GEN_LAGUERRE - Fatal error!");
            Console.WriteLine("  The quadrature region file should have the");
            Console.WriteLine("  same number of values on each line as the");
            Console.WriteLine("  abscissa file does.");
            return;
        }

        if (point_num2 != 2)
        {
            Console.WriteLine("");
            Console.WriteLine("INT_EXACTNESS_GEN_LAGUERDRE - Fatal error!");
            Console.WriteLine("  The quadrature region file should have two lines.");
            return;
        }

        r = typeMethods.r8mat_data_read(quad_r_filename, dim_num, point_num2);
        //
        //  Print the input quadrature rule.
        //
        a = r[0];

        Console.WriteLine("");
        Console.WriteLine("  The quadrature rule to be tested is");
        Console.WriteLine("  a generalized Gauss-Laguerre rule");
        Console.WriteLine("  ORDER = " + order + "");
        Console.WriteLine("  with A =     " + a + "");
        Console.WriteLine("  and  ALPHA = " + alpha + "");
        Console.WriteLine("");
        switch (option)
        {
            case 0:
                Console.WriteLine("  Standard rule:");
                Console.WriteLine("    Integral ( A <= x < +oo ) x^alpha exp(-x) f(x) dx");
                Console.WriteLine("    is to be approximated by");
                Console.WriteLine("    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).");
                break;
            default:
                Console.WriteLine("  Modified rule:");
                Console.WriteLine("    Integral ( A <= x < +oo )         f(x) dx");
                Console.WriteLine("    is to be approximated by");
                Console.WriteLine("    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).");
                break;
        }

        Console.WriteLine("");
        Console.WriteLine("  Weights W:");
        Console.WriteLine("");
        for (i = 0; i < order; i++)
        {
            Console.WriteLine("  w[" + i.ToString().PadLeft(2)
                                     + "] = " + w[i].ToString("0.################").PadLeft(24) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Abscissas X:");
        Console.WriteLine("");
        for (i = 0; i < order; i++)
        {
            Console.WriteLine("  x[" + i.ToString().PadLeft(2)
                                     + "] = " + x[i].ToString("0.################").PadLeft(24) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Region R:");
        Console.WriteLine("");

        for (i = 0; i < 2; i++)
        {
            Console.WriteLine("  r[" + i.ToString().PadLeft(2)
                                     + "] = " + r[i].ToString("0.################").PadLeft(24) + "");
        }

        //
        //  Supposing the input rule is defined on [A,+oo),
        //  rescale the weights, and translate the abscissas,
        //  so our rule is defined on [0,+oo).
        //
        volume = Math.Exp(-a);
        for (i = 0; i < order; i++)
        {
            w[i] /= volume;
        }

        for (i = 0; i < order; i++)
        {
            x[i] -= a;
        }

        //
        //  Explore the monomials.
        //
        Console.WriteLine("");
        Console.WriteLine("  A generalized Gauss-Laguerre rule would be able to exactly");
        Console.WriteLine("  integrate monomials up to and including degree = " +
                          (2 * order - 1) + "");
        Console.WriteLine("");
        Console.WriteLine("          Error          Degree");
        Console.WriteLine("");

        for (degree = 0; degree <= degree_max; degree++)
        {
            quad_error = QuadratureRule.monomial_quadrature_gen_laguerre(degree, alpha, order,
                option, w, x);

            Console.WriteLine("  " + quad_error.ToString("0.################").PadLeft(24)
                                   + "  " + degree.ToString().PadLeft(2) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("INT_EXACTNESS_GEN_LAGUERRE:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}