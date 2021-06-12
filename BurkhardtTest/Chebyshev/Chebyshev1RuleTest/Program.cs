using System;
using Burkardt;

namespace Chebyshev1RuleTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for CHEBYSHEV1_RULE.
            //
            //  Discussion:
            //
            //    This program computes a standard Gauss-Chebyshev type 1 quadrature rule
            //    and writes it to a file.
            //
            //    The user specifies:
            //    * the ORDER (number of points) in the rule
            //    * A, the left endpoint;
            //    * B, the right endpoint;
            //    * FILENAME, the root name of the output files.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 February 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            double alpha;
            double b;
            double beta;
            string filename;
            int kind;
            int order;
            double[] r;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("CHEBYSHEV1_RULE");
            Console.WriteLine("");
            Console.WriteLine("  Compute a Gauss-Chebyshev type 1 rule for approximating");
            Console.WriteLine("");
            Console.WriteLine("    Integral ( A <= x <= B ) f(x) / sqrt ( ( x - A ) * ( B - x ) ) dx");
            Console.WriteLine("");
            Console.WriteLine("  of order ORDER.");
            Console.WriteLine("");
            Console.WriteLine("  The user specifies ORDER, A, B and FILENAME.");
            Console.WriteLine("");
            Console.WriteLine("  Order is the number of points.");
            Console.WriteLine("");
            Console.WriteLine("  A is the left endpoint.");
            Console.WriteLine("");
            Console.WriteLine("  B is the right endpoint.");
            Console.WriteLine("");
            Console.WriteLine("  FILENAME is used to generate 3 files:");
            Console.WriteLine("");
            Console.WriteLine("    filename_w.txt - the weight file");
            Console.WriteLine("    filename_x.txt - the abscissa file.");
            Console.WriteLine("    filename_r.txt - the region file.");
            //
            //  Initialize parameters;
            //
            alpha = 0.0;
            beta = 0.0;
            //
            //  Get ORDER.
            //
            try
            {
                order = Convert.ToInt32(args[0]);
            }
            catch (Exception)
            {
                Console.WriteLine("");
                Console.WriteLine("  Enter the value of ORDER (1 or greater)");
                order = Convert.ToInt32(Console.ReadLine());
            }

            //
            //  Get A.
            //
            try
            {
                a = Convert.ToDouble(args[1]);
            }
            catch (Exception)
            {
                Console.WriteLine("");
                Console.WriteLine("  Enter the value of A:");
                a = Convert.ToDouble(Console.ReadLine());
            }

            //
            //  Get B.
            //
            try
            {
                b = Convert.ToDouble(args[2]);
            }
            catch (Exception)
            {
                Console.WriteLine("");
                Console.WriteLine("  Enter the value of B:");
                b = Convert.ToDouble(Console.ReadLine());
            }

            //
            //  Get FILENAME:
            //
            try
            {
                filename = args[3];
            }
            catch (Exception)
            {
                Console.WriteLine("");
                Console.WriteLine("  Enter FILENAME, the \"root name\" of the quadrature files).");
                filename = Console.ReadLine();
            }

            //
            //  Input summary.
            //
            Console.WriteLine("");
            Console.WriteLine("  ORDER = " + order + "");
            Console.WriteLine("  A = " + a + "");
            Console.WriteLine("  B = " + b + "");
            Console.WriteLine("  FILENAME = \"" + filename + "\".");
            //
            //  Construct the rule.
            //
            w = new double[order];
            x = new double[order];

            kind = 2;
            CQGF.cgqf(order, kind, alpha, beta, a, b, x, w);
            //
            //  Write the rule.
            //
            r = new double[2];
            r[0] = a;
            r[1] = b;

            QuadratureRule.rule_write(order, filename, x, w, r);

            Console.WriteLine("");
            Console.WriteLine("CHEBYSHEV1_RULE:");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }
    }
}