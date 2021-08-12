using System;
using Burkardt;
using Burkardt.Quadrature;

namespace GegenbauerQuadratureRuleTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for GEGENBAUER_RULE.
            //
            //  Discussion:
            //
            //    This program computes a standard Gauss-Gegenbauer quadrature rule
            //    and writes it to a file.
            //
            //    The user specifies:
            //    * the ORDER (number of points) in the rule
            //    * the ALPHA parameter,
            //    * A, the left endpoint.
            //    * B, the right endpoint.
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
            Console.WriteLine("");
            Console.WriteLine("GEGENBAUER_RULE");
            Console.WriteLine("");
            Console.WriteLine("  Compute a Gauss-Gegenbauer quadrature rule for approximating");
            Console.WriteLine("    Integral ( A <= x <= B ) ((x-A)(B-X))^ALPHA f(x) dx");
            Console.WriteLine("  of order ORDER.");
            Console.WriteLine("");
            Console.WriteLine("  The user specifies ORDER, ALPHA, A, B, and FILENAME.");
            Console.WriteLine("");
            Console.WriteLine("  ORDER is the number of points:");
            Console.WriteLine("  ALPHA is the exponent:");
            Console.WriteLine("  A is the left endpoint:");
            Console.WriteLine("  B is the right endpoint:");
            Console.WriteLine("  FILENAME is used to generate 3 files:");
            Console.WriteLine("    filename_w.txt - the weight file");
            Console.WriteLine("    filename_x.txt - the abscissa file.");
            Console.WriteLine("    filename_r.txt - the region file.");
            //
            //  Initialize parameters;
            //
            beta = 0.0;
            //
            //  Get ORDER.
            //
            try
            {
                order = Convert.ToInt32(args[0]);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("  Enter the value of ORDER (1 or greater)");
                order = Convert.ToInt32(Console.ReadLine());
            }

            //
            //  Get ALPHA.
            //
            try
            {
                alpha = Convert.ToDouble(args[1]);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("  ALPHA is the exponent of ((x-A)(B-X)) in the integral:");
                Console.WriteLine("  Note that -1.0 < ALPHA is required.");
                Console.WriteLine("  Enter the value of ALPHA:");
                alpha = Convert.ToDouble(Console.ReadLine());
            }

            //
            //  Get A.
            //
            try
            {
                a = Convert.ToDouble(args[2]);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("  Enter the left endpoint A:");
                a = Convert.ToDouble(Console.ReadLine());
            }

            //
            //  Get B.
            //
            try
            {
                b = Convert.ToDouble(args[3]);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("  Enter the right endpoint B:");
                b = Convert.ToDouble(Console.ReadLine());
            }

            //
            //  Get FILENAME.
            //
            try
            {
                filename = args[4];
            }
            catch
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
            Console.WriteLine("  ALPHA = " + alpha + "");
            Console.WriteLine("  A = " + a + "");
            Console.WriteLine("  B = " + b + "");
            Console.WriteLine("  FILENAME = \"" + filename + "\".");
            //
            //  Construct the rule.
            //
            w = new double[order];
            x = new double[order];

            kind = 3;
            CGQF.cgqf(order, kind, alpha, beta, a, b, ref x, ref w);
            //
            //  Write the rule.
            //
            r = new double[2];
            r[0] = a;
            r[1] = b;

            QuadratureRule.rule_write(order, filename, x, w, r);

            Console.WriteLine("");
            Console.WriteLine("GEGENBAUER_RULE:");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }
    }
}