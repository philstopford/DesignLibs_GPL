using System;
using Burkardt;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace GenGaussLaguerreRuleTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for GEN_LAGUERRE_RULE.
            //
            //  Discussion:
            //
            //    This program computes a standard or exponentially weighted 
            //    Gauss-Laguerre quadrature rule and writes it to a file.
            //
            //    The user specifies:
            //    * the ORDER (number of points) in the rule;
            //    * ALPHA, the exponent of |X|;
            //    * A, the left endpoint;
            //    * B, the scale factor in the exponential;
            //    * FILENAME, the root name of the output files.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 February 2010
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
            Console.WriteLine("GEN_LAGUERRE_RULE");
            Console.WriteLine("");
            Console.WriteLine("  Compute a generalized Gauss-Laguerre rule for approximating");
            Console.WriteLine("    Integral ( a <= x < +oo ) |x-a|^ALPHA exp(-B*x(x-a)) f(x) dx");
            Console.WriteLine("  of order ORDER.");
            Console.WriteLine("");
            Console.WriteLine("  The user specifies ORDER, ALPHA, A, B, and FILENAME.");
            Console.WriteLine("");
            Console.WriteLine("  ORDER is the number of points in the rule:");
            Console.WriteLine("  ALPHA is the exponent of |X|:");
            Console.WriteLine("  A is the left endpoint (typically 0).");
            Console.WriteLine("  B is the exponential scale factor, typically 1:");
            Console.WriteLine("  FILENAME is used  to generate 3 files:");
            Console.WriteLine("  * filename_w.txt - the weight file");
            Console.WriteLine("  * filename_x.txt - the abscissa file.");
            Console.WriteLine("  * filename_r.txt - the region file.");
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
                Console.WriteLine("  Enter the value of ALPHA.");
                Console.WriteLine("  ( -1.0 < ALPHA is required.)");
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
                Console.WriteLine("  Enter the left endpoint A.");
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
                Console.WriteLine("  Enter the value of B");
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
                Console.WriteLine("  Enter FILENAME the \"root name\" of the quadrature files).");
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
            Console.WriteLine("  FILENAME \"" + filename + "\".");
            //
            //  Construct the rule.
            //
            w = new double[order];
            x = new double[order];

            kind = 5;
            CGQF.cgqf(order, kind, alpha, beta, a, b, ref x, ref w);
            //
            //  Write the rule.
            //
            r = new double[2];
            r[0] = a;
            r[1] = typeMethods.r8_huge();

            QuadratureRule.rule_write(order, filename, x, w, r);

            Console.WriteLine("");
            Console.WriteLine("GEN_LAGUERRE_RULE:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}