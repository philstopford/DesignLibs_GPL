using System;
using Burkardt;
using Burkardt.Legendre;

namespace LegendreRuleFastTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for LEGENDRE_RULE.
            //
            //  Discussion:
            //
            //    This program computes a standard Gauss-Legendre quadrature rule
            //    and writes it to a file.
            //
            //  Usage:
            //
            //    legendre_rule_fast n a b
            //
            //    where
            //
            //    * n is the number of points in the rule;
            //    * a is the left endpoint;
            //    * b is the right endpoint.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 October 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            double b;
            int n;

            Console.WriteLine("");
            Console.WriteLine("LEGENDRE_RULE");
            Console.WriteLine("");
            Console.WriteLine("  Compute a Gauss-Legendre rule for approximating");
            Console.WriteLine("");
            Console.WriteLine("    Integral ( a <= x <= b ) f(x) dx");
            Console.WriteLine("");
            Console.WriteLine("  of order N.");
            Console.WriteLine("");
            Console.WriteLine("  The computed rule is written to 3 files:");
            Console.WriteLine("");
            Console.WriteLine("  * leg_oN_w.txt - the weight file");
            Console.WriteLine("  * leg_oN_x.txt - the abscissa file.");
            Console.WriteLine("  * leg_oN_r.txt - the region file.");
            //
            //  Get N.
            //
            try
            {
                n = Convert.ToInt32(args[0]);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("  Enter N:");
                n = Convert.ToInt32(Console.ReadLine());
            }

            Console.WriteLine("");
            Console.WriteLine("  N = " + n + "");
            //
            //  Get A:
            //
            try
            {
                a = Convert.ToDouble(args[1]);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("  Enter A.");
                a = Convert.ToDouble(Console.ReadLine());
            }

            Console.WriteLine("");
            Console.WriteLine("  A = " + a + "");
            //
            //  Get B:
            //
            try
            {
                b = Convert.ToDouble(args[2]);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("  Enter B.");
                b = Convert.ToDouble(Console.ReadLine());
            }

            Console.WriteLine("");
            Console.WriteLine("  B = " + b + "");
            //
            //  Construct the rule and output it.
            //
            QuadratureRuleFast.legendre_handle(n, a, b);

            Console.WriteLine("");
            Console.WriteLine("LEGENDRE_RULE_FAST:");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }
    }
}