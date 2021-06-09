using System;
using Burkardt;

namespace ClenshawCurtisTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for CCN_RULE.
            //
            //  Discussion:
            //
            //    This program computes a nested Clenshaw Curtis quadrature rule
            //    and writes it to a file.
            //
            //    The user specifies:
            //    * N, the number of points in the rule;
            //    * A, the left endpoint;
            //    * B, the right endpoint;
            //    * FILENAME, which defines the output filenames.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 March 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            double b;
            string filename;
            int n;
            double[] r;
            double[] w;
            double[] x;
            double x_max;
            double x_min;

            Console.WriteLine("");
            Console.WriteLine("CCN_RULE");
            Console.WriteLine("  Compute one of a family of nested Clenshaw Curtis rules");
            Console.WriteLine("  for approximating");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) dx");
            Console.WriteLine("  of order N.");
            Console.WriteLine("");
            Console.WriteLine("  The user specifies N, A, B and FILENAME.");
            Console.WriteLine("");
            Console.WriteLine("  N is the number of points.");
            Console.WriteLine("  A is the left endpoint.");
            Console.WriteLine("  B is the right endpoint.");
            Console.WriteLine("  FILENAME is used to generate 3 files:");
            Console.WriteLine("    filename_w.txt - the weight file");
            Console.WriteLine("    filename_x.txt - the abscissa file.");
            Console.WriteLine("    filename_r.txt - the region file.");
            //
            //  Get N.
            //
            try
            {
                n = Convert.ToInt32(args[0]);
            }
            catch (Exception)
            {
                Console.WriteLine("");
                Console.WriteLine("  Enter the value of N (1 or greater)");
                string tmp = "";
                tmp = Console.ReadLine();
                n = Convert.ToInt32(tmp);
            }

            //
            //  Get A.
            //
            try
            {
                a = Convert.ToInt32(args[1]);
            }
            catch (Exception)
            {
                Console.WriteLine("");
                Console.WriteLine("  Enter the left endpoint A:");
                string tmp = "";
                tmp = Console.ReadLine();
                a = Convert.ToInt32(tmp);
            }

            //
            //  Get B.
            //
            try
            {
                b = Convert.ToInt32(args[2]);
            }
            catch (Exception)
            {
                Console.WriteLine("");
                Console.WriteLine("  Enter the right endpoint B:");
                string tmp = "";
                tmp = Console.ReadLine();
                b = Convert.ToInt32(tmp);
            }

            //
            //  Get FILENAME:
            //
            try
            {
                filename = (args[3]);
            }
            catch (Exception)
            {
                Console.WriteLine("");
                Console.WriteLine("  Enter FILENAME, the \"root name\" of the quadrature files.");
                filename = Console.ReadLine();
            }

            //
            //  Input summary.
            //
            Console.WriteLine("");
            Console.WriteLine("  N = " + n + "");
            Console.WriteLine("  A = " + a + "");
            Console.WriteLine("  B = " + b + "");
            Console.WriteLine("  FILENAME = \"" + filename + "\".");
            //
            //  Construct the rule.
            //
            r = new double[2];

            r[0] = a;
            r[1] = b;

            x = ClenshawCurtis.ccn_compute_points_new(n);

            x_min = -1.0;
            x_max = +1.0;
            w = ClenshawCurtis.nc_compute_new(n, x_min, x_max, x);
            //
            //  Rescale the rule.
            //
            ClenshawCurtis.rescale(a, b, n, ref x, ref w);
            //
            //  Output the rule.
            //
            ClenshawCurtis.rule_write(n, filename, x, w, r);

            Console.WriteLine("");
            Console.WriteLine("CCN_RULE:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}