using System;
using Burkardt.Table;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.MonteCarlo
{
    public class MonteCarlo
    {
        public static void calc(int m, int n, int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for MONTE_CARLO_RULE.
        //
        //  Discussion:
        //
        //    MONTE_CARLO_RULE generates N points in the M-dimensional unit hypercube,
        //    and writes out files so that the data can be regarded as a quadrature rule.
        //
        //  Usage:
        //
        //    monte_carlo_rule m n seed
        //
        //    where
        //
        //    * M, the spatial dimension,
        //    * N, the number of points to generate,
        //    * SEED, the seed, a positive integer.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 March 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            //
            //  Input summary.
            //
            Console.WriteLine();
            Console.WriteLine("  M = " + m);
            Console.WriteLine("  N = " + n);
            Console.WriteLine("  S = " + seed);
            //
            //  Construct the rule.
            //
            double[] r = new double[m*2];

            for (int i = 0; i < m; i++ )
            {
                r[i+0*m] = 0.0;
                r[i+1*m] = 1.0;
            }

            double[] w = new double[n];

            for (int i = 0; i < n; i++ )
            {
                w[i] = 1.0 / n;
            }

            double[] x = UniformRNG.r8mat_uniform_01 ( m, n, ref seed );
            //
            //  Output the rule.
            //
            string filename_r = "mc_d" + m + "_n" + n + "_s" + seed + "_r.txt";

            string filename_w = "mc_d" + m + "_n" + n + "_s" + seed + "_w.txt";

            string filename_x = "mc_d" + m + "_n" + n + "_s" + seed + "_x.txt";

            Console.WriteLine();
            Console.WriteLine("  Region file will be   " + filename_r);
            Console.WriteLine("  Weight file will be   " + filename_w);
            Console.WriteLine("  Abscissa file will be " + filename_x);

            typeMethods.r8mat_write ( filename_r, m, 2, r );
            typeMethods.r8mat_write ( filename_w, 1, n, w );
            typeMethods.r8mat_write ( filename_x, m, n, x );

            Console.WriteLine();
            Console.WriteLine("MONTE_CARLO_RULE:");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine();
        }
    }
}