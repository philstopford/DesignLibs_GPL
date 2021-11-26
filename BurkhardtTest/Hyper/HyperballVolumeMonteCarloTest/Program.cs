using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace HyperballVolumeMonteCarloTest;

using Montecarlo = Burkardt.HyperGeometry.Hyperball.MonteCarlo;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for HYPERBALL_VOLUME_MONTE_CARLO.
        //
        //  Discussion:
        //
        //    DIM_NUM = 6 is a reasonable test.
        //
        //    N_LOG2_MAX = 25 puts a strain on the system, since we generate that
        //    many temporary points at once.  To solve bigger problems, it would
        //    be better to compute the new points in batches whose maximum size
        //    is limited.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim_num;
        double exact = 0;
        int n = 0;
        int n_more = 0;
        int n_log2;
        const int n_log2_max = 25;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("HYPERBALL_VOLUME_MONTE_CARLO:");
        Console.WriteLine("  Use a Monte Carlo approach to estimate the volume of");
        Console.WriteLine("  the unit hyperball in M dimensions.");
        //
        //  Get the quadrature file root name:
        //
        try
        {
            dim_num = Convert.ToInt32(args[0]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("HYPERBALL_VOLUME_MONTE_CARLO:");
            Console.WriteLine("  Enter the spatial dimension:");

            dim_num = Convert.ToInt32(Console.ReadLine());
        }

        //
        //  Get the random number seed, if supplied.
        //
        try
        {
            seed = Convert.ToInt32(args[1]);
        }
        catch
        {
            seed = 123456789;
            Console.WriteLine("");
            Console.WriteLine("HYPERBALL_VOLUME_MONTE_CARLO:");
            Console.WriteLine("  Using default seed for random number generator.");
        }

        //
        //  Report user input.
        //
        Console.WriteLine("");
        Console.WriteLine("  The spatial dimension is  " + dim_num + "");
        Console.WriteLine("  The random number seed is " + seed + "");
        //
        //  Begin computation.
        //
        Console.WriteLine("");
        Console.WriteLine("    Log(N)         N      Estimate         Error");
        Console.WriteLine("");

        double quad = 0.0;
        double volume = Math.Pow(2.0, dim_num);

        for (n_log2 = 0; n_log2 <= n_log2_max; n_log2++)
        {
            switch (n_log2)
            {
                case 0:
                    quad = 0.0;
                    n_more = 1;
                    n = 0;
                    break;
                case 1:
                    n_more = 1;
                    break;
                default:
                    n_more = 2 * n_more;
                    break;
            }

            double[] x = UniformRNG.r8mat_uniform_01_new(dim_num, n_more, ref seed);
            //
            //  Rescale X from [0,1] to [-1,1].
            //
            int j;
            for (j = 0; j < n_more; j++)
            {
                int i;
                for (i = 0; i < dim_num; i++)
                {
                    x[i + j * dim_num] = 2.0 * x[i + j * dim_num] - 1.0;
                }
            }

            double[] fx = Montecarlo.hyperball01_indicator(dim_num, n_more, x);

            double quad_more = typeMethods.r8vec_sum(n_more, fx);

            //
            //  Incorporate the new data into the totals.
            //
            n += n_more;
            quad += quad_more;

            double estimate = volume * quad / n;
            exact = Montecarlo.hyperball01_volume(dim_num);
            double error = Math.Abs(exact - estimate);
            Console.WriteLine("  " + n_log2.ToString().PadLeft(8)
                                   + "  " + n.ToString().PadLeft(8)
                                   + "  " + estimate.ToString("0.##########").PadLeft(16)
                                   + "  " + error.ToString("0.##").PadLeft(16) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("        oo        oo"
                          + "  " + exact.ToString("0.##########").PadLeft(16)
                          + "  " + 0.0.ToString("0.##").PadLeft(10) + "");

        Console.WriteLine("");
        Console.WriteLine("HYPERBALL_VOLUME_MONTE_CARLO:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}