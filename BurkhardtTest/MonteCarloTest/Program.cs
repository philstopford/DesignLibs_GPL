using System;
using Burkardt.MonteCarlo;

namespace MonteCarloTest;

internal class Program
{
    private static void Main(string[] args)
    {
        Console.WriteLine();
        Console.WriteLine("MONTE_CARLO_RULE");
        Console.WriteLine("");
        Console.WriteLine("  Compute the abscissas and weights of a quadrature rule");
        Console.WriteLine("  that is simply a Monte Carlo sampling.");
        Console.WriteLine("");
        Console.WriteLine("  The program requests input values from the user:");
        Console.WriteLine("");
        Console.WriteLine("  * M, the spatial dimension,");
        Console.WriteLine("  * N, the number of points to generate,");
        Console.WriteLine("  * SEED, a positive integer.");
        Console.WriteLine("");
        Console.WriteLine("  Output from the program includes");
        Console.WriteLine("  a set of 3 files that define the quadrature rule.");
        Console.WriteLine("");
        Console.WriteLine("    (1) \"mc_m?_n?_s?_r.txt\", the ranges;");
        Console.WriteLine("    (2) \"mc_m?_n?_s?_w.txt\", the weights;");
        Console.WriteLine("    (3) \"mc_m?_n?_s?_x.txt\", the abscissas.");

        args = args.Length switch
        {
            0 => new[] {"3", "2", "12345"},
            _ => args
        };

        test01(args);
    }

    private static void test01(string[] args)
    {
        //
        //  Get the spatial dimension M.
        //
        int m, n, seed;
        try
        {
            m = Convert.ToInt32(args[0]);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the spatial dimension M (1 or greater)");
            m = Convert.ToInt32(Console.Read());
        }

        m = Math.Max(1, m);

        //
        //  Get the number of points N.
        //
        try
        {
            n = Convert.ToInt32(args[1]);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the number of points N (1 or greater):");
            n = Convert.ToInt32(Console.Read());
        }

        n = Math.Max(1, n);

        //
        //  Get the seed S.
        //
        try
        {
            seed = Convert.ToInt32(args[2]);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the seed S (1 or greater):");
            seed = Convert.ToInt32(Console.Read());
        }
        seed = Math.Max(1, seed);

        MonteCarlo.calc(m, n, seed);
    }
}