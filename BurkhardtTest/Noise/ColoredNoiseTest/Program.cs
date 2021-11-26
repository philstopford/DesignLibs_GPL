using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using Burkardt.Noise;
using Burkardt.Types;

namespace ColoredNoiseTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for COLORED_NOISE_TEST.
        //
        //  Discussion:
        //
        //    COLORED_NOISE_TEST tests the COLORED_NOISE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("COLORED_NOISE_TEST");
        Console.WriteLine("  Test the COLORED_NOISE library.");

        const int n = 128;
        const double q_d = 1.0;
        const int seed_init = 123456789;

        for (i = 0; i <= 8; i++)
        {
            double alpha = 0.25 * i;
            test01(n, q_d, alpha, seed_init);
        }

        Console.WriteLine("");
        Console.WriteLine("COLORED_NOISE_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(int n, double q_d, double alpha, int seed_init)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 calls F_ALPHA with particular parameters.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements of the sequence to generate.
        //
        //    Input, double Q_D, the variance of the sequence.
        //
        //    Input, double ALPHA, the exponent of the power law.
        //
        //    Input, int SEED_INIT, the initial seed for the random number generator.
        //
    {
        int i;
        List<string> output_unit = new();
        typeMethods.r8NormalData data = new();

        string output_filename = "alpha_" + alpha.ToString("0000.##") + ".txt";
        //
        //  Report parameters.
        //
        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  Generating " + n + " sample points.");
        Console.WriteLine("  1/F^ALPHA noise has ALPHA = " + alpha + "");
        Console.WriteLine("  Variance is " + q_d + "");
        Console.WriteLine("  Initial random number seed = " + seed_init + "");

        int seed = seed_init;

        double[] x = Colored.f_alpha(n, q_d, alpha, ref data, ref seed);
        //
        //  Print no more than 10 entries of the data.
        //
        typeMethods.r8vec_print_part(n, x, 10, "  Noise sample:");
        //
        //  Write the data to a file.
        //

        for (i = 0; i < n; i++)
        {
            output_unit.Add(x[i].ToString(CultureInfo.InvariantCulture));
        }

        File.WriteAllLines(output_filename, output_unit);

        Console.WriteLine("  Data written to file \"" + output_filename + "\"");
    }
}