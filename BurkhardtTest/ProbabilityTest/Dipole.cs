using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void dipole_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    DIPOLE_CDF_TEST tests DIPOLE_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
    {
        const int TEST_NUM = 3;

        int seed = 123456789;
        int test_i;
        double[] test_a =  {
                0.0, Math.PI / 4.0, Math.PI / 2.0
            }
            ;
        double[] test_b =  {
                1.0, 0.5, 0.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("DIPOLE_CDF_TEST");
        Console.WriteLine("  DIPOLE_CDF evaluates the Dipole CDF;");
        Console.WriteLine("  DIPOLE_CDF_INV inverts the Dipole CDF.");
        Console.WriteLine("  DIPOLE_PDF evaluates the Dipole PDF;");

        for (test_i = 0; test_i < TEST_NUM; test_i++)
        {
            double a = test_a[test_i];
            double b = test_b[test_i];

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =      " + a + "");
            Console.WriteLine("  PDF parameter B =      " + b + "");

            if (!Dipole.dipole_check(a, b))
            {
                Console.WriteLine("");
                Console.WriteLine("DIPOLE_CDF_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("       X            PDF           CDF            CDF_INV");
            Console.WriteLine("");

            int i;
            for (i = 1; i <= 10; i++)
            {
                double x = Dipole.dipole_sample(a, b, ref seed);
                double pdf = Dipole.dipole_pdf(x, a, b);
                double cdf = Dipole.dipole_cdf(x, a, b);
                double x2 = Dipole.dipole_cdf_inv(cdf, a, b);

                Console.WriteLine("  "
                                  + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }

        }
    }

    private static void dipole_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    DIPOLE_SAMPLE_TEST tests DIPOLE_SAMPLE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
    {
        const int SAMPLE_NUM = 1000;
        const int TEST_NUM = 3;

        int seed = 123456789;
        double[] test_a =  {
                0.0, Math.PI / 4.0, Math.PI / 2.0
            }
            ;
        double[] test_b =  {
                1.0, 0.5, 0.0
            }
            ;
        int test_i;
        double[] x = new double[SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("DIPOLE_SAMPLE_TEST");
        Console.WriteLine("  DIPOLE_SAMPLE samples the Dipole distribution;");

        for (test_i = 0; test_i < TEST_NUM; test_i++)
        {
            double a = test_a[test_i];
            double b = test_b[test_i];

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =      " + a + "");
            Console.WriteLine("  PDF parameter B =      " + b + "");

            if (!Dipole.dipole_check(a, b))
            {
                Console.WriteLine("");
                Console.WriteLine("DIPOLE_SAMPLE_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            int i;
            for (i = 0; i < SAMPLE_NUM; i++)
            {
                x[i] = Dipole.dipole_sample(a, b, ref seed);
            }

            double mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
            double variance = typeMethods.r8vec_variance(SAMPLE_NUM, x);
            double xmax = typeMethods.r8vec_max(SAMPLE_NUM, x);
            double xmin = typeMethods.r8vec_min(SAMPLE_NUM, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
            Console.WriteLine("  Sample mean =     " + mean + "");
            Console.WriteLine("  Sample variance = " + variance + "");
            Console.WriteLine("  Sample maximum =  " + xmax + "");
            Console.WriteLine("  Sample minimum =  " + xmin + "");
        }
    }

}