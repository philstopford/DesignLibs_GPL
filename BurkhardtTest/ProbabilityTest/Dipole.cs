using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
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
        int TEST_NUM = 3;

        double a;
        double b;
        double cdf;
        int i;
        double pdf;
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
        double x;
        double x2;

        Console.WriteLine("");
        Console.WriteLine("DIPOLE_CDF_TEST");
        Console.WriteLine("  DIPOLE_CDF evaluates the Dipole CDF;");
        Console.WriteLine("  DIPOLE_CDF_INV inverts the Dipole CDF.");
        Console.WriteLine("  DIPOLE_PDF evaluates the Dipole PDF;");

        for (test_i = 0; test_i < TEST_NUM; test_i++)
        {
            a = test_a[test_i];
            b = test_b[test_i];

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

            for (i = 1; i <= 10; i++)
            {
                x = Dipole.dipole_sample(a, b, ref seed);
                pdf = Dipole.dipole_pdf(x, a, b);
                cdf = Dipole.dipole_cdf(x, a, b);
                x2 = Dipole.dipole_cdf_inv(cdf, a, b);

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + pdf.ToString().PadLeft(12) + "  "
                                  + cdf.ToString().PadLeft(12) + "  "
                                  + x2.ToString().PadLeft(12) + "");
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
        int SAMPLE_NUM = 1000;
        int TEST_NUM = 3;

        double a;
        double b;
        int i;
        double mean;
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
        double variance;
        double[] x = new double[SAMPLE_NUM];
        double xmax;
        double xmin;

        Console.WriteLine("");
        Console.WriteLine("DIPOLE_SAMPLE_TEST");
        Console.WriteLine("  DIPOLE_SAMPLE samples the Dipole distribution;");

        for (test_i = 0; test_i < TEST_NUM; test_i++)
        {
            a = test_a[test_i];
            b = test_b[test_i];

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

            for (i = 0; i < SAMPLE_NUM; i++)
            {
                x[i] = Dipole.dipole_sample(a, b, ref seed);
            }

            mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
            variance = typeMethods.r8vec_variance(SAMPLE_NUM, x);
            xmax = typeMethods.r8vec_max(SAMPLE_NUM, x);
            xmin = typeMethods.r8vec_min(SAMPLE_NUM, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
            Console.WriteLine("  Sample mean =     " + mean + "");
            Console.WriteLine("  Sample variance = " + variance + "");
            Console.WriteLine("  Sample maximum =  " + xmax + "");
            Console.WriteLine("  Sample minimum =  " + xmin + "");
        }
    }

}