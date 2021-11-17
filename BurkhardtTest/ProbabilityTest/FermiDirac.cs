using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
{
    private static void fermi_dirac_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    FERMI_DIRAC_SAMPLE_TEST tests FERMI_DIRAC_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        int i;
        double mean;
        int sample_num = 10000;
        int seed;
        int test;
        double u;
        double[] u_test =  {
                1.0, 2.0, 4.0, 8.0, 16.0,
                32.0, 1.0
            }
            ;
        double v;
        double[] v_test =  {
                1.0, 1.0, 1.0, 1.0, 1.0,
                1.0, 0.25
            }
            ;
        double variance;
        double[] x = new double[10000];
        double xmax;
        double xmin;

        Console.WriteLine("");
        Console.WriteLine("FERMI_DIRAC_SAMPLE_TEST");
        Console.WriteLine("  FERMI_DIRAC_SAMPLE samples the Fermi Dirac distribution.");

        for (test = 0; test < 7; test++)
        {
            u = u_test[test];
            v = v_test[test];
            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("  U =          " + u + "");
            Console.WriteLine("  V =          " + v + "");

            for (i = 0; i < sample_num; i++)
            {
                x[i] = FermiDirac.fermi_dirac_sample(u, v, ref seed);
            }

            mean = typeMethods.r8vec_mean(sample_num, x);
            variance = typeMethods.r8vec_variance(sample_num, x);
            xmax = typeMethods.r8vec_max(sample_num, x);
            xmin = typeMethods.r8vec_min(sample_num, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + mean + "");
            Console.WriteLine("  Sample mean =     " + mean + "");
            Console.WriteLine("  Sample variance = " + variance + "");
            Console.WriteLine("  Sample maximum =  " + xmax + "");
            Console.WriteLine("  Sample minimum =  " + xmin + "");
        }
    }
}