using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest
{
    partial class Program
    {
        static void disk_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    DISK_SAMPLE_TEST tests DISK_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 March 2016
//
//  Author:
//
//    John Burkardt
//
        {
            int SAMPLE_NUM = 1000;

            double a;
            double b;
            double c;
            int j;
            double[] mean;
            int seed = 123456789;
            double variance;
            double[] x = new double [2 * SAMPLE_NUM];
            double[] xmax;
            double[] xmin;
            double[] y;

            Console.WriteLine("");
            Console.WriteLine("DISK_SAMPLE");
            Console.WriteLine("  DISK_MEAN returns the Disk mean;");
            Console.WriteLine("  DISK_SAMPLE samples the Disk distribution;");
            Console.WriteLine("  DISK_VARIANCE returns the Disk variance;");

            a = 10.0;
            b = 4.0;
            c = 3.0;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =      " + a + "");
            Console.WriteLine("  PDF parameter B =      " + b + "");
            Console.WriteLine("  PDF parameter C =      " + c + "");

            mean = Disk.disk_mean(a, b, c);
            variance = Disk.disk_variance(a, b, c);

            Console.WriteLine("");
            Console.WriteLine("  Disk mean ="
                              + "  " + mean[0].ToString().PadLeft(12)
                              + "  " + mean[1].ToString().PadLeft(12) + "");
            Console.WriteLine("  Disk variance = " + variance.ToString().PadLeft(12) + "");


            for (j = 0; j < SAMPLE_NUM; j++)
            {
                y = Disk.disk_sample(a, b, c, ref seed);
                x[0 + j * 2] = y[0];
                x[1 + j * 2] = y[1];
            }

            variance = 0.0;
            for (j = 0; j < SAMPLE_NUM; j++)
            {
                variance = variance + Math.Pow(x[0 + j * 2] - a, 2) + Math.Pow(x[1 + j * 2] - b, 2);
            }

            variance = variance / (double) (SAMPLE_NUM);


            mean = typeMethods.r8row_mean(2, SAMPLE_NUM, x);
            xmax = typeMethods.r8row_max(2, SAMPLE_NUM, x);
            xmin = typeMethods.r8row_min(2, SAMPLE_NUM, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
            Console.WriteLine("  Sample mean =     "
                              + mean[0].ToString().PadLeft(12) + "  "
                              + mean[1].ToString().PadLeft(12) + "");
            Console.WriteLine("  Sample variance = "
                              + variance.ToString().PadLeft(12) + "");
            Console.WriteLine("  Sample maximum =  "
                              + xmax[0].ToString().PadLeft(12) + "  "
                              + xmax[1].ToString().PadLeft(12) + "");
            Console.WriteLine("  Sample minimum =  "
                              + xmin[0].ToString().PadLeft(12) + "  "
                              + xmin[1].ToString().PadLeft(12) + "");

        }
    }
}