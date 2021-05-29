using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace Burkardt.ProbabilityTest
{
    partial class Program
    {
        static void circular_normal_01_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    Circular.circular_NORMAL_01_SAMPLE_TEST tests Circular.circular_NORMAL_01_SAMPLE.
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

            int j;
            double[] mean;
            int seed = 123456789;
            double[] variance;
            double[] x = new double [2 * SAMPLE_NUM];
            double[] xmax;
            double[] xmin;
            double[] y;

            Console.WriteLine("");
            Console.WriteLine("Circular.circular_NORMAL_01_SAMPLE_TEST");
            Console.WriteLine("  Circular.circular_NORMAL_01_MEAN computes the Circular Normal 01 mean;");
            Console.WriteLine("  Circular.circular_NORMAL_01_SAMPLE samples the Circular Normal 01 distribution;");
            Console.WriteLine("  Circular.circular_NORMAL_01_VARIANCE computes the Circular Normal 01 variance.");

            mean = Circular.circular_normal_01_mean();
            variance = Circular.circular_normal_01_variance();

            Console.WriteLine("");
            Console.WriteLine("  PDF mean =     "
                              + mean[0].ToString().PadLeft(12) + "  "
                              + mean[1].ToString().PadLeft(12) + "");
            Console.WriteLine("  PDF variance = "
                              + variance[0].ToString().PadLeft(12) + "  "
                              + variance[1].ToString().PadLeft(12) + "");

            for (j = 0; j < SAMPLE_NUM; j++)
            {
                y = Circular.circular_normal_01_sample(ref seed);
                x[0 + j * 2] = y[0];
                x[1 + j * 2] = y[1];
            }

            mean = typeMethods.r8row_mean(2, SAMPLE_NUM, x);
            variance = typeMethods.r8row_variance(2, SAMPLE_NUM, x);
            xmax = typeMethods.r8row_max(2, SAMPLE_NUM, x);
            xmin = typeMethods.r8row_min(2, SAMPLE_NUM, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
            Console.WriteLine("  Sample mean =     "
                              + mean[0].ToString().PadLeft(12) + "  "
                              + mean[1].ToString().PadLeft(12) + "");
            Console.WriteLine("  Sample variance = "
                              + variance[0].ToString().PadLeft(12) + "  "
                              + variance[1].ToString().PadLeft(12) + "");
            Console.WriteLine("  Sample maximum =  "
                              + xmax[0].ToString().PadLeft(12) + "  "
                              + xmax[1].ToString().PadLeft(12) + "");
            Console.WriteLine("  Sample minimum =  "
                              + xmin[0].ToString().PadLeft(12) + "  "
                              + xmin[1].ToString().PadLeft(12) + "");

        }

        static void circular_normal_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    Circular.circular_NORMAL_SAMPLE_TEST tests Circular.circular_NORMAL_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2012
//
//  Author:
//
//    John Burkardt
//
        {
            int SAMPLE_NUM = 1000;

            double[] a = new double [2];
            double b;
            int j;
            double[] mean;
            int seed = 123456789;
            double[] variance;
            double[] x = new double [2 * SAMPLE_NUM];
            double[] xmax;
            double[] xmin;
            double[] y;

            a[0] = 1.0;
            a[1] = 5.0;
            b = 0.75;

            Console.WriteLine("");
            Console.WriteLine("Circular.circular_NORMAL_SAMPLE_TEST");
            Console.WriteLine("  Circular.circular_NORMAL_MEAN computes the Circular Normal mean;");
            Console.WriteLine("  Circular.circular_NORMAL_SAMPLE samples the Circular Normal distribution;");
            Console.WriteLine("  Circular.circular_NORMAL_VARIANCE computes the Circular Normal variance.");

            mean = Circular.circular_normal_mean(a, b);
            variance = Circular.circular_normal_variance(a, b);

            Console.WriteLine("");
            Console.WriteLine("  PDF mean =     "
                              + mean[0].ToString().PadLeft(12) + "  "
                              + mean[1].ToString().PadLeft(12) + "");
            Console.WriteLine("  PDF variance = "
                              + variance[0].ToString().PadLeft(12) + "  "
                              + variance[1].ToString().PadLeft(12) + "");

            for (j = 0; j < SAMPLE_NUM; j++)
            {
                y = Circular.circular_normal_sample(a, b, ref seed);
                x[0 + j * 2] = y[0];
                x[1 + j * 2] = y[1];
            }

            mean = typeMethods.r8row_mean(2, SAMPLE_NUM, x);
            variance = typeMethods.r8row_variance(2, SAMPLE_NUM, x);
            xmax = typeMethods.r8row_max(2, SAMPLE_NUM, x);
            xmin = typeMethods.r8row_min(2, SAMPLE_NUM, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
            Console.WriteLine("  Sample mean =     "
                              + mean[0].ToString().PadLeft(12) + "  "
                              + mean[1].ToString().PadLeft(12) + "");
            Console.WriteLine("  Sample variance = "
                              + variance[0].ToString().PadLeft(12) + "  "
                              + variance[1].ToString().PadLeft(12) + "");
            Console.WriteLine("  Sample maximum =  "
                              + xmax[0].ToString().PadLeft(12) + "  "
                              + xmax[1].ToString().PadLeft(12) + "");
            Console.WriteLine("  Sample minimum =  "
                              + xmin[0].ToString().PadLeft(12) + "  "
                              + xmin[1].ToString().PadLeft(12) + "");

        }

    }
}