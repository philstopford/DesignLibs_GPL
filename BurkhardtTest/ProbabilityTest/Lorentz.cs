using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
{
    private static void lorentz_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    LORENTZ_CDF_TEST tests LORENTZ_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 April 2016
//
//  Author:
//
//    John Burkardt
//

    {
        double cdf;
        int i;
        double pdf;
        int seed = 123456789;
        double x;
        double x2;

        Console.WriteLine("");
        Console.WriteLine("LORENTZ_CDF_TEST");
        Console.WriteLine("  LORENTZ_CDF evaluates the Lorentz CDF;");
        Console.WriteLine("  LORENTZ_CDF_INV inverts the Lorentz CDF.");
        Console.WriteLine("  LORENTZ_PDF evaluates the Lorentz PDF;");

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            x = Lorentz.lorentz_sample(ref seed);
            pdf = Lorentz.lorentz_pdf(x);
            cdf = Lorentz.lorentz_cdf(x);
            x2 = Lorentz.lorentz_cdf_inv(cdf);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    private static void lorentz_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    LORENTZ_SAMPLE_TEST tests LORENTZ_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        int SAMPLE_NUM = 1000;

        int i;
        double mean;
        int seed = 123456789;
        double variance;
        double[] x = new double [SAMPLE_NUM];
        double xmax;
        double xmin;

        Console.WriteLine("");
        Console.WriteLine("LORENTZ_SAMPLE_TEST");
        Console.WriteLine("  LORENTZ_MEAN computes the Lorentz mean;");
        Console.WriteLine("  LORENTZ_SAMPLE samples the Lorentz distribution;");
        Console.WriteLine("  LORENTZ_VARIANCE computes the Lorentz variance.");

        mean = Lorentz.lorentz_mean();
        variance = Lorentz.lorentz_variance();

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Lorentz.lorentz_sample(ref seed);
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