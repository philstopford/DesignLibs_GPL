using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
{
    private static void anglit_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    ANGLIT_CDF_TEST tests ANGLIT_CDF, ANGLIT_CDF_INV, ANGLIT_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
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
        Console.WriteLine("ANGLIT_CDF_TEST");
        Console.WriteLine("  ANGLIT_CDF evaluates the Anglit CDF;");
        Console.WriteLine("  ANGLIT_CDF_INV inverts the Anglit CDF.");
        Console.WriteLine("  ANGLIT_PDF evaluates the Anglit PDF;");

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            x = Anglit.anglit_sample(ref seed);
            pdf = Anglit.anglit_pdf(x);
            cdf = Anglit.anglit_cdf(x);
            x2 = Anglit.anglit_cdf_inv(cdf);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void anglit_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    ANGLIT_SAMPLE_TEST tests ANGLIT_MEAN, ANGLIT_SAMPLE, ANGLIT_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
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
        double[] x = new double[SAMPLE_NUM];
        double xmax;
        double xmin;

        Console.WriteLine("");
        Console.WriteLine("ANGLIT_SAMPLE_TEST");
        Console.WriteLine("  ANGLIT_MEAN computes the Anglit mean;");
        Console.WriteLine("  ANGLIT_SAMPLE samples the Anglit distribution;");
        Console.WriteLine("  ANGLIT_VARIANCE computes the Anglit variance.");

        mean = Anglit.anglit_mean();
        variance = Anglit.anglit_variance();

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Anglit.anglit_sample(ref seed);
        }

        mean = typeMethods.r8vec_mean(i, x);
        variance = typeMethods.r8vec_variance(i, x);
        xmax = typeMethods.r8vec_max(i, x);
        xmin = typeMethods.r8vec_min(i, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
        Console.WriteLine("  Sample mean =     " + mean + "");
        Console.WriteLine("  Sample variance = " + variance + "");
        Console.WriteLine("  Sample maximum =  " + xmax + "");
        Console.WriteLine("  Sample minimum =  " + xmin + "");

    }
        
}