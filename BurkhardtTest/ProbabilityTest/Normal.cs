using System;
using System.Globalization;
using Burkardt.CDFLib;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void normal_01_cdf_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    Normal.normal_01_cdf_TEST tests Normal.normal_01_cdf.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("Normal.normal_01_cdf_TEST");
        Console.WriteLine("  Normal.normal_01_cdf evaluates the Normal 01 CDF;");
        Console.WriteLine("  Normal.normal_01_cdf_INV inverts the Normal 01 CDF.");
        Console.WriteLine("  Normal.normal_01_pdf evaluates the Normal 01 PDF;");

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Normal.normal_01_sample(ref seed);
            double pdf = Normal.normal_01_pdf(x);
            double cdf = Normal.normal_01_cdf(x);
            double x2 = Normal.normal_01_cdf_inv(cdf);

            Console.WriteLine("  "
                              + x.ToString("0.################").PadLeft(24) + "  "
                              + pdf.ToString("0.######").PadLeft(12) + "  "
                              + cdf.ToString("0.####").PadLeft(12) + "  "
                              + x2.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void normal_01_samples_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    Normal.normal_01_sample_TEST tests Normal.normal_01_sampleS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int SAMPLE_NUM = 1000;

        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("Normal.normal_01_sampleS_TEST");
        Console.WriteLine("  Normal.normal_01_mean computes the Normal 01 mean;");
        Console.WriteLine("  Normal.normal_01_sampleS samples the Normal 01 distribution;");
        Console.WriteLine("  Normal.normal_01_variance computes the Normal 01 variance;");

        double mean = Normal.normal_01_mean();
        double variance = Normal.normal_01_variance();

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        double[] x = Normal.normal_01_samples(SAMPLE_NUM, ref seed);

        mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
        variance = typeMethods.r8vec_variance(SAMPLE_NUM, x);
        double xmax = typeMethods.r8vec_max(SAMPLE_NUM, x);
        double xmin = typeMethods.r8vec_min(SAMPLE_NUM, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
        Console.WriteLine("  Sample mean =     " + mean + "");
        Console.WriteLine("  Sample variance = " + variance + "");
        Console.WriteLine("  Sample maximum =  " + xmax + "");
        Console.WriteLine("  Sample minimum =  " + xmin + "");

    }

    private static void normal_cdf_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_CDF_TEST tests NORMAL_CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("NORMAL_CDF_TEST");
        Console.WriteLine("  NORMAL_CDF evaluates the Normal CDF;");
        Console.WriteLine("  NORMAL_CDF_INV inverts the Normal CDF.");
        Console.WriteLine("  NORMAL_PDF evaluates the Normal PDF;");

        const double mu = 100.0;
        const double sigma = 15.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter MU =    " + mu + "");
        Console.WriteLine("  PDF parameter SIGMA = " + sigma + "");

        if (!Normal.normal_check(mu, sigma))
        {
            Console.WriteLine("");
            Console.WriteLine("NORMAL_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Normal.normal_sample(mu, sigma, ref seed);
            double pdf = Normal.normal_pdf(x, mu, sigma);
            double cdf = CDF.normal_cdf(x, mu, sigma);
            double x2 = CDF.normal_cdf_inv(cdf, mu, sigma);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    private static void normal_samples_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_SAMPLES_TEST tests NORMAL_SAMPLES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int SAMPLE_NUM = 1000;

        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("NORMAL_SAMPLES_TEST");
        Console.WriteLine("  NORMAL_MEAN computes the Normal mean;");
        Console.WriteLine("  NORMAL_SAMPLES samples the Normal distribution;");
        Console.WriteLine("  NORMAL_VARIANCE computes the Normal variance;");

        const double mu = 100.0;
        const double sigma = 15.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter MU=     " + mu + "");
        Console.WriteLine("  PDF parameter SIGMA = " + sigma + "");

        if (!Normal.normal_check(mu, sigma))
        {
            Console.WriteLine("");
            Console.WriteLine("NORMAL_SAMPLES_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Normal.normal_mean(mu, sigma);
        double variance = Normal.normal_variance(mu, sigma);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        double[] x = Normal.normal_samples(SAMPLE_NUM, mu, sigma, ref seed);

        mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
        variance = typeMethods.r8vec_variance(SAMPLE_NUM, x);
        double xmax = typeMethods.r8vec_max(SAMPLE_NUM, x);
        double xmin = typeMethods.r8vec_min(SAMPLE_NUM, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
        Console.WriteLine("  Sample mean =     " + mean + "");
        Console.WriteLine("  Sample variance = " + variance + "");
        Console.WriteLine("  Sample maximum =  " + xmax + "");
        Console.WriteLine("  Sample minimum =  " + xmin + "");

    }

    private static void normal_truncated_ab_cdf_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_AB_CDF_TEST tests NORMAL_TRUNCATED_AB_CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        const double a = 50.0;
        const double b = 150.0;
        const double mu = 100.0;
        const double s = 25.0;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("NORMAL_TRUNCATED_AB_CDF_TEST");
        Console.WriteLine("  NORMAL_TRUNCATED_AB_CDF evaluates the Normal Truncated AB CDF.");
        Console.WriteLine("  NORMAL_TRUNCATED_AB_CDF_INV inverts the Normal Truncated AB CDF.");
        Console.WriteLine("  NORMAL_TRUNCATED_AB_PDF evaluates the Normal Truncated AB PDF.");
        Console.WriteLine("");
        Console.WriteLine("  The parent normal distribution has");
        Console.WriteLine("    mean =               " + mu + "");
        Console.WriteLine("    standard deviation = " + s + "");
        Console.WriteLine("  The parent distribution is truncated to");
        Console.WriteLine("  the interval [" + a + "," + b + "]");

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Truncated.normal_truncated_ab_sample(mu, s, a, b, ref seed);

            double pdf = Truncated.normal_truncated_ab_pdf(x, mu, s, a, b);

            double cdf = CDF.normal_truncated_ab_cdf(x, mu, s, a, b);

            double x2 = CDF.normal_truncated_ab_cdf_inv(cdf, mu, s, a, b);

            Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + x2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void normal_truncated_ab_sample_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_AB_SAMPLE_TEST tests NORMAL_TRUNCATED_AB_SAMPLE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int sample_num = 1000;
        double[] x = new double[sample_num];

        const double a = 50.0;
        const double b = 150.0;
        const double mu = 100.0;
        const double s = 25.0;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("NORMAL_TRUNCATED_AB_SAMPLE_TEST");
        Console.WriteLine("  NORMAL_TRUNCATED_AB_MEAN computes the Normal Truncated AB mean;");
        Console.WriteLine("  NORMAL_TRUNCATED_AB_SAMPLE samples the Normal Truncated AB distribution;");
        Console.WriteLine("  NORMAL_TRUNCATED_AB_VARIANCE computes the Normal Truncated AB variance.");
        Console.WriteLine("");
        Console.WriteLine("  The parent normal distribution has");
        Console.WriteLine("    mean =               " + mu + "");
        Console.WriteLine("    standard deviation = " + s + "");
        Console.WriteLine("  The parent distribution is truncated to");
        Console.WriteLine("  the interval [" + a + "," + b + "]");

        double mean = Truncated.normal_truncated_ab_mean(mu, s, a, b);

        double variance = Truncated.normal_truncated_ab_variance(mu, s, a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean      =               " + mean + "");
        Console.WriteLine("  PDF variance =                " + variance + "");

        for (i = 0; i < sample_num; i++)
        {
            x[i] = Truncated.normal_truncated_ab_sample(mu, s, a, b, ref seed);
        }

        mean = typeMethods.r8vec_mean(sample_num, x);
        variance = typeMethods.r8vec_variance(sample_num, x);
        double xmax = typeMethods.r8vec_max(sample_num, x);
        double xmin = typeMethods.r8vec_min(sample_num, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + sample_num + "");
        Console.WriteLine("  Sample mean =     " + mean + "");
        Console.WriteLine("  Sample variance = " + variance + "");
        Console.WriteLine("  Sample maximum =  " + xmax + "");
        Console.WriteLine("  Sample minimum =  " + xmin + "");

    }

    private static void normal_truncated_a_cdf_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_A_CDF_TEST tests NORMAL_TRUNCATED_A_CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        const double a = 50.0;
        const double mu = 100.0;
        const double s = 25.0;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("NORMAL_TRUNCATED_A_CDF_TEST");
        Console.WriteLine("  NORMAL_TRUNCATED_A_CDF evaluates the Normal Truncated A CDF.");
        Console.WriteLine("  NORMAL_TRUNCATED_A_CDF_INV inverts the Normal Truncated A CDF.");
        Console.WriteLine("  NORMAL_TRUNCATED_A_PDF evaluates the Normal Truncated A PDF.");
        Console.WriteLine("");
        Console.WriteLine("  The parent normal distribution has");
        Console.WriteLine("    mean =               " + mu + "");
        Console.WriteLine("    standard deviation = " + s + "");
        Console.WriteLine("  The parent distribution is truncated to");
        Console.WriteLine("  the interval [" + a + ",+oo]");

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Truncated.normal_truncated_a_sample(mu, s, a, ref seed);

            double pdf = Truncated.normal_truncated_a_pdf(x, mu, s, a);

            double cdf = CDF.normal_truncated_a_cdf(x, mu, s, a);

            double x2 = CDF.normal_truncated_a_cdf_inv(cdf, mu, s, a);

            Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + x2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

    }

    private static void normal_truncated_a_sample_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_A_SAMPLE_TEST tests NORMAL_TRUNCATED_A_SAMPLE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int sample_num = 1000;
        double[] x = new double[sample_num];

        const double a = 50.0;
        const double mu = 100.0;
        const double s = 25.0;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("NORMAL_TRUNCATED_A_SAMPLE_TEST");
        Console.WriteLine("  NORMAL_TRUNCATED_A_MEAN computes the Normal Truncated A mean;");
        Console.WriteLine("  NORMAL_TRUNCATED_A_SAMPLE samples the Normal Truncated A distribution;");
        Console.WriteLine("  NORMAL_TRUNCATED_A_VARIANCE computes the Normal Truncated A variance.");
        Console.WriteLine("");
        Console.WriteLine("  The parent normal distribution has");
        Console.WriteLine("    mean =               " + mu + "");
        Console.WriteLine("    standard deviation = " + s + "");
        Console.WriteLine("  The parent distribution is truncated to");
        Console.WriteLine("  the interval [" + a + ",+oo]");

        double mean = Truncated.normal_truncated_a_mean(mu, s, a);

        double variance = Truncated.normal_truncated_a_variance(mu, s, a);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean      =               " + mean + "");
        Console.WriteLine("  PDF variance =                " + variance + "");

        for (i = 0; i < sample_num; i++)
        {
            x[i] = Truncated.normal_truncated_a_sample(mu, s, a, ref seed);
        }

        mean = typeMethods.r8vec_mean(sample_num, x);
        variance = typeMethods.r8vec_variance(sample_num, x);
        double xmax = typeMethods.r8vec_max(sample_num, x);
        double xmin = typeMethods.r8vec_min(sample_num, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + sample_num + "");
        Console.WriteLine("  Sample mean =     " + mean + "");
        Console.WriteLine("  Sample variance = " + variance + "");
        Console.WriteLine("  Sample maximum =  " + xmax + "");
        Console.WriteLine("  Sample minimum =  " + xmin + "");

    }

    private static void normal_truncated_b_cdf_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_B_CDF_TEST tests NORMAL_TRUNCATED_B_CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        const double b = 150.0;
        const double mu = 100.0;
        const double s = 25.0;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("NORMAL_TRUNCATED_B_CDF_TEST");
        Console.WriteLine("  NORMAL_TRUNCATED_B_CDF evaluates the Normal Truncated B CDF.");
        Console.WriteLine("  NORMAL_TRUNCATED_B_CDF_INV inverts the Normal Truncated B CDF.");
        Console.WriteLine("  NORMAL_TRUNCATED_B_PDF evaluates the Normal Truncated B PDF.");
        Console.WriteLine("");
        Console.WriteLine("  The parent normal distribution has");
        Console.WriteLine("    mean =               " + mu + "");
        Console.WriteLine("    standard deviation = " + s + "");
        Console.WriteLine("  The parent distribution is truncated to");
        Console.WriteLine("  the interval [-oo," + b + "]");

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Truncated.normal_truncated_b_sample(mu, s, b, ref seed);

            double pdf = Truncated.normal_truncated_b_pdf(x, mu, s, b);

            double cdf = CDF.normal_truncated_b_cdf(x, mu, s, b);

            double x2 = CDF.normal_truncated_b_cdf_inv(cdf, mu, s, b);

            Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + x2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void normal_truncated_b_sample_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_B_SAMPLE_TEST tests NORMAL_TRUNCATED_B_SAMPLE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int sample_num = 1000;
        double[] x = new double[sample_num];

        const double b = 150.0;
        const double mu = 100.0;
        const double s = 25.0;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("NORMAL_TRUNCATED_B_SAMPLE_TEST");
        Console.WriteLine("  NORMAL_TRUNCATED_B_MEAN computes the Normal Truncated B mean;");
        Console.WriteLine("  NORMAL_TRUNCATED_B_SAMPLE samples the Normal Truncated B distribution;");
        Console.WriteLine("  NORMAL_TRUNCATED_B_VARIANCE computes the Normal Truncated B variance.");
        Console.WriteLine("");
        Console.WriteLine("  The parent normal distribution has");
        Console.WriteLine("    mean =               " + mu + "");
        Console.WriteLine("    standard deviation = " + s + "");
        Console.WriteLine("  The parent distribution is truncated to");
        Console.WriteLine("  the interval [-oo," + b + "]");

        double mean = Truncated.normal_truncated_b_mean(mu, s, b);

        double variance = Truncated.normal_truncated_b_variance(mu, s, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean      =               " + mean + "");
        Console.WriteLine("  PDF variance =                " + variance + "");

        for (i = 0; i < sample_num; i++)
        {
            x[i] = Truncated.normal_truncated_b_sample(mu, s, b, ref seed);
        }

        mean = typeMethods.r8vec_mean(sample_num, x);
        variance = typeMethods.r8vec_variance(sample_num, x);
        double xmax = typeMethods.r8vec_max(sample_num, x);
        double xmin = typeMethods.r8vec_min(sample_num, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + sample_num + "");
        Console.WriteLine("  Sample mean =     " + mean + "");
        Console.WriteLine("  Sample variance = " + variance + "");
        Console.WriteLine("  Sample maximum =  " + xmax + "");
        Console.WriteLine("  Sample minimum =  " + xmin + "");

    }

}