using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest
{
    partial class Program
    {
        static void normal_01_cdf_test()

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
            double cdf;
            int i;
            double pdf;
            int seed = 123456789;
            double x;
            double x2;

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
                x = Normal.normal_01_sample(ref seed);
                pdf = Normal.normal_01_pdf(x);
                cdf = Normal.normal_01_cdf(x);
                x2 = Normal.normal_01_cdf_inv(cdf);

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + pdf.ToString("0.######").PadLeft(12) + "  "
                                  + cdf.ToString("0.####").PadLeft(12) + "  "
                                  + x2.ToString("0.################").PadLeft(24) + "");
            }
        }

        static void normal_01_samples_test()

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
            int SAMPLE_NUM = 1000;

            double mean;
            int seed = 123456789;
            double variance;
            double[] x;
            double xmax;
            double xmin;

            Console.WriteLine("");
            Console.WriteLine("Normal.normal_01_sampleS_TEST");
            Console.WriteLine("  Normal.normal_01_mean computes the Normal 01 mean;");
            Console.WriteLine("  Normal.normal_01_sampleS samples the Normal 01 distribution;");
            Console.WriteLine("  Normal.normal_01_variance computes the Normal 01 variance;");

            mean = Normal.normal_01_mean();
            variance = Normal.normal_01_variance();

            Console.WriteLine("");
            Console.WriteLine("  PDF mean =     " + mean + "");
            Console.WriteLine("  PDF variance = " + variance + "");

            x = Normal.normal_01_samples(SAMPLE_NUM, ref seed);

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

        static void normal_cdf_test()

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
            double cdf;
            int i;
            double mu;
            double pdf;
            int seed = 123456789;
            double sigma;
            double x;
            double x2;

            Console.WriteLine("");
            Console.WriteLine("NORMAL_CDF_TEST");
            Console.WriteLine("  NORMAL_CDF evaluates the Normal CDF;");
            Console.WriteLine("  NORMAL_CDF_INV inverts the Normal CDF.");
            Console.WriteLine("  NORMAL_PDF evaluates the Normal PDF;");

            mu = 100.0;
            sigma = 15.0;

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
                x = Normal.normal_sample(mu, sigma, ref seed);
                pdf = Normal.normal_pdf(x, mu, sigma);
                cdf = Normal.normal_cdf(x, mu, sigma);
                x2 = Normal.normal_cdf_inv(cdf, mu, sigma);

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + pdf.ToString().PadLeft(12) + "  "
                                  + cdf.ToString().PadLeft(12) + "  "
                                  + x2.ToString().PadLeft(12) + "");
            }

        }

        static void normal_samples_test()

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
            int SAMPLE_NUM = 1000;

            double mean;
            double mu;
            int seed = 123456789;
            double sigma;
            double variance;
            double[] x;
            double xmax;
            double xmin;

            Console.WriteLine("");
            Console.WriteLine("NORMAL_SAMPLES_TEST");
            Console.WriteLine("  NORMAL_MEAN computes the Normal mean;");
            Console.WriteLine("  NORMAL_SAMPLES samples the Normal distribution;");
            Console.WriteLine("  NORMAL_VARIANCE computes the Normal variance;");

            mu = 100.0;
            sigma = 15.0;

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

            mean = Normal.normal_mean(mu, sigma);
            variance = Normal.normal_variance(mu, sigma);

            Console.WriteLine("");
            Console.WriteLine("  PDF mean =     " + mean + "");
            Console.WriteLine("  PDF variance = " + variance + "");

            x = Normal.normal_samples(SAMPLE_NUM, mu, sigma, ref seed);

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

        static void normal_truncated_ab_cdf_test()

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
            double a;
            double b;
            double cdf;
            int i;
            double mu;
            double pdf;
            double s;
            int seed;
            double x;
            double x2;

            a = 50.0;
            b = 150.0;
            mu = 100.0;
            s = 25.0;
            seed = 123456789;

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
                x = Normal.normal_truncated_ab_sample(mu, s, a, b, ref seed);

                pdf = Normal.normal_truncated_ab_pdf(x, mu, s, a, b);

                cdf = Normal.normal_truncated_ab_cdf(x, mu, s, a, b);

                x2 = Normal.normal_truncated_ab_cdf_inv(cdf, mu, s, a, b);

                Console.WriteLine("  " + x.ToString().PadLeft(14)
                                  + "  " + pdf.ToString().PadLeft(14)
                                  + "  " + cdf.ToString().PadLeft(14)
                                  + "  " + x2.ToString().PadLeft(14) + "");
            }

            return;
        }

        static void normal_truncated_ab_sample_test()

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
            double a;
            double b;
            int i;
            double mean;
            double mu;
            double s;
            int sample_num = 1000;
            int seed;
            double variance;
            double[] x = new double[sample_num];
            double xmax;
            double xmin;

            a = 50.0;
            b = 150.0;
            mu = 100.0;
            s = 25.0;
            seed = 123456789;

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

            mean = Normal.normal_truncated_ab_mean(mu, s, a, b);

            variance = Normal.normal_truncated_ab_variance(mu, s, a, b);

            Console.WriteLine("");
            Console.WriteLine("  PDF mean      =               " + mean + "");
            Console.WriteLine("  PDF variance =                " + variance + "");

            for (i = 0; i < sample_num; i++)
            {
                x[i] = Normal.normal_truncated_ab_sample(mu, s, a, b, ref seed);
            }

            mean = typeMethods.r8vec_mean(sample_num, x);
            variance = typeMethods.r8vec_variance(sample_num, x);
            xmax = typeMethods.r8vec_max(sample_num, x);
            xmin = typeMethods.r8vec_min(sample_num, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + sample_num + "");
            Console.WriteLine("  Sample mean =     " + mean + "");
            Console.WriteLine("  Sample variance = " + variance + "");
            Console.WriteLine("  Sample maximum =  " + xmax + "");
            Console.WriteLine("  Sample minimum =  " + xmin + "");

        }

        static void normal_truncated_a_cdf_test()

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
            double a;
            double cdf;
            int i;
            double mu;
            double pdf;
            double s;
            int seed;
            double x;
            double x2;

            a = 50.0;
            mu = 100.0;
            s = 25.0;
            seed = 123456789;

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
                x = Normal.normal_truncated_a_sample(mu, s, a, ref seed);

                pdf = Normal.normal_truncated_a_pdf(x, mu, s, a);

                cdf = Normal.normal_truncated_a_cdf(x, mu, s, a);

                x2 = Normal.normal_truncated_a_cdf_inv(cdf, mu, s, a);

                Console.WriteLine("  " + x.ToString().PadLeft(14)
                                  + "  " + pdf.ToString().PadLeft(14)
                                  + "  " + cdf.ToString().PadLeft(14)
                                  + "  " + x2.ToString().PadLeft(14) + "");
            }

        }

        static void normal_truncated_a_sample_test()

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
            double a;
            int i;
            double mean;
            double mu;
            double s;
            int sample_num = 1000;
            int seed;
            double variance;
            double[] x = new double[sample_num];
            double xmax;
            double xmin;

            a = 50.0;
            mu = 100.0;
            s = 25.0;
            seed = 123456789;

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

            mean = Normal.normal_truncated_a_mean(mu, s, a);

            variance = Normal.normal_truncated_a_variance(mu, s, a);

            Console.WriteLine("");
            Console.WriteLine("  PDF mean      =               " + mean + "");
            Console.WriteLine("  PDF variance =                " + variance + "");

            for (i = 0; i < sample_num; i++)
            {
                x[i] = Normal.normal_truncated_a_sample(mu, s, a, ref seed);
            }

            mean = typeMethods.r8vec_mean(sample_num, x);
            variance = typeMethods.r8vec_variance(sample_num, x);
            xmax = typeMethods.r8vec_max(sample_num, x);
            xmin = typeMethods.r8vec_min(sample_num, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + sample_num + "");
            Console.WriteLine("  Sample mean =     " + mean + "");
            Console.WriteLine("  Sample variance = " + variance + "");
            Console.WriteLine("  Sample maximum =  " + xmax + "");
            Console.WriteLine("  Sample minimum =  " + xmin + "");

        }

        static void normal_truncated_b_cdf_test()

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
            double b;
            double cdf;
            int i;
            double mu;
            double pdf;
            double s;
            int seed;
            double x;
            double x2;

            b = 150.0;
            mu = 100.0;
            s = 25.0;
            seed = 123456789;

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
                x = Normal.normal_truncated_b_sample(mu, s, b, ref seed);

                pdf = Normal.normal_truncated_b_pdf(x, mu, s, b);

                cdf = Normal.normal_truncated_b_cdf(x, mu, s, b);

                x2 = Normal.normal_truncated_b_cdf_inv(cdf, mu, s, b);

                Console.WriteLine("  " + x.ToString().PadLeft(14)
                                  + "  " + pdf.ToString().PadLeft(14)
                                  + "  " + cdf.ToString().PadLeft(14)
                                  + "  " + x2.ToString().PadLeft(14) + "");
            }

            return;
        }

        static void normal_truncated_b_sample_test()

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
            double b;
            int i;
            double mean;
            double mu;
            double s;
            int sample_num = 1000;
            int seed;
            double variance;
            double[] x = new double[sample_num];
            double xmax;
            double xmin;

            b = 150.0;
            mu = 100.0;
            s = 25.0;
            seed = 123456789;

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

            mean = Normal.normal_truncated_b_mean(mu, s, b);

            variance = Normal.normal_truncated_b_variance(mu, s, b);

            Console.WriteLine("");
            Console.WriteLine("  PDF mean      =               " + mean + "");
            Console.WriteLine("  PDF variance =                " + variance + "");

            for (i = 0; i < sample_num; i++)
            {
                x[i] = Normal.normal_truncated_b_sample(mu, s, b, ref seed);
            }

            mean = typeMethods.r8vec_mean(sample_num, x);
            variance = typeMethods.r8vec_variance(sample_num, x);
            xmax = typeMethods.r8vec_max(sample_num, x);
            xmin = typeMethods.r8vec_min(sample_num, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + sample_num + "");
            Console.WriteLine("  Sample mean =     " + mean + "");
            Console.WriteLine("  Sample variance = " + variance + "");
            Console.WriteLine("  Sample maximum =  " + xmax + "");
            Console.WriteLine("  Sample minimum =  " + xmin + "");

        }

    }
}