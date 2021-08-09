using System;
using Burkardt.CDFLib;
using Burkardt.Probability;
using Burkardt.Types;
using Burkardt.Uniform;

namespace TruncatedNormalTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TRUNCATED_NORMAL_TEST.
            //
            //  Discussion:
            //
            //    TRUNCATED_NORMAL_TEST tests the TRUNCATED_NORMAL library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_TEST");
            Console.WriteLine("  Test the TRUNCATED_NORMAL library.");
            //
            //  Utilities.
            //
            i4_uniform_ab_test();
            r8_choose_test();
            r8_factorial2_test();
            r8_mop_test();
            r8_uniform_01_test();
            r8poly_print_test();
            r8poly_value_horner_test();
            r8vec_linspace_new_test();
            r8vec_print_test();
            //
            //  Library.
            //
            normal_01_cdf_test();
            normal_01_cdf_inv_test();
            normal_01_mean_test();
            normal_01_moment_test();
            normal_01_pdf_test();
            normal_01_sample_test();
            normal_01_variance_test();

            normal_ms_cdf_test();
            normal_ms_cdf_inv_test();
            normal_ms_mean_test();
            normal_ms_moment_test();
            normal_ms_moment_central_test();
            normal_ms_pdf_test();
            normal_ms_sample_test();
            normal_ms_variance_test();

            truncated_normal_a_cdf_test();
            truncated_normal_a_cdf_inv_test();
            truncated_normal_a_mean_test();
            truncated_normal_a_moment_test();
            truncated_normal_a_pdf_test();
            truncated_normal_a_sample_test();
            truncated_normal_a_variance_test();

            truncated_normal_ab_cdf_test();
            truncated_normal_ab_cdf_inv_test();
            truncated_normal_ab_mean_test();
            truncated_normal_ab_moment_test();
            truncated_normal_ab_pdf_test();
            truncated_normal_ab_sample_test();
            truncated_normal_ab_variance_test();

            truncated_normal_b_cdf_test();
            truncated_normal_b_cdf_inv_test();
            truncated_normal_b_mean_test();
            truncated_normal_b_moment_test();
            truncated_normal_b_pdf_test();
            truncated_normal_b_sample_test();
            truncated_normal_b_variance_test();
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void i4_uniform_ab_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_UNIFORM_AB_TEST tests I4_UNIFORM_AB.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int a = -100;
            int b = 200;
            int i;
            int j;
            int seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("I4_UNIFORM_AB_TEST");
            Console.WriteLine("  I4_UNIFORM_AB computes pseudorandom values");
            Console.WriteLine("  in an interval [A,B].");

            Console.WriteLine("");
            Console.WriteLine("  The lower endpoint A = " + a + "");
            Console.WriteLine("  The upper endpoint B = " + b + "");
            Console.WriteLine("  The initial seed is " + seed + "");
            Console.WriteLine("");

            for (i = 1; i <= 20; i++)
            {
                j = UniformRNG.i4_uniform_ab(a, b, ref seed);

                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                       + "  " + j.ToString().PadLeft(8) + "");
            }

        }

        static void normal_01_cdf_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_01_CDF_TEST tests NORMAL_01_CDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 February 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double cdf1 = 0;
            double cdf2;
            int n_data;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("NORMAL_01_CDF_TEST");
            Console.WriteLine("  NORMAL_01_CDF evaluates the Normal 01 CDF;");
            Console.WriteLine("");
            Console.WriteLine("       X              CDF                       CDF");
            Console.WriteLine("                     (exact)                   (computed)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                TestValues.Normal.normal_01_cdf_values(ref n_data, ref x, ref cdf1);

                if (n_data == 0)
                {
                    break;
                }

                cdf2 = CDF.normal_01_cdf(x);

                Console.WriteLine("  " + x.ToString().PadLeft(14)
                                       + "  " + cdf1.ToString().PadLeft(24)
                                       + "  " + cdf2.ToString().PadLeft(24) + "");
            }

        }

        static void normal_01_cdf_inv_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_01_CDF_INV_TEST tests NORMAL_01_CDF_INV.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 February 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double cdf = 0;
            int n_data;
            double x1 = 0;
            double x2;

            Console.WriteLine("");
            Console.WriteLine("NORMAL_01_CDF_INV_TEST");
            Console.WriteLine("  NORMAL_01_CDF_INV inverts the Normal 01 CDF;");
            Console.WriteLine("");
            Console.WriteLine("      CDF             X                         X");
            Console.WriteLine("                     (exact)                   (computed)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                TestValues.Normal.normal_01_cdf_values(ref n_data, ref x1, ref cdf);

                if (n_data == 0)
                {
                    break;
                }

                x2 = CDF.normal_01_cdf_inv(cdf);

                Console.WriteLine("  " + cdf.ToString().PadLeft(14)
                                       + "  " + x1.ToString().PadLeft(24)
                                       + "  " + x2.ToString().PadLeft(24) + "");
            }

        }

        static void normal_01_mean_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_01_MEAN_TEST tests NORMAL_01_MEAN.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            double mean;
            int sample_num;
            int seed = 123456789;
            double[] x;
            double xmax;
            double xmin;

            Console.WriteLine("");
            Console.WriteLine("NORMAL_01_MEAN_TEST");
            Console.WriteLine("  NORMAL_01_MEAN computes the Normal 01 mean;");

            mean = Normal.normal_01_mean();

            Console.WriteLine("");
            Console.WriteLine("  PDF mean =     " + mean + "");

            sample_num = 1000;
            x = new double[sample_num];

            for (i = 0; i < sample_num; i++)
            {
                x[i] = Normal.normal_01_sample(ref seed);
            }

            mean = typeMethods.r8vec_mean(sample_num, x);
            xmax = typeMethods.r8vec_max(sample_num, x);
            xmin = typeMethods.r8vec_min(sample_num, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + sample_num + "");
            Console.WriteLine("  Sample mean =     " + mean + "");
            Console.WriteLine("  Sample maximum =  " + xmax + "");
            Console.WriteLine("  Sample minimum =  " + xmin + "");


        }

        static void normal_01_moment_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_01_MOMENT_TEST tests NORMAL_01_MOMENT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double moment;
            int order;

            Console.WriteLine("");
            Console.WriteLine("NORMAL_01_MOMENT_TEST");
            Console.WriteLine("  NORMAL_01_MOMENT evaluates Normal 01 moments;");
            Console.WriteLine("");
            Console.WriteLine("      Order              Moment");
            Console.WriteLine("");

            for (order = 0; order <= 10; order++)
            {
                moment = Normal.normal_01_moment(order);

                Console.WriteLine("  " + order.ToString().PadLeft(14)
                                       + "  " + moment.ToString().PadLeft(14) + "");
            }

        }

        static void normal_01_pdf_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_01_PDF_TEST tests NORMAL_01_PDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            double pdf;
            double x;

            Console.WriteLine("");
            Console.WriteLine("NORMAL_01_PDF_TEST");
            Console.WriteLine("  NORMAL_01_PDF evaluates the Normal 01 PDF;");
            Console.WriteLine("");
            Console.WriteLine("       X              PDF");
            Console.WriteLine("");

            for (i = -20; i <= 20; i++)
            {
                x = (double)(i) / 10.0;

                pdf = Normal.normal_01_pdf(x);

                Console.WriteLine("  " + x.ToString().PadLeft(14)
                                       + "  " + pdf.ToString().PadLeft(14) + "");
            }

        }

        static void normal_01_sample_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_01_SAMPLE_TEST tests NORMAL_01_SAMPLE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int seed;
            double x;

            Console.WriteLine("");
            Console.WriteLine("NORMAL_01_SAMPLE_TEST");
            Console.WriteLine("  NORMAL_01_SAMPLE returns samples from the normal");
            Console.WriteLine("  distribution with mean 0 and standard deviation 1.");
            Console.WriteLine("");

            seed = 123456789;

            for (i = 1; i <= 10; i++)
            {
                x = Normal.normal_01_sample(ref seed);
                Console.WriteLine("  " + i.ToString().PadLeft(4)
                                       + "  " + x.ToString().PadLeft(14) + "");
            }

        }

        static void normal_01_variance_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_01_VARIANCE_TEST tests NORMAL_01_VARIANCE;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int sample_num;
            int seed = 123456789;
            double variance;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("NORMAL_01_VARIANCE_TEST");
            Console.WriteLine("  NORMAL_01_VARIANCE computes the Normal 01 variance;");

            variance = Normal.normal_01_variance();

            Console.WriteLine("");
            Console.WriteLine("  PDF variance = " + variance + "");

            sample_num = 1000;

            x = new double[sample_num];

            for (i = 0; i < sample_num; i++)
            {
                x[i] = Normal.normal_01_sample(ref seed);
            }

            variance = typeMethods.r8vec_variance(sample_num, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + sample_num + "");
            Console.WriteLine("  Sample variance = " + variance + "");


        }

        static void normal_ms_cdf_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_MS_CDF_TEST tests NORMAL_MS_CDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double cdf;
            int i;
            double mu;
            double sigma;
            double x;

            Console.WriteLine("");
            Console.WriteLine("NORMAL_MS_CDF_TEST");
            Console.WriteLine("  NORMAL_MS_CDF evaluates the Normal MS CDF;");

            mu = 100.0;
            sigma = 15.0;

            Console.WriteLine("");
            Console.WriteLine("  Parameter MU = " + mu + "");
            Console.WriteLine("  Parameteter SIGMA = " + sigma + "");
            Console.WriteLine("");
            Console.WriteLine("       X              CDF");
            Console.WriteLine("");

            for (i = -20; i <= +20; i++)
            {
                x = mu + sigma * (double)(i) / 10.0;
                cdf = CDF.normal_ms_cdf(x, mu, sigma);
                Console.WriteLine("  " + x.ToString().PadLeft(14)
                                       + "  " + cdf.ToString().PadLeft(24) + "");
            }

        }

        static void normal_ms_cdf_inv_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_MS_CDF_INV_TEST tests NORMAL_MS_CDF_INV.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double cdf;
            int i;
            double mu;
            double sigma;
            double x;
            double x2;

            Console.WriteLine("");
            Console.WriteLine("NORMAL_MS_CDF_INV_TEST");
            Console.WriteLine("  NORMAL_MS_CDF_INV inverts the Normal MS CDF;");

            mu = 100.0;
            sigma = 15.0;

            Console.WriteLine("");
            Console.WriteLine("  Parameter MU = " + mu + "");
            Console.WriteLine("  Parameteter SIGMA = " + sigma + "");

            Console.WriteLine("");
            Console.WriteLine("       X            CDF           CDF_INV");
            Console.WriteLine("");

            for (i = -20; i <= +20; i++)
            {
                x = mu + sigma * (double)(i) / 10.0;
                cdf = CDF.normal_ms_cdf(x, mu, sigma);
                x2 = CDF.normal_ms_cdf_inv(cdf, mu, sigma);
                Console.WriteLine("  " + x.ToString().PadLeft(14)
                                       + "  " + cdf.ToString().PadLeft(14)
                                       + "  " + x2.ToString().PadLeft(14) + "");
            }

        }

        static void normal_ms_mean_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_MS_MEAN_TEST tests NORMAL_MS_MEAN.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            double mean;
            double mu;
            int sample_num;
            int seed = 123456789;
            double sigma;
            double[] x;
            double xmax;
            double xmin;

            Console.WriteLine("");
            Console.WriteLine("NORMAL_MS_MEAN_TEST");
            Console.WriteLine("  NORMAL_MS_MEAN computes the Normal MS mean.");

            mu = 100.0;
            sigma = 15.0;

            Console.WriteLine("");
            Console.WriteLine("  Parameter MU = " + mu + "");
            Console.WriteLine("  Parameteter SIGMA = " + sigma + "");

            mean = Normal.normal_ms_mean(mu, sigma);

            Console.WriteLine("");
            Console.WriteLine("  PDF mean = " + mean + "");

            sample_num = 1000;
            x = new double[sample_num];

            for (i = 0; i < sample_num; i++)
            {
                x[i] = Normal.normal_ms_sample(mu, sigma, ref seed);
            }

            mean = typeMethods.r8vec_mean(sample_num, x);
            xmax = typeMethods.r8vec_max(sample_num, x);
            xmin = typeMethods.r8vec_min(sample_num, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + sample_num + "");
            Console.WriteLine("  Sample mean =     " + mean + "");
            Console.WriteLine("  Sample maximum =  " + xmax + "");
            Console.WriteLine("  Sample minimum =  " + xmin + "");


        }

        static void normal_ms_moment_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_MOMENT_MS_TEST tests NORMAL_MS_MOMENT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 August 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double moment1;
            double moment2;
            double mu;
            double[] mu_test = { 0.0, 2.0, 10.0, 0.0 };
            int order;
            double sigma;
            double[] sigma_test = { 1.0, 1.0, 2.0, 2.0 };
            int test;
            int test_num = 4;

            Console.WriteLine("");
            Console.WriteLine("NORMAL_MOMENT_MS_TEST");
            Console.WriteLine("  NORMAL_MS_MOMENT evaluates the moments of the Normal MS distribution.");

            for (test = 0; test < test_num; test++)
            {
                mu = mu_test[test];
                sigma = sigma_test[test];
                Console.WriteLine("");
                Console.WriteLine("  Mu = " + mu
                                            + "  Sigma = " + sigma + "");
                Console.WriteLine(" Order  Moment");
                Console.WriteLine("");
                for (order = 0; order <= 8; order++)
                {
                    moment1 = Normal.normal_ms_moment(order, mu, sigma);
                    moment2 = Normal.normal_ms_moment_values(order, mu, sigma);
                    Console.WriteLine("  " + order.ToString().PadLeft(2)
                                           + "  " + moment1.ToString().PadLeft(14)
                                           + "  " + moment2.ToString().PadLeft(14) + "");
                }
            }

        }

        static void normal_ms_moment_central_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_MS_MOMENT_CENTRAL_TEST tests NORMAL_MS_MOMENT_CENTRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 August 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double moment1;
            double moment2;
            double mu;
            double[] mu_test = { 0.0, 2.0, 10.0, 0.0 };
            int order;
            double sigma;
            double[] sigma_test = { 1.0, 1.0, 2.0, 2.0 };
            int test;
            int test_num = 4;

            Console.WriteLine("");
            Console.WriteLine("NORMAL_MS_MOMENT_CENTRAL_TEST");
            Console.WriteLine("  NORMAL_MS_MOMENT_CENTRAL evaluates the central moments of the");
            Console.WriteLine("  Normal MS distribution.");

            for (test = 0; test < test_num; test++)
            {
                mu = mu_test[test];
                sigma = sigma_test[test];
                Console.WriteLine("");
                Console.WriteLine("  Mu = " + mu
                                            + "  Sigma = " + sigma + "");
                Console.WriteLine(" Order  Moment");
                Console.WriteLine("");
                for (order = 0; order <= 8; order++)
                {
                    moment1 = Normal.normal_ms_moment_central(order, mu, sigma);
                    moment2 = Normal.normal_ms_moment_central_values(order, mu, sigma);
                    Console.WriteLine("  " + order.ToString().PadLeft(2)
                                           + "  " + moment1.ToString().PadLeft(14)
                                           + "  " + moment2.ToString().PadLeft(14) + "");
                }
            }

        }

        static void normal_ms_pdf_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_MS_PDF_TEST tests NORMAL_MS_PDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            double mu;
            double pdf;
            double sigma;
            double x;

            Console.WriteLine("");
            Console.WriteLine("NORMAL_MS_PDF_TEST");
            Console.WriteLine("  NORMAL_MS_PDF evaluates the Normal MS PDF;");

            mu = 100.0;
            sigma = 15.0;

            Console.WriteLine("");
            Console.WriteLine("  Parameter MU = " + mu + "");
            Console.WriteLine("  Parameteter SIGMA = " + sigma + "");

            Console.WriteLine("");
            Console.WriteLine("       X              PDF");
            Console.WriteLine("");

            for (i = -20; i <= +20; i++)
            {
                x = mu + sigma * (double)(i) / 10.0;
                pdf = Normal.normal_ms_pdf(mu, sigma, x);
                Console.WriteLine("  " + x.ToString().PadLeft(14)
                                       + "  " + pdf.ToString().PadLeft(24) + "");
            }

        }

        static void normal_ms_sample_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_MS_SAMPLE_TEST tests NORMAL_MS_SAMPLE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            double mu;
            int seed;
            double sigma;
            double x;

            Console.WriteLine("");
            Console.WriteLine("NORMAL_MS_SAMPLE_TEST");
            Console.WriteLine("  NORMAL_MS_SAMPLE returns samples from the Normal MS PDF.");

            mu = 100.0;
            sigma = 15.0;

            Console.WriteLine("");
            Console.WriteLine("  Parameter MU = " + mu + "");
            Console.WriteLine("  Parameteter SIGMA = " + sigma + "");

            Console.WriteLine("");

            seed = 123456789;

            for (i = 1; i <= 10; i++)
            {
                x = Normal.normal_ms_sample(mu, sigma, ref seed);
                Console.WriteLine("  " + i.ToString().PadLeft(4)
                                       + "  " + x.ToString().PadLeft(14) + "");
            }

        }

        static void normal_ms_variance_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_MS_VARIANCE_TEST tests NORMAL_MS_VARIANCE;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            double mu;
            int sample_num;
            int seed = 123456789;
            double sigma;
            double variance;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("NORMAL_MS_VARIANCE_TEST");
            Console.WriteLine("  NORMAL_MS_VARIANCE computes the Normal MS variance;");

            mu = 100.0;
            sigma = 15.0;

            Console.WriteLine("");
            Console.WriteLine("  Parameter MU = " + mu + "");
            Console.WriteLine("  Parameteter SIGMA = " + sigma + "");

            variance = Normal.normal_ms_variance(mu, sigma);

            Console.WriteLine("");
            Console.WriteLine("  PDF variance = " + variance + "");

            sample_num = 1000;
            x = new double[sample_num];

            for (i = 0; i < sample_num; i++)
            {
                x[i] = Normal.normal_ms_sample(mu, sigma, ref seed);
            }

            variance = typeMethods.r8vec_variance(sample_num, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + sample_num + "");
            Console.WriteLine("  Sample variance = " + variance + "");


        }

        static void r8_choose_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_CHOOSE_TEST tests R8_CHOOSE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    26 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double cnk;
            int k;
            int n;

            Console.WriteLine("");
            Console.WriteLine("R8_CHOOSE_TEST");
            Console.WriteLine("  R8_CHOOSE evaluates C(N,K).");
            Console.WriteLine("");
            Console.WriteLine("         N         K       CNK");

            for (n = 0; n <= 5; n++)
            {
                Console.WriteLine("");
                for (k = 0; k <= n; k++)
                {
                    cnk = typeMethods.r8_choose(n, k);
                    Console.WriteLine(n.ToString().PadLeft(10) + "  "
                                                               + k.ToString().PadLeft(8) + "  "
                                                               + cnk.ToString().PadLeft(14) + "");
                }
            }

        }

        static void r8_factorial2_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_FACTORIAL2_TEST tests R8_FACTORIAL2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 February 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double f1 = 0;
            double f2;
            int n = 0;
            int n_data;

            Console.WriteLine("");
            Console.WriteLine("R8_FACTORIAL2_TEST");
            Console.WriteLine("  R8_FACTORIAL2 evaluates the double factorial function.");
            Console.WriteLine("");
            Console.WriteLine("    N                Exact                  Computed");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                typeMethods.r8_factorial2_values(ref n_data, ref n, ref f1);

                if (n_data == 0)
                {
                    break;
                }

                f2 = typeMethods.r8_factorial2(n);

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(4) + "  "
                                  + f1.ToString("0.################").PadLeft(24) + "  "
                                  + f2.ToString("0.################").PadLeft(24) + "");
            }

        }

        static void r8_mop_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_MOP_TEST tests R8_MOP.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 December 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i4;
            int i4_max;
            int i4_min;
            double r8;
            int seed = 123456789;
            int test;

            Console.WriteLine("");
            Console.WriteLine("R8_MOP_TEST");
            Console.WriteLine("  R8_MOP evaluates (-1.0)^I4 as an R8.");
            Console.WriteLine("");
            Console.WriteLine("    I4  R8_MOP(I4)");
            Console.WriteLine("");

            i4_min = -100;
            i4_max = +100;

            for (test = 1; test <= 10; test++)
            {
                i4 = UniformRNG.i4_uniform_ab(i4_min, i4_max, ref seed);
                r8 = typeMethods.r8_mop(i4);
                Console.WriteLine("  "
                                  + i4.ToString().PadLeft(4) + "  "
                                  + r8.ToString().PadLeft(4) + "");
            }

        }

        static void r8_uniform_01_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_UNIFORM_01_TEST tests R8_UNIFORM_01.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 1000;

            int i;
            double max;
            double mean;
            double min;
            int seed = 123456789;
            double[] x = new double[N];
            double variance;

            Console.WriteLine("");
            Console.WriteLine("R8_UNIFORM_01_TEST");
            Console.WriteLine("  R8_UNIFORM_01 samples a uniform random distribution in [0,1].");
            Console.WriteLine("  distributed random numbers.");
            Console.WriteLine("  Using initial random number seed = " + seed + "");

            for (i = 0; i < N; i++)
            {
                x[i] = UniformRNG.r8_uniform_01(ref seed);
            }

            Console.WriteLine("");
            Console.WriteLine("  First few values:");
            Console.WriteLine("");
            for (i = 0; i < 10; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(6)
                                       + "  " + x[i].ToString().PadLeft(14) + "");
            }

            min = typeMethods.r8vec_min(N, x);
            max = typeMethods.r8vec_max(N, x);
            mean = typeMethods.r8vec_mean(N, x);
            variance = typeMethods.r8vec_variance(N, x);

            Console.WriteLine("");
            Console.WriteLine("  Number of samples was " + N + "");
            Console.WriteLine("  Minimum value was " + min + "");
            Console.WriteLine("  Maximum value was " + max + "");
            Console.WriteLine("  Average value was " + mean + "");
            Console.WriteLine("  Variance was      " + variance + "");
        }

        static void r8poly_print_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_PRINT_TEST tests R8POLY_PRINT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 January 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] c = { 2.0, -3.4, 56.0, 0.0, 0.78, 9.0 };
            int m = 5;

            Console.WriteLine("");
            Console.WriteLine("R8POLY_PRINT_TEST");
            Console.WriteLine("  R8POLY_PRINT prints an R8POLY.");

            typeMethods.r8poly_print(m, c, "  The R8POLY:");

        }

        static void r8poly_value_horner_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_VALUE_HORNER_TEST tests R8POLY_VALUE_HORNER.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 January 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] c = { 24.0, -50.0, +35.0, -10.0, 1.0 };
            int i;
            int m = 4;
            int n = 16;
            double p;
            double[] x;
            double x_hi;
            double x_lo;

            Console.WriteLine("");
            Console.WriteLine("R8POLY_VALUE_HORNER_TEST");
            Console.WriteLine("  R8POLY_VALUE_HORNER evaluates a polynomial at");
            Console.WriteLine("  one point, using Horner's method.");

            typeMethods.r8poly_print(m, c, "  The polynomial coefficients:");

            x_lo = 0.0;
            x_hi = 5.0;
            x = typeMethods.r8vec_linspace_new(n, x_lo, x_hi);

            Console.WriteLine("");
            Console.WriteLine("   I    X    P(X)");
            Console.WriteLine("");

            for (i = 0; i < n; i++)
            {
                p = typeMethods.r8poly_value_horner(m, c, x[i]);
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                                       + "  " + x[i].ToString().PadLeft(8)
                                       + "  " + p.ToString().PadLeft(14) + "");
            }


        }

        static void r8vec_linspace_new_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_LINSPACE_NEW_TEST tests R8VEC_LINSPACE_NEW.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            double b;
            int n = 5;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("R8VEC_LINSPACE_NEW_TEST");
            Console.WriteLine("  For a R8VEC:");
            Console.WriteLine("  R8VEC_LINSPACE_NEW: evenly spaced points between A and B;");

            a = 10.0;
            b = 20.0;

            x = typeMethods.r8vec_linspace_new(n, a, b);
            typeMethods.r8vec_print(n, x, "  r8vec_linspace ( 5, 10, 20 )");

        }

        static void r8vec_print_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_PRINT_TEST tests R8VEC_PRINT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] a = { 123.456, 0.000005, -1.0E+06, 3.14159265 };
            int n = 4;

            Console.WriteLine("");
            Console.WriteLine("TEST1335");
            Console.WriteLine("  R8VEC_PRINT prints an R8VEC.");

            typeMethods.r8vec_print(n, a, "  The R8VEC:");

        }

        static void truncated_normal_a_cdf_test()

            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    TRUNCATED_NORMAL_A_CDF_TEST tests TRUNCATED_NORMAL_A_CDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 September 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double cdf1 = 0;
            double cdf2;
            double mu = 0;
            int n_data;
            double sigma = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_A_CDF_TEST:");
            Console.WriteLine("  TRUNCATED_NORMAL_A_CDF evaluates");
            Console.WriteLine("  the lower Truncated Normal Cumulative Density Function.");
            Console.WriteLine("");
            Console.WriteLine("        MU       S         A         X        CDF1           CDF2");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                TestValues.Truncated.truncated_normal_a_cdf_values(ref n_data, ref mu, ref sigma, ref a, ref x,
                    ref cdf1);

                if (n_data == 0)
                {
                    break;
                }

                cdf2 = CDF.truncated_normal_a_cdf(x, mu, sigma, a);

                Console.WriteLine("  " + mu.ToString().PadLeft(8)
                                       + "  " + sigma.ToString().PadLeft(8)
                                       + "  " + a.ToString().PadLeft(8)
                                       + "  " + x.ToString().PadLeft(8)
                                       + "  " + cdf1.ToString("0.################").PadLeft(24)
                                       + "  " + cdf2.ToString("0.################").PadLeft(24) + "");
            }
        }

        static void truncated_normal_a_cdf_inv_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_A_CDF_INV_TEST tests TRUNCATED_NORMAL_A_CDF_INV.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 March 2015
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
            int sample_num = 10;
            int seed;
            double sigma;
            double x;
            double x2;

            a = 50.0;
            mu = 100.0;
            sigma = 25.0;
            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_A_CDF_INV_TEST:");
            Console.WriteLine("  TRUNCATED_NORMAL_A_CDF_INV inverts the CDF of");
            Console.WriteLine("  the lower Truncated Normal Distribution.");
            Console.WriteLine("");
            Console.WriteLine("  The parent normal distribution has");
            Console.WriteLine("    mean =               " + mu + "");
            Console.WriteLine("    standard deviation = " + sigma + "");
            Console.WriteLine("  The parent distribution is truncated to");
            Console.WriteLine("  the interval [" + a + ",+oo)");
            Console.WriteLine("");
            Console.WriteLine("       X            CDF           CDF_INV");
            Console.WriteLine("");

            for (i = 0; i < sample_num; i++)
            {
                x = Truncated.truncated_normal_a_sample(mu, sigma, a, ref seed);
                cdf = CDF.truncated_normal_a_cdf(x, mu, sigma, a);
                x2 = CDF.truncated_normal_a_cdf_inv(cdf, mu, sigma, a);
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                                       + "  " + x.ToString().PadLeft(14)
                                       + "  " + cdf.ToString().PadLeft(14)
                                       + "  " + x2.ToString().PadLeft(14) + "");
            }

        }

        static void truncated_normal_a_mean_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_A_MEAN_TEST tests TRUNCATED_NORMAL_A_MEAN.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 March 2015
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
            int sample_num = 1000;
            int seed;
            double sigma;
            double[] x;
            double xmax;
            double xmin;

            a = 50.0;
            mu = 100.0;
            sigma = 25.0;
            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_A_MEAN_TEST");
            Console.WriteLine("  TRUNCATED_NORMAL_A_MEAN computes the mean");
            Console.WriteLine("  of the lower Truncated Normal Distribution.");
            Console.WriteLine("");
            Console.WriteLine("  The parent normal distribution has");
            Console.WriteLine("    mean =               " + mu + "");
            Console.WriteLine("    standard deviation = " + sigma + "");
            Console.WriteLine("  The parent distribution is truncated to");
            Console.WriteLine("  the interval [" + a + ",+oo)");

            mean = Truncated.truncated_normal_a_mean(mu, sigma, a);

            Console.WriteLine("");
            Console.WriteLine("  PDF mean = " + mean + "");

            x = new double[sample_num];

            for (i = 0; i < sample_num; i++)
            {
                x[i] = Truncated.truncated_normal_a_sample(mu, sigma, a, ref seed);
            }

            mean = typeMethods.r8vec_mean(sample_num, x);
            xmax = typeMethods.r8vec_max(sample_num, x);
            xmin = typeMethods.r8vec_min(sample_num, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + sample_num + "");
            Console.WriteLine("  Sample mean =     " + mean + "");
            Console.WriteLine("  Sample maximum =  " + xmax + "");
            Console.WriteLine("  Sample minimum =  " + xmin + "");


        }

        static void truncated_normal_a_moment_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_A_MOMENT_TEST tests TRUNCATED_NORMAL_A_MOMENT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 September 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            double[] a_test =
            {
                0.0, -10.0, 10.0, -10.0, 10.0, -10.0
            };
            double moment;
            double mu;
            double[] mu_test =
            {
                0.0, 0.0, 0.0, 0.0, 0.0, -5.0
            };
            int order;
            double sigma;
            double[] sigma_test =
            {
                1.0, 1.0, 1.0, 2.0, 2.0, 1.0
            };
            int test;
            int test_num;

            test_num = 6;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_A_MOMENT_TEST");
            Console.WriteLine("  TRUNCATED_NORMAL_A_MOMENT evaluates the moments");
            Console.WriteLine("  of the Lower Truncated Normal Distribution.");

            for (test = 0; test < test_num; test++)
            {
                mu = mu_test[test];
                sigma = sigma_test[test];
                a = a_test[test];
                Console.WriteLine("");
                Console.WriteLine("  Test = " + test
                                              + ", Mu = " + mu
                                              + ", Sigma = " + sigma
                                              + ", A = " + a + "");
                Console.WriteLine(" Order  Moment");
                Console.WriteLine("");
                for (order = 0; order <= 8; order++)
                {
                    moment = Truncated.truncated_normal_a_moment(order, mu, sigma, a);
                    Console.WriteLine("  " + order.ToString().PadLeft(2)
                                           + "  " + moment.ToString().PadLeft(14) + "");
                }
            }
        }

        static void truncated_normal_a_pdf_test()

            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    TRUNCATED_NORMAL_A_PDF_TEST tests TRUNCATED_NORMAL_A_PDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 September 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double mu = 0;
            int n_data;
            double pdf1 = 0;
            double pdf2;
            double sigma = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_A_PDF_TEST:");
            Console.WriteLine("  TRUNCATED_NORMAL_A_PDF evaluates the PDF of");
            Console.WriteLine("  the lower Truncated Normal Distribution.");
            Console.WriteLine("");
            Console.WriteLine("        MU       S         A         X        PDF1        PDF2");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                TestValues.Truncated.truncated_normal_a_pdf_values(ref n_data, ref mu, ref sigma, ref a, ref x,
                    ref pdf1);

                if (n_data == 0)
                {
                    break;
                }

                pdf2 = Truncated.truncated_normal_a_pdf(x, mu, sigma, a);

                Console.WriteLine("  " + mu.ToString().PadLeft(8)
                                       + "  " + sigma.ToString().PadLeft(8)
                                       + "  " + a.ToString().PadLeft(8)
                                       + "  " + x.ToString().PadLeft(8)
                                       + "  " + pdf1.ToString("0.################").PadLeft(24)
                                       + "  " + pdf2.ToString("0.################").PadLeft(24) + "");
            }
        }

        static void truncated_normal_a_sample_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_A_SAMPLE_TEST tests TRUNCATED_NORMAL_A_SAMPLE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            int i;
            double mu;
            int sample_num = 10;
            int seed;
            double sigma;
            double x;

            a = 50.0;
            mu = 100.0;
            sigma = 25.0;
            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_A_SAMPLE_TEST:");
            Console.WriteLine("  TRUNCATED_NORMAL_A_SAMPLE samples");
            Console.WriteLine("  the lower Truncated Normal Distribution.");
            Console.WriteLine("");
            Console.WriteLine("  The parent normal distribution has");
            Console.WriteLine("    mean =               " + mu + "");
            Console.WriteLine("    standard deviation = " + sigma + "");
            Console.WriteLine("  The parent distribution is truncated to");
            Console.WriteLine("  the interval [" + a + ",+oo)");
            Console.WriteLine("");

            for (i = 0; i < sample_num; i++)
            {
                x = Truncated.truncated_normal_a_sample(mu, sigma, a, ref seed);
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                                       + "  " + x.ToString().PadLeft(14) + "");
            }

        }

        static void truncated_normal_a_variance_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_A_VARIANCE_TEST tests TRUNCATED_NORMAL_A_VARIANCE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            int i;
            double mu;
            int sample_num = 1000;
            int seed;
            double sigma;
            double variance;
            double[] x;

            a = 50.0;
            mu = 100.0;
            sigma = 25.0;
            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_A_VARIANCE_TEST");
            Console.WriteLine("  TRUNCATED_NORMAL_A_VARIANCE computes the variance");
            Console.WriteLine("  of the lower Truncated Normal Distribution.");
            Console.WriteLine("");
            Console.WriteLine("  The parent normal distribution has");
            Console.WriteLine("    mean =               " + mu + "");
            Console.WriteLine("    standard deviation = " + sigma + "");
            Console.WriteLine("  The parent distribution is truncated to");
            Console.WriteLine("  the interval [" + a + ",+oo)");

            variance = Truncated.truncated_normal_a_variance(mu, sigma, a);

            Console.WriteLine("");
            Console.WriteLine("  PDF variance = " + variance + "");

            x = new double[sample_num];

            for (i = 0; i < sample_num; i++)
            {
                x[i] = Truncated.truncated_normal_a_sample(mu, sigma, a, ref seed);
            }

            variance = typeMethods.r8vec_variance(sample_num, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + sample_num + "");
            Console.WriteLine("  Sample variance = " + variance + "");


        }

        static void truncated_normal_ab_cdf_test()

            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    TRUNCATED_NORMAL_AB_CDF_TEST tests TRUNCATED_NORMAL_AB_CDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 September 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double b = 0;
            double cdf1 = 0;
            double cdf2 = 0;
            double mu = 0;
            int n_data;
            double sigma = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_AB_CDF_TEST:");
            Console.WriteLine("  TRUNCATED_NORMAL_AB_CDF evaluates");
            Console.WriteLine("  the Truncated Normal Cumulative Density Function.");
            Console.WriteLine("");
            Console.WriteLine("        MU       S         A         B         X        CDF1           CDF2");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                TestValues.Truncated.truncated_normal_ab_cdf_values(ref n_data, ref mu, ref sigma, ref a, ref b, ref x,
                    ref cdf1);

                if (n_data == 0)
                {
                    break;
                }

                cdf2 = CDF.truncated_normal_ab_cdf(x, mu, sigma, a, b);

                Console.WriteLine("  " + mu.ToString().PadLeft(8)
                                       + "  " + sigma.ToString().PadLeft(8)
                                       + "  " + a.ToString().PadLeft(8)
                                       + "  " + b.ToString().PadLeft(8)
                                       + "  " + x.ToString().PadLeft(8)
                                       + "  " + cdf1.ToString("0.################").PadLeft(24)
                                       + "  " + cdf2.ToString("0.################").PadLeft(24) + "");
            }
        }

        static void truncated_normal_ab_cdf_inv_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_AB_CDF_INV_TEST tests TRUNCATED_NORMAL_AB_CDF_INV.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 March 2015
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
            int sample_num = 10;
            int seed;
            double sigma;
            double x;
            double x2;

            a = 50.0;
            b = 150.0;
            mu = 100.0;
            sigma = 25.0;
            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_AB_CDF_INV_TEST:");
            Console.WriteLine("  TRUNCATED_NORMAL_AB_CDF_INV inverts the CDF of");
            Console.WriteLine("  the Truncated Normal Distribution.");
            Console.WriteLine("");
            Console.WriteLine("  The parent normal distribution has");
            Console.WriteLine("    mean =               " + mu + "");
            Console.WriteLine("    standard deviation = " + sigma + "");
            Console.WriteLine("  The parent distribution is truncated to");
            Console.WriteLine("  the interval [" + a + "," + b + "]");
            Console.WriteLine("");
            Console.WriteLine("       X            CDF           CDF_INV");
            Console.WriteLine("");

            for (i = 0; i < sample_num; i++)
            {
                x = Truncated.truncated_normal_ab_sample(mu, sigma, a, b, ref seed);
                cdf = CDF.truncated_normal_ab_cdf(x, mu, sigma, a, b);
                x2 = CDF.truncated_normal_ab_cdf_inv(cdf, mu, sigma, a, b);
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                                       + "  " + x.ToString().PadLeft(14)
                                       + "  " + cdf.ToString().PadLeft(14)
                                       + "  " + x2.ToString().PadLeft(14) + "");
            }

        }

        static void truncated_normal_ab_mean_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_AB_MEAN_TEST tests TRUNCATED_NORMAL_AB_MEAN.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 March 2015
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
            int sample_num = 1000;
            int seed;
            double sigma;
            double[] x;
            double xmax;
            double xmin;

            a = 50.0;
            b = 150.0;
            mu = 100.0;
            sigma = 25.0;
            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_AB_MEAN_TEST");
            Console.WriteLine("  TRUNCATED_NORMAL_AB_MEAN computes the mean");
            Console.WriteLine("  of the Truncated Normal Distribution.");
            Console.WriteLine("");
            Console.WriteLine("  The parent normal distribution has");
            Console.WriteLine("    mean =               " + mu + "");
            Console.WriteLine("    standard deviation = " + sigma + "");
            Console.WriteLine("  The parent distribution is truncated to");
            Console.WriteLine("  the interval [" + a + "," + b + "]");

            mean = Truncated.truncated_normal_ab_mean(mu, sigma, a, b);

            Console.WriteLine("");
            Console.WriteLine("  PDF mean = " + mean + "");

            x = new double[sample_num];

            for (i = 0; i < sample_num; i++)
            {
                x[i] = Truncated.truncated_normal_ab_sample(mu, sigma, a, b, ref seed);
            }

            mean = typeMethods.r8vec_mean(sample_num, x);
            xmax = typeMethods.r8vec_max(sample_num, x);
            xmin = typeMethods.r8vec_min(sample_num, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + sample_num + "");
            Console.WriteLine("  Sample mean =     " + mean + "");
            Console.WriteLine("  Sample maximum =  " + xmax + "");
            Console.WriteLine("  Sample minimum =  " + xmin + "");


        }

        static void truncated_normal_ab_moment_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_AB_MOMENT_TEST tests TRUNCATED_NORMAL_AB_MOMENT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 September 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            double[] a_test =
            {
                -1.0, 0.0, -1.0, -1.0, 0.0, 0.5, -2.0, -4.0, 4.0
            };
            double b;
            double[] b_test =
            {
                1.0, 1.0, 0.0, 1.0, 2.0, 2.0, 2.0, 4.0, 7.0
            };
            double moment;
            double mu;
            double[] mu_test =
            {
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 5.0
            };
            int order;
            double sigma;
            double[] sigma_test =
            {
                1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 0.5
            };
            int test;
            int test_num;

            test_num = 9;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_AB_MOMENT_TEST");
            Console.WriteLine("  TRUNCATED_NORMAL_AB_MOMENT evaluates the moments");
            Console.WriteLine("  of the Truncated Normal PDF:");

            for (test = 0; test < test_num; test++)
            {
                mu = mu_test[test];
                sigma = sigma_test[test];
                a = a_test[test];
                b = b_test[test];
                Console.WriteLine("");
                Console.WriteLine("  Test = " + test
                                              + ", Mu = " + mu
                                              + ", Sigma = " + sigma
                                              + ", A = " + a
                                              + ", B = " + b + "");
                Console.WriteLine(" Order  Moment");
                Console.WriteLine("");
                for (order = 0; order <= 8; order++)
                {
                    moment = Truncated.truncated_normal_ab_moment(order, mu, sigma, a, b);
                    Console.WriteLine("  " + order.ToString().PadLeft(2)
                                           + "  " + moment.ToString().PadLeft(14) + "");
                }
            }
        }

        static void truncated_normal_ab_pdf_test()

            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    TRUNCATED_NORMAL_AB_PDF_TEST tests TRUNCATED_NORMAL_AB_PDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 September 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double b = 0;
            double mu = 0;
            int n_data;
            double pdf1 = 0;
            double pdf2;
            double sigma = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_AB_PDF_TEST:");
            Console.WriteLine("  TRUNCATED_NORMAL_AB_PDF evaluates");
            Console.WriteLine("  the Truncated Normal Probability Density Function.");
            Console.WriteLine("");
            Console.WriteLine("        MU       S         A         B         X        PDF1        PDF2");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                TestValues.Truncated.truncated_normal_ab_pdf_values(ref n_data, ref mu, ref sigma, ref a, ref b, ref x,
                    ref pdf1);

                if (n_data == 0)
                {
                    break;
                }

                pdf2 = Truncated.truncated_normal_ab_pdf(x, mu, sigma, a, b);

                Console.WriteLine("  " + mu.ToString().PadLeft(8)
                                       + "  " + sigma.ToString().PadLeft(8)
                                       + "  " + a.ToString().PadLeft(8)
                                       + "  " + b.ToString().PadLeft(8)
                                       + "  " + x.ToString().PadLeft(8)
                                       + "  " + pdf1.ToString("0.################").PadLeft(24)
                                       + "  " + pdf2.ToString("0.################").PadLeft(24) + "");
            }
        }

        static void truncated_normal_ab_sample_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_AB_SAMPLE_TEST tests TRUNCATED_NORMAL_AB_SAMPLE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            double b;
            int i;
            double mu;
            int sample_num = 10;
            int seed;
            double sigma;
            double x;

            a = 50.0;
            b = 150.0;
            mu = 100.0;
            sigma = 25.0;
            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_AB_SAMPLE_TEST:");
            Console.WriteLine("  TRUNCATED_NORMAL_AB_SAMPLE samples");
            Console.WriteLine("  the Truncated Normal Distribution.");
            Console.WriteLine("");
            Console.WriteLine("  The parent normal distribution has");
            Console.WriteLine("    mean =               " + mu + "");
            Console.WriteLine("    standard deviation = " + sigma + "");
            Console.WriteLine("  The parent distribution is truncated to");
            Console.WriteLine("  the interval [" + a + "," + b + "]");
            Console.WriteLine("");

            for (i = 0; i < sample_num; i++)
            {
                x = Truncated.truncated_normal_ab_sample(mu, sigma, a, b, ref seed);
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                                       + "  " + x.ToString().PadLeft(14) + "");
            }

        }

        static void truncated_normal_ab_variance_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_AB_VARIANCE_TEST tests TRUNCATED_NORMAL_AB_VARIANCE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            double b;
            int i;
            double mu;
            int sample_num = 1000;
            int seed;
            double sigma;
            double variance;
            double[] x;

            a = 50.0;
            b = 150.0;
            mu = 100.0;
            sigma = 25.0;
            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_AB_VARIANCE_TEST");
            Console.WriteLine("  TRUNCATED_NORMAL_AB_VARIANCE computes the variance");
            Console.WriteLine("  of the Truncated Normal Distribution.");
            Console.WriteLine("");
            Console.WriteLine("  The parent normal distribution has");
            Console.WriteLine("    mean =               " + mu + "");
            Console.WriteLine("    standard deviation = " + sigma + "");
            Console.WriteLine("  The parent distribution is truncated to");
            Console.WriteLine("  the interval [" + a + "," + b + "]");

            variance = Truncated.truncated_normal_ab_variance(mu, sigma, a, b);

            Console.WriteLine("");
            Console.WriteLine("  PDF variance = " + variance + "");

            x = new double[sample_num];

            for (i = 0; i < sample_num; i++)
            {
                x[i] = Truncated.truncated_normal_ab_sample(mu, sigma, a, b, ref seed);
            }

            variance = typeMethods.r8vec_variance(sample_num, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + sample_num + "");
            Console.WriteLine("  Sample variance = " + variance + "");


        }

        static void truncated_normal_b_cdf_test()

            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    TRUNCATED_NORMAL_B_CDF_TEST tests TRUNCATED_NORMAL_B_CDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 September 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double b = 0;
            double cdf1 = 0;
            double cdf2 = 0;
            double mu = 0;
            int n_data;
            double sigma = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_B_CDF_TEST:");
            Console.WriteLine("  TRUNCATED_NORMAL_B_CDF evaluates");
            Console.WriteLine("  the upper Truncated Normal Cumulative Density Function.");
            Console.WriteLine("");
            Console.WriteLine("        MU       S         B         X        CDF1           CDF2");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                TestValues.Truncated.truncated_normal_b_cdf_values(ref n_data, ref mu, ref sigma, ref b, ref x,
                    ref cdf1);

                if (n_data == 0)
                {
                    break;
                }

                cdf2 = CDF.truncated_normal_b_cdf(x, mu, sigma, b);

                Console.WriteLine("  " + mu.ToString().PadLeft(8)
                                       + "  " + sigma.ToString().PadLeft(8)
                                       + "  " + b.ToString().PadLeft(8)
                                       + "  " + x.ToString().PadLeft(8)
                                       + "  " + cdf1.ToString("0.################").PadLeft(24)
                                       + "  " + cdf2.ToString("0.################").PadLeft(24) + "");
            }
        }

        static void truncated_normal_b_cdf_inv_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_B_CDF_INV_TEST tests TRUNCATED_NORMAL_B_CDF_INV.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 March 2015
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
            int sample_num = 10;
            int seed;
            double sigma;
            double x;
            double x2;

            b = 150.0;
            mu = 100.0;
            sigma = 25.0;
            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_B_CDF_INV_TEST:");
            Console.WriteLine("  TRUNCATED_NORMAL_B_CDF_INV inverts the CDF of");
            Console.WriteLine("  the upper Truncated Normal Distribution.");
            Console.WriteLine("");
            Console.WriteLine("  The parent normal distribution has");
            Console.WriteLine("    mean =               " + mu + "");
            Console.WriteLine("    standard deviation = " + sigma + "");
            Console.WriteLine("  The parent distribution is truncated to");
            Console.WriteLine("  the interval (-oo," + b + "]");
            Console.WriteLine("");
            Console.WriteLine("       X            CDF           CDF_INV");
            Console.WriteLine("");

            for (i = 0; i < sample_num; i++)
            {
                x = Truncated.truncated_normal_b_sample(mu, sigma, b, ref seed);
                cdf = CDF.truncated_normal_b_cdf(x, mu, sigma, b);
                x2 = CDF.truncated_normal_b_cdf_inv(cdf, mu, sigma, b);
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                                       + "  " + x.ToString().PadLeft(14)
                                       + "  " + cdf.ToString().PadLeft(14)
                                       + "  " + x2.ToString().PadLeft(14) + "");
            }

        }

        static void truncated_normal_b_mean_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_B_MEAN_TEST tests TRUNCATED_NORMAL_B_MEAN.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 March 2015
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
            int sample_num = 1000;
            int seed;
            double sigma;
            double[] x;
            double xmax;
            double xmin;

            b = 150.0;
            mu = 100.0;
            sigma = 25.0;
            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_B_MEAN_TEST");
            Console.WriteLine("  TRUNCATED_NORMAL_B_MEAN computes the mean");
            Console.WriteLine("  of the upper Truncated Normal Distribution.");
            Console.WriteLine("");
            Console.WriteLine("  The parent normal distribution has");
            Console.WriteLine("    mean =               " + mu + "");
            Console.WriteLine("    standard deviation = " + sigma + "");
            Console.WriteLine("  The parent distribution is truncated to");
            Console.WriteLine("  the interval (-oo," + b + "]");

            mean = Truncated.truncated_normal_b_mean(mu, sigma, b);

            Console.WriteLine("");
            Console.WriteLine("  PDF mean = " + mean + "");

            x = new double[sample_num];

            for (i = 0; i < sample_num; i++)
            {
                x[i] = Truncated.truncated_normal_b_sample(mu, sigma, b, ref seed);
            }

            mean = typeMethods.r8vec_mean(sample_num, x);
            xmax = typeMethods.r8vec_max(sample_num, x);
            xmin = typeMethods.r8vec_min(sample_num, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + sample_num + "");
            Console.WriteLine("  Sample mean =     " + mean + "");
            Console.WriteLine("  Sample maximum =  " + xmax + "");
            Console.WriteLine("  Sample minimum =  " + xmin + "");


        }

        static void truncated_normal_b_moment_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_B_MOMENT_TEST tests TRUNCATED_NORMAL_B_MOMENT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 September 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double b;
            double[] b_test =
            {
                0.0, 10.0, -10.0, 10.0, -10.0, 10.0
            };
            double moment;
            double mu;
            double[] mu_test =
            {
                0.0, 0.0, 0.0, 0.0, 0.0, 5.0
            };
            int order;
            double sigma;
            double[] sigma_test =
            {
                1.0, 1.0, 1.0, 2.0, 2.0, 1.0
            };
            int test;
            int test_num;

            test_num = 6;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_B_MOMENT_TEST");
            Console.WriteLine("  For the Upper Truncated Normal PDF:");
            Console.WriteLine("  TRUNCATED_NORMAL_B_MOMENT evaluates the moments.");

            for (test = 0; test < test_num; test++)
            {
                mu = mu_test[test];
                sigma = sigma_test[test];
                b = b_test[test];
                Console.WriteLine("");
                Console.WriteLine("  Test = " + test
                                              + ", Mu = " + mu
                                              + ", Sigma = " + sigma
                                              + ", B = " + b + "");
                Console.WriteLine(" Order  Moment");
                Console.WriteLine("");
                for (order = 0; order <= 8; order++)
                {
                    moment = Truncated.truncated_normal_b_moment(order, mu, sigma, b);
                    Console.WriteLine("  " + order.ToString().PadLeft(2)
                                           + "  " + moment.ToString().PadLeft(14) + "");
                }
            }
        }

        static void truncated_normal_b_pdf_test()

            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    TRUNCATED_NORMAL_B_PDF_TEST tests TRUNCATED_NORMAL_B_PDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 September 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double b = 0;
            double mu = 0;
            int n_data;
            double pdf1 = 0;
            double pdf2;
            double sigma = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_B_PDF_TEST:");
            Console.WriteLine("  TRUNCATED_NORMAL_B_PDF evaluates");
            Console.WriteLine("  the upper Truncated Normal Distribution.");
            Console.WriteLine("");
            Console.WriteLine("        MU       S         B         X        PDF1        PDF2");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                TestValues.Truncated.truncated_normal_b_pdf_values(ref n_data, ref mu, ref sigma, ref b, ref x,
                    ref pdf1);

                if (n_data == 0)
                {
                    break;
                }

                pdf2 = Truncated.truncated_normal_b_pdf(x, mu, sigma, b);

                Console.WriteLine("  " + mu.ToString().PadLeft(8)
                                       + "  " + sigma.ToString().PadLeft(8)
                                       + "  " + b.ToString().PadLeft(8)
                                       + "  " + x.ToString().PadLeft(8)
                                       + "  " + pdf1.ToString("0.################").PadLeft(24)
                                       + "  " + pdf2.ToString("0.################").PadLeft(24) + "");
            }
        }

        static void truncated_normal_b_sample_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_B_SAMPLE_TEST tests TRUNCATED_NORMAL_B_SAMPLE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double b;
            int i;
            double mu;
            int sample_num = 10;
            int seed;
            double sigma;
            double x;

            b = 150.0;
            mu = 100.0;
            sigma = 25.0;
            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_B_SAMPLE_TEST:");
            Console.WriteLine("  TRUNCATED_NORMAL_B_SAMPLE samples");
            Console.WriteLine("  the upper Truncated Normal Distribution.");
            Console.WriteLine("");
            Console.WriteLine("  The parent normal distribution has");
            Console.WriteLine("    mean =               " + mu + "");
            Console.WriteLine("    standard deviation = " + sigma + "");
            Console.WriteLine("  The parent distribution is truncated to");
            Console.WriteLine("  the interval (-oo," + b + "]");
            Console.WriteLine("");

            for (i = 0; i < sample_num; i++)
            {
                x = Truncated.truncated_normal_b_sample(mu, sigma, b, ref seed);
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                                       + "  " + x.ToString().PadLeft(14) + "");
            }

        }

        static void truncated_normal_b_variance_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_B_VARIANCE_TEST tests TRUNCATED_NORMAL_B_VARIANCE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double b;
            int i;
            double mu;
            int sample_num = 1000;
            int seed;
            double sigma;
            double variance;
            double[] x;

            b = 150.0;
            mu = 100.0;
            sigma = 25.0;
            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_B_VARIANCE_TEST");
            Console.WriteLine("  TRUNCATED_NORMAL_B_VARIANCE computes the variance");
            Console.WriteLine("  of the upper Truncated Normal Distribution.");
            Console.WriteLine("");
            Console.WriteLine("  The parent normal distribution has");
            Console.WriteLine("    mean =               " + mu + "");
            Console.WriteLine("    standard deviation = " + sigma + "");
            Console.WriteLine("  The parent distribution is truncated to");
            Console.WriteLine("  the interval (-oo," + b + "]");

            variance = Truncated.truncated_normal_b_variance(mu, sigma, b);

            Console.WriteLine("");
            Console.WriteLine("  PDF variance = " + variance + "");

            x = new double[sample_num];

            for (i = 0; i < sample_num; i++)
            {
                x[i] = Truncated.truncated_normal_b_sample(mu, sigma, b, ref seed);
            }

            variance = typeMethods.r8vec_variance(sample_num, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + sample_num + "");
            Console.WriteLine("  Sample variance = " + variance + "");
        }
    }
}