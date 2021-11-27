using System;
using System.Globalization;
using Burkardt.AppliedStatistics;
using Burkardt.Types;
using Burkardt.Uniform;

namespace ASA266Test;

internal static class Program
{
    private static void Main()
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA266_TEST.
//
//  Discussion:
//
//    ASA266_TEST tests the ASA266 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 June 2013
//
//  Author:
//
//    John Burkardt
//
    {
        Console.WriteLine("");
        Console.WriteLine("ASA266_TEST:");
        Console.WriteLine("  Test the ASA266 library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();
        test07();
        test08();
        test085();
        test09();
        test10();

        Console.WriteLine("");
        Console.WriteLine("ASA266_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");

    }

    private static void test01()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests ALNORM, NORMP, NPROB.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double ccdf2 = 0;
        double ccdf3 = 0;
        double cdf2 = 0;
        double cdf3 = 0;
        const int ntest = 16;
        double pdf2 = 0;
        double pdf3 = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  ALNORM,");
        Console.WriteLine("  NORMP, and");
        Console.WriteLine("  NPROB are routines that compute the cumulative");
        Console.WriteLine("  density function for the normal distribution.");
        Console.WriteLine("");
        Console.WriteLine("  X  CDF1  1-CDF1");
        Console.WriteLine("     CDF2  1-CDF2  PDF2");
        Console.WriteLine("     CDF3  1-CDF3  PDF3");

        for (int i = 0; i < ntest; i++)
        {
            double x = 3.0 * i / (ntest - 1);

            bool upper = false;
            double cdf1 = Algorithms.alnorm(x, upper);

            upper = true;
            double ccdf1 = Algorithms.alnorm(x, upper);

            Algorithms.normp(x, ref cdf2, ref ccdf2, ref pdf2);

            Algorithms.nprob(x, ref cdf3, ref ccdf3, ref pdf3);

            Console.WriteLine("");
            Console.WriteLine(x.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + cdf1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + ccdf1.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            Console.WriteLine("              "
                              + cdf2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + ccdf2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + pdf2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            Console.WriteLine("              "
                              + cdf3.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + ccdf3.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + pdf3.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void test02()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests PPND, PPND16.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int ntest = 9;

        int ifault = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  PPND,");
        Console.WriteLine("  PPND16 compute the percentage ");
        Console.WriteLine("  points of the normal distribution.");
        Console.WriteLine("");
        Console.WriteLine("           CDF     PPND(CDF)   PPND16(CDF)");
        Console.WriteLine("");

        for (i = 1; i <= ntest; i++)
        {
            double cdf = i / (double) (ntest + 1);
            double x1 = Algorithms.ppnd(cdf, ref ifault);
            double x2 = Algorithms.ppnd16(cdf, ref ifault);
            Console.WriteLine(cdf.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + x1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

    }

    private static void test03()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests DIGAMMA, R8_PSI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int ntest = 10;
        int ifault = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  digamma(X) = d ( Log ( Gamma ( X ) ) ) / dX.");
        Console.WriteLine("");
        Console.WriteLine("  DIGAMMA and");
        Console.WriteLine("  R8_PSI compute the digamma function:");
        Console.WriteLine("");
        Console.WriteLine("             X       DIGAMMA        R8_PSI");
        Console.WriteLine("");

        for (int i = 1; i <= ntest; i++)
        {
            double x = i / (double) ntest;

            Console.WriteLine(x.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + Algorithms.digamma(x, ref ifault).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + typeMethods.r8_psi(x).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void test04()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests TRIGAMMA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int ifault = 0;
        const int ntest = 10;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  TRIGAMMA computes the trigamma function:");
        Console.WriteLine("    trigamma(X) = d^2 ( Log ( Gamma ( X ) ) ) / dX^2.");
        Console.WriteLine("");
        Console.WriteLine("             X       TRIGAMMA");
        Console.WriteLine("");

        for (int i = 1; i <= ntest; i++)
        {
            double x = i / (double) ntest;

            double t = Algorithms.trigamma(x, ref ifault);
            Console.WriteLine(x.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + t.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void test05()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests ALNGAM, ALOGAM, R8_GAMMA_LOG, LNGAMMA;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int ntest = 10;

        int ifault = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  ALNGAM");
        Console.WriteLine("  ALOGAM,");
        Console.WriteLine("  R8_GAMMA_LOG, and");
        Console.WriteLine("  LNGAMMA compute the logarithm of the gamma function.");
        Console.WriteLine("");
        Console.WriteLine("             X        ALNGAM        ALOGAM    R8_GAMMA_LOG     LNGAMMA");
        Console.WriteLine("");

        for (int i = 1; i <= ntest; i++)
        {
            double x = i / (double) ntest;

            double log1 = Algorithms.alngam(x, ref ifault);
            double log2 = Algorithms.alogam(x, ref ifault);
            double log3 = typeMethods.r8_gamma_log(x);
            double log4 = Algorithms.lngamma(x, ref ifault);

            Console.WriteLine(x.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + log1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + log2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + log3.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + log4.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void test06()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests GAMAIN, GAMMDS, GAMMAD.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int ntest = 10;

        int ifault = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  GAMAIN, ");
        Console.WriteLine("  GAMMDS and ");
        Console.WriteLine("  GAMMAD compute the incomplete Gamma integral.");
        Console.WriteLine("");
        Console.WriteLine("             X             P        GAMMDS        GAMMAD        GAMAIN");
        Console.WriteLine("");

        for (int i = 1; i <= ntest; i++)
        {
            double x = i / (double) ntest;

            Console.WriteLine("");
            for (int j = 1; j <= ntest; j++)
            {
                double p = j / (double) ntest;
                double g1 = Algorithms.gammds(x, p, ref ifault);
                if (ifault != 0)
                {
                    g1 = -99.0;
                }

                double g2 = Algorithms.gammad(x, p, ref ifault);
                if (ifault != 0)
                {
                    g2 = -99.0;
                }

                double g3 = Algorithms.gamain(x, p, ref ifault);
                if (ifault != 0)
                {
                    g3 = -99.0;
                }

                Console.WriteLine(x.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                  + p.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                  + g1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                  + g2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                  + g3.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
    }

    private static void test07()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests PPCHI2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int nitest = 9;
        const int njtest = 9;

        int ifault = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  PPCHI2 computes the percentage points");
        Console.WriteLine("  of the chi squared distribution.");
        Console.WriteLine("");
        Console.WriteLine("      CDF      PPCHI2(CDF)");
        Console.WriteLine("");

        for (int j = 1; j <= njtest; j++)
        {
            double v = j;

            Console.WriteLine("");
            Console.WriteLine("  For Chi^2 parameter value = " + v + "");
            Console.WriteLine("");

            for (int i = 1; i <= nitest; i++)
            {
                double cdf = i / (double) (nitest + 1);
                double gg = Algorithms.alngam(v / 2.0, ref ifault);
                double x1 = Algorithms.ppchi2(cdf, v, gg, ref ifault);
                Console.WriteLine(cdf.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                  + x1.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
    }

    private static void test08()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests DIRICHLET_ESTIMATE, DIRICHLET_MEAN, DIRICHLET_VARIANCE.
        //
        //  Discussion:
        //
        //    Canned data is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int ELEM_NUM = 3;
        const int SAMPLE_NUM = 23;

        double aminus;
        double aplus;
        double eps = 0;
        int ifault = 0;
        int niter = 0;
        double rlogl = 0;
        double s = 0;
        double vari;
        double[] x =  {
            0.178, 0.162, 0.083, 0.087, 0.078, 0.040, 0.049, 0.100, 0.075, 0.084,
            0.060, 0.089, 0.050, 0.073, 0.064, 0.085, 0.094, 0.014, 0.060, 0.031,
            0.025, 0.045, 0.0195,
            0.346, 0.307, 0.448, 0.474, 0.503, 0.456, 0.363, 0.317, 0.394, 0.445,
            0.435, 0.418, 0.485, 0.378, 0.562, 0.465, 0.388, 0.449, 0.544, 0.569,
            0.491, 0.613, 0.526,
            0.476, 0.531, 0.469, 0.439, 0.419, 0.504, 0.588, 0.583, 0.531, 0.471,
            0.505, 0.493, 0.465, 0.549, 0.374, 0.450, 0.518, 0.537, 0.396, 0.400,
            0.484, 0.342, 0.4545
        };

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  For samples of a Dirichlet PDF,");
        Console.WriteLine("  DIRICHLET_ESTIMATE estimates the parameters.");
        Console.WriteLine("  DIRICHLET_MEAN finds the means;");
        Console.WriteLine("  DIRICHLET_VARIANCE finds the variances;");

        typeMethods.r8mat_print(SAMPLE_NUM, ELEM_NUM, x, "  Sampled data:");
        //
        //  Compute the observed averages.
        //
        double[] mean = typeMethods.r8col_mean(SAMPLE_NUM, ELEM_NUM, x);

        double[] variance = typeMethods.r8col_variance(SAMPLE_NUM, ELEM_NUM, x);

        Console.WriteLine("");
        Console.WriteLine("  Observed means, variances are:");
        Console.WriteLine("");
        for (int elem_i = 0; elem_i < ELEM_NUM; elem_i++)
        {
            Console.WriteLine(elem_i.ToString().PadLeft(6)
                              + mean[elem_i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + variance[elem_i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        const int init = 1;
        double[] alpha = new double[ELEM_NUM];
        double[] g = new double[ELEM_NUM];
        double[] v = new double[ELEM_NUM * ELEM_NUM];

        Algorithms.dirichlet_estimate(ELEM_NUM, SAMPLE_NUM, x, SAMPLE_NUM,
            init, ref alpha, ref rlogl, ref v, ref g, ref niter, ref s, ref eps, ref ifault);

        if (ifault != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("WARNING!");
            Console.WriteLine("  DIRICHLET_ESTIMATE error code:");
            Console.WriteLine("  IFAULT = " + ifault + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Index, Estimate, Lower Limit, Upper Limit:");
        Console.WriteLine("");

        for (int elem_i = 0; elem_i < ELEM_NUM; elem_i++)
        {
            vari = v[elem_i + elem_i * ELEM_NUM];
            aminus = alpha[elem_i] - 1.96 * Math.Sqrt(vari);
            aplus = alpha[elem_i] + 1.96 * Math.Sqrt(vari);
            Console.WriteLine(elem_i.ToString().PadLeft(6)
                              + alpha[elem_i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + aminus.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + aplus.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        mean = Algorithms.dirichlet_mean(ELEM_NUM, alpha);

        variance = Algorithms.dirichlet_variance(ELEM_NUM, alpha);

        Console.WriteLine("");
        Console.WriteLine("  Expected means, variances are:");
        Console.WriteLine("");
        for (int elem_i = 0; elem_i < ELEM_NUM; elem_i++)
        {
            Console.WriteLine(elem_i.ToString().PadLeft(6)
                              + mean[elem_i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + variance[elem_i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        double alpha_sum = typeMethods.r8vec_sum(ELEM_NUM, alpha);

        Console.WriteLine("");
        Console.WriteLine("  Alpha sum is " + alpha_sum + "");
        Console.WriteLine("");
        Console.WriteLine("  NORMALIZED VALUES:");
        Console.WriteLine("  Index, Estimate, Lower Limit, Upper Limit:");
        Console.WriteLine("");

        for (int elem_i = 0; elem_i < ELEM_NUM; elem_i++)
        {
            vari = v[elem_i + elem_i * ELEM_NUM];
            aminus = (alpha[elem_i] - 1.96 * Math.Sqrt(vari)) / alpha_sum;
            aplus = (alpha[elem_i] + 1.96 * Math.Sqrt(vari)) / alpha_sum;
            Console.WriteLine(elem_i.ToString().PadLeft(6)
                              + (alpha[elem_i] / alpha_sum).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + aminus.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + aplus.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Log likelihood function = " + rlogl + "");

    }

    private static void test085()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST085 tests GAMMA_SAMPLE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int rep_num = 5;
        const int test_num = 5;

        typeMethods.r8NormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST085");
        Console.WriteLine("  GAMMA_SAMPLE samples a Gamma distribution.");

        int seed = 123456789;

        for (int test = 1; test <= test_num; test++)
        {
            double a = UniformRNG.r8_uniform_ab(0.1, 2.0, ref seed);
            double b = UniformRNG.r8_uniform_ab(0.1, 2.0, ref seed);
            Console.WriteLine("");
            Console.WriteLine("  A = " + a + ", B = " + b + "");
            for (int rep = 1; rep <= rep_num; rep++)
            {
                double x = Algorithms.gamma_sample(a, b, ref data, ref seed);
                Console.WriteLine("  " + rep.ToString().PadLeft(2)
                                       + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
    }

    private static void test09()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests DIRICHLET_ESTIMATE, _MEAN, _VARIANCE, _SAMPLE.
        //
        //  Discussion:
        //
        //    Data is generated by sampling a distribution with known parameters.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int ELEM_NUM = 3;
        const int SAMPLE_NUM = 1000;

        double[] alpha =  {
                3.22, 20.38, 21.68
            }
            ;
        double aminus;
        double aplus;
        double eps = 0;
        int ifault = 0;
        int niter = 0;
        double rlogl = 0;
        double s = 0;
        double vari;

        int seed = 123456789;
        typeMethods.r8NormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  For a Dirichlet distribution,");
        Console.WriteLine("  DIRICHLET_SAMPLE samples;");
        Console.WriteLine("  DIRICHLET_MEAN finds the means;");
        Console.WriteLine("  DIRICHLET_VARIANCE finds the variances;");
        Console.WriteLine("  DIRICHLET_ESTIMATE estimates the parameters.");
        //
        //  Report.
        //
        typeMethods.r8vec_print(ELEM_NUM, alpha, "  Distribution parameters:");

        double[] mean = Algorithms.dirichlet_mean(ELEM_NUM, alpha);

        double[] variance = Algorithms.dirichlet_variance(ELEM_NUM, alpha);

        Console.WriteLine("");
        Console.WriteLine("  Distribution means, variances are:");
        Console.WriteLine("");
        for (int elem_i = 0; elem_i < ELEM_NUM; elem_i++)
        {
            Console.WriteLine(elem_i.ToString().PadLeft(6)
                              + mean[elem_i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + variance[elem_i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  Sample the distribution.
        //
        double[] x_sample = new double[SAMPLE_NUM * ELEM_NUM];

        Console.WriteLine("");
        Console.WriteLine("  Number of samples is " + SAMPLE_NUM + "");

        for (int sample_i = 0; sample_i < SAMPLE_NUM; sample_i++)
        {
            double[] x = Algorithms.dirichlet_sample(ELEM_NUM, alpha, ref data, ref seed);

            for (int elem_i = 0; elem_i < ELEM_NUM; elem_i++)
            {
                x_sample[sample_i + elem_i * SAMPLE_NUM] = x[elem_i];
            }
        }

        //
        //  Print some results.
        //
        Console.WriteLine("");
        Console.WriteLine("  First few samples:");
        Console.WriteLine("");

        for (int sample_i = 0; sample_i < Math.Min(SAMPLE_NUM, 10); sample_i++)
        {
            string cout = sample_i.ToString().PadLeft(6);
            for (int elem_i = 0; elem_i < ELEM_NUM; elem_i++)
            {
                cout += x_sample[sample_i + elem_i * SAMPLE_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }

        //
        //  Compute means, variances.
        //
        mean = typeMethods.r8col_mean(SAMPLE_NUM, ELEM_NUM, x_sample);

        variance = typeMethods.r8col_variance(SAMPLE_NUM, ELEM_NUM, x_sample);

        Console.WriteLine("");
        Console.WriteLine("  Observed means, variances are:");
        Console.WriteLine("");
        for (int elem_i = 0; elem_i < ELEM_NUM; elem_i++)
        {
            Console.WriteLine(elem_i.ToString().PadLeft(6)
                              + mean[elem_i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + variance[elem_i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  Destroy the values of ALPHA.
        //
        for (int elem_i = 0; elem_i < ELEM_NUM; elem_i++)
        {
            alpha[elem_i] = 0.0;
        }

        //
        //  Try to recover the values of ALPHA.
        //
        const int init = 1;
        double[] v = new double[ELEM_NUM * ELEM_NUM];
        double[] g = new double[ELEM_NUM];

        Algorithms.dirichlet_estimate(ELEM_NUM, SAMPLE_NUM, x_sample, SAMPLE_NUM,
            init, ref alpha, ref rlogl, ref v, ref g, ref niter, ref s, ref eps, ref ifault);

        if (ifault != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("Warning!");
            Console.WriteLine("  DIRICHLET_ESTIMATE error code:");
            Console.WriteLine("  IFAULT = " + ifault + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Index, Estimate, Lower Limit, Upper Limit:");
        Console.WriteLine("");

        for (int elem_i = 0; elem_i < ELEM_NUM; elem_i++)
        {
            vari = v[elem_i + elem_i * ELEM_NUM];
            aminus = alpha[elem_i] - 1.96 * Math.Sqrt(vari);
            aplus = alpha[elem_i] + 1.96 * Math.Sqrt(vari);
            Console.WriteLine(elem_i.ToString().PadLeft(6)
                              + alpha[elem_i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + aminus.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + aplus.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        double alpha_sum = typeMethods.r8vec_sum(ELEM_NUM, alpha);

        Console.WriteLine("");
        Console.WriteLine("  Alpha sum is " + alpha_sum + "");
        Console.WriteLine("");
        Console.WriteLine("  NORMALIZED VALUES:");
        Console.WriteLine("  Index, Estimate, Lower Limit, Upper Limit:");
        Console.WriteLine("");

        for (int elem_i = 0; elem_i < ELEM_NUM; elem_i++)
        {
            vari = v[elem_i + elem_i * ELEM_NUM];
            aminus = (alpha[elem_i] - 1.96 * Math.Sqrt(vari)) / alpha_sum;
            aplus = (alpha[elem_i] + 1.96 * Math.Sqrt(vari)) / alpha_sum;
            Console.WriteLine(elem_i.ToString().PadLeft(6)
                              + (alpha[elem_i] / alpha_sum).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + aminus.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                              + aplus.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Log likelikhood function = " + rlogl + "");

    }

    private static void test10()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tests DIRICHLET_MIX_SAMPLE, _MEAN, _VARIANCE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int COMP_NUM = 3;
        const int COMP_MAX = 3;
        const int ELEM_NUM = 3;
        const int SAMPLE_NUM = 200;

        double[] a = new double[ELEM_NUM];
        double[] alpha =  {
                0.05, 0.85, 0.00,
                0.20, 0.10, 0.50,
                0.75, 0.05, 0.50
            }
            ;
        double[] comp_weight = {
                3.0, 2.0, 1.0
            }
            ;
        double[] mean;
        double[] variance;

        int seed = 123456789;

        typeMethods.r8NormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  For a Dirichlet mixture distribution,");
        Console.WriteLine("  DIRICHLET_MIX_SAMPLE samples;");
        Console.WriteLine("  DIRICHLET_MIX_MEAN computes means;");
        Console.WriteLine("  DIRICHLET_MIX_VARIANCE computes variances.");
        //
        //  Report.
        //
        typeMethods.r8vec_print(COMP_NUM, comp_weight, "  Component weight:");

        Console.WriteLine("");
        Console.WriteLine("  Component  Parameters Means Variances");
        for (int comp_i = 0; comp_i < COMP_NUM; comp_i++)
        {
            Console.WriteLine("");
            string cout = comp_i.ToString().PadLeft(6);
            for (int elem_i = 0; elem_i < ELEM_NUM; elem_i++)
            {
                a[elem_i] = alpha[comp_i + elem_i * COMP_MAX];
            }

            mean = Algorithms.dirichlet_mean(ELEM_NUM, a);
            variance = Algorithms.dirichlet_variance(ELEM_NUM, a);
            for (int elem_i = 0; elem_i < ELEM_NUM; elem_i++)
            {
                Console.WriteLine(cout + elem_i.ToString().PadLeft(6)
                                       + "  " + alpha[comp_i + elem_i * COMP_MAX].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + mean[elem_i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + variance[elem_i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            }

        }

        mean = Algorithms.dirichlet_mix_mean(COMP_MAX, COMP_NUM, ELEM_NUM, alpha,
            comp_weight);

        typeMethods.r8vec_print(ELEM_NUM, mean, "  Element means:");
        //
        //  Sample the distribution.
        //
        int[] comp_sample = new int[SAMPLE_NUM];
        double[] x_sample = new double[ELEM_NUM * SAMPLE_NUM];
        Console.WriteLine("");
        Console.WriteLine("  Number of samples is " + SAMPLE_NUM + "");

        for (int sample_i = 0; sample_i < SAMPLE_NUM; sample_i++)
        {
            int comp_i = 0;
            double[] x = Algorithms.dirichlet_mix_sample(COMP_MAX, COMP_NUM, ELEM_NUM, alpha,
                comp_weight, ref data, ref seed, ref comp_i);

            for (int elem_i = 0; elem_i < ELEM_NUM; elem_i++)
            {
                x_sample[elem_i + sample_i * ELEM_NUM] = x[elem_i];
            }

            comp_sample[sample_i] = comp_i;

        }

        //
        //  Print some results.
        //
        Console.WriteLine("");
        Console.WriteLine("  First few samples:");
        Console.WriteLine("");
        Console.WriteLine("  Sample  Component  X");
        Console.WriteLine("");

        for (int sample_i = 0; sample_i < Math.Min(SAMPLE_NUM, 10); sample_i++)
        {
            string cout = "  " + sample_i.ToString().PadLeft(2)
                               + "  " + comp_sample[sample_i].ToString().PadLeft(2);
            for (int elem_i = 0; elem_i < ELEM_NUM; elem_i++)
            {
                cout += "  " + x_sample[elem_i + sample_i * ELEM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(10);
            }

            Console.WriteLine(cout);
        }
        //
        //  Compute the observed averages.
        //
        mean = typeMethods.r8col_mean(SAMPLE_NUM, ELEM_NUM, x_sample);

        variance = typeMethods.r8col_variance(SAMPLE_NUM, ELEM_NUM, x_sample);

        Console.WriteLine("");
        Console.WriteLine("  Element  Observed mean, variance");
        Console.WriteLine("");
        for (int elem_i = 0; elem_i < ELEM_NUM; elem_i++)
        {
            Console.WriteLine(elem_i.ToString().PadLeft(6)
                              + "  " + mean[elem_i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                              + "  " + variance[elem_i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

    }


}