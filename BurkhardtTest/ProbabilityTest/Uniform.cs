﻿using System;
using System.Globalization;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void uniform_01_order_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_01_ORDER_SAMPLE_TEST tests UNIFORM_01_ORDER_SAMPLE;
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
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("UNIFORM_01_ORDER_SAMPLE_TEST");
        Console.WriteLine("  For the Uniform 01 Order PDF:");
        Console.WriteLine("  UNIFORM_ORDER_SAMPLE samples.");

        int n = 10;
        double[] x = Burkardt.Probability.Uniform.uniform_01_order_sample(n, ref seed);

        typeMethods.r8vec_print(n, x, "  Ordered sample:");
    }

    private static void uniform_nsphere_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_NSPHERE_SAMPLE_TEST tests UNIFORM_NSPHERE_SAMPLE;
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
        int seed = 123456789;

        const int n = 3;

        Console.WriteLine("");
        Console.WriteLine("UNIFORM_NSPHERE_SAMPLE_TEST");
        Console.WriteLine("  UNIFORM_NSPHERE_SAMPLE samples the Uniform Nsphere distribution.");

        Console.WriteLine("");
        Console.WriteLine("  Dimension N of sphere = " + n + "");
        Console.WriteLine("");
        Console.WriteLine("  Points on the sphere:");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double[] x = Burkardt.Probability.Uniform.uniform_nsphere_sample(n, ref seed);
            string cout = "  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  ";
            int j;
            for (j = 0; j < n; j++)
            {
                cout += x[j].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
            }

            Console.WriteLine(cout);

        }

    }

    private static void uniform_01_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_01_CDF_TEST tests UNIFORM_01_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("UNIFORM_01_CDF_TEST");

        Console.WriteLine("  UNIFORM_01_CDF evaluates the Uniform 01 CDF;");
        Console.WriteLine("  UNIFORM_01_CDF_INV inverts the Uniform 01 CDF.");
        Console.WriteLine("  UNIFORM_01_PDF evaluates the Uniform 01 PDF;");

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Burkardt.Probability.Uniform.uniform_01_sample(ref seed);
            double pdf = Burkardt.Probability.Uniform.uniform_01_pdf(x);
            double cdf = Burkardt.Probability.Uniform.uniform_01_cdf(x);
            double x2 = Burkardt.Probability.Uniform.uniform_01_cdf_inv(cdf);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    private static void uniform_01_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_01_SAMPLE_TEST tests UNIFORM_01_SAMPLE.
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
        const int SAMPLE_NUM = 1000;

        int i;
        int seed = 123456789;
        double[] x = new double [SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("UNIFORM_01_SAMPLE_TEST");
        Console.WriteLine("  UNIFORM_01_MEAN computes the Uniform 01 mean;");
        Console.WriteLine("  UNIFORM_01_SAMPLE samples the Uniform 01 distribution;");
        Console.WriteLine("  UNIFORM_01_VARIANCE computes the Uniform 01 variance.");

        double mean = Burkardt.Probability.Uniform.uniform_01_mean();
        double variance = Burkardt.Probability.Uniform.uniform_01_variance();

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Burkardt.Probability.Uniform.uniform_01_sample(ref seed);
        }

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

    private static void uniform_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_CDF_TEST tests UNIFORM_CDF.
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
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("UNIFORM_CDF_TEST");
        Console.WriteLine("  UNIFORM_CDF evaluates the Uniform CDF;");
        Console.WriteLine("  UNIFORM_CDF_INV inverts the Uniform CDF.");
        Console.WriteLine("  UNIFORM_PDF evaluates the Uniform PDF;");

        double a = 1.0;
        double b = 10.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Burkardt.Probability.Uniform.uniform_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("UNIFORM_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Burkardt.Probability.Uniform.uniform_sample(a, b, ref seed);
            double pdf = Burkardt.Probability.Uniform.uniform_pdf(x, a, b);
            double cdf = Burkardt.Probability.Uniform.uniform_cdf(x, a, b);
            double x2 = Burkardt.Probability.Uniform.uniform_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void uniform_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_SAMPLE_TEST tests UNIFORM_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        const int SAMPLE_NUM = 1000;

        int i;
        int seed = 123456789;
        double[] x = new double [SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("UNIFORM_SAMPLE_TEST");
        Console.WriteLine("  UNIFORM_MEAN computes the Uniform mean;");
        Console.WriteLine("  UNIFORM_SAMPLE samples the Uniform distribution;");
        Console.WriteLine("  UNIFORM_VARIANCE computes the Uniform variance;");

        double a = 1.0;
        double b = 10.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Burkardt.Probability.Uniform.uniform_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("UNIFORM_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Burkardt.Probability.Uniform.uniform_mean(a, b);
        double variance = Burkardt.Probability.Uniform.uniform_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Burkardt.Probability.Uniform.uniform_sample(a, b, ref seed);
        }

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

    private static void uniform_discrete_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_DISCRETE_CDF_TEST tests UNIFORM_DISCRETE_CDF.
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
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("UNIFORM_DISCRETE_CDF_TEST");
        Console.WriteLine("  UNIFORM_DISCRETE_CDF evaluates the Uniform Discrete CDF;");
        Console.WriteLine("  UNIFORM_DISCRETE_CDF_INV inverts the Uniform Discrete CDF.");
        Console.WriteLine("  UNIFORM_DISCRETE_PDF evaluates the Uniform Discrete PDF;");

        int a = 1;
        int b = 6;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Burkardt.Probability.Uniform.uniform_discrete_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("UNIFORM_DISCRETE_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            int x = Burkardt.Probability.Uniform.uniform_discrete_sample(a, b, ref seed);
            double pdf = Burkardt.Probability.Uniform.uniform_discrete_pdf(x, a, b);
            double cdf = Burkardt.Probability.Uniform.uniform_discrete_cdf(x, a, b);
            int x2 = Burkardt.Probability.Uniform.uniform_discrete_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void uniform_discrete_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_DISCRETE_SAMPLE_TEST tests UNIFORM_DISCRETE_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        const int SAMPLE_NUM = 1000;

        int i;
        int seed = 123456789;
        int[] x = new int[SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("UNIFORM_DISCRETE_SAMPLE_TEST");
        Console.WriteLine("  UNIFORM_DISCRETE_MEAN computes the Uniform Discrete mean;");
        Console.WriteLine("  UNIFORM_DISCRETE_SAMPLE samples the Uniform Discrete distribution;");
        Console.WriteLine("  UNIFORM_DISCRETE_VARIANCE computes the Uniform Discrete variance;");

        int a = 1;
        int b = 6;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Burkardt.Probability.Uniform.uniform_discrete_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("UNIFORM_DISCRETE_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Burkardt.Probability.Uniform.uniform_discrete_mean(a, b);
        double variance = Burkardt.Probability.Uniform.uniform_discrete_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Burkardt.Probability.Uniform.uniform_discrete_sample(a, b, ref seed);
        }

        mean = typeMethods.i4vec_mean(SAMPLE_NUM, x);
        variance = typeMethods.i4vec_variance(SAMPLE_NUM, x);
        int xmax = typeMethods.i4vec_max(SAMPLE_NUM, x);
        int xmin = typeMethods.i4vec_min(SAMPLE_NUM, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
        Console.WriteLine("  Sample mean =     " + mean + "");
        Console.WriteLine("  Sample variance = " + variance + "");
        Console.WriteLine("  Sample maximum =  " + xmax + "");
        Console.WriteLine("  Sample minimum =  " + xmin + "");

    }

}